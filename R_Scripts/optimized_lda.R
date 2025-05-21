#' Optimized LDA Analyzer
#'
#' An R6 class for performing optimized Latent Dirichlet Allocation (LDA) analysis 
#' with GO term integration and network visualization capabilities.
#'
#' @field master_df The master data frame containing gene and GO term information
#' @field dtm Document-term matrix (internal use)
#' @field lda_model The fitted LDA model (internal use)
#' 
#' @importFrom R6 R6Class
#' @import tidyverse topicmodels ldatuning future furr ggraph igraph tm AnnotationDbi org.EcK12.eg.db widyr
#' @export
library(R6)
library(tidyverse)
library(topicmodels)
library(ldatuning)
library(future)
library(furrr)
library(ggraph)
library(igraph)
library(tm)
library(AnnotationDbi)
library(org.EcK12.eg.db)
library(widyr)

OptimizedLDAAnalyzer <- R6Class(
  "OptimizedLDAAnalyzer",
  public = list(
    combined_network = NULL,
    doc_doc_edges = NULL,
    doc_network = NULL,
    doc_network_edges = NULL,
    doc_network_nodes = NULL,
    doc_network_nodes_expanded = NULL,
    doc_nodes = NULL,
    doc_topics = NULL,
    doc_topic_edges = NULL,
    dtm = NULL,
    expanded_edges = NULL,
    g_edges = NULL,
    g_nodes = NULL,
    genotype_edges = NULL,
    go_corpus_col = NULL,
    go_map = NULL,
    k_metrics = NULL,
    lda_model = NULL,
    master_df = NULL,
    multiplex_edges = NULL,
    multiplex_network = NULL,
    multiplex_nodes = NULL,
    geno_graph = NULL,
    optimal_k = NULL,
    terms = NULL,
    term_term_edges = NULL,
    topic_network = NULL,
    top_terms = NULL,
    topic_edges = NULL,
    topic_doc_network = NULL,
    topic_gene_matrix = NULL,
    topic_labels = NULL,
    topic_sim = NULL,
    topic_terms = NULL,

    # Constructor with input validation
    #' @description 
    #' Initialize a new OptimizedLDAAnalyzer object
    #' 
    #' @param master_df Data frame containing required columns. The required columns are identity (a unique ideentifier in the form gene_strain where strain is 2 letters),
    #' genotype (an identifier for the experimental condition, two letter (i.e WT)), Gene (the gene name, i.e. cyoA), and GO (the GO term, i.e GO:0000123).
    #' @param go_corpus_col Name of column containing GO term strings
    #' 
    #' @examples
    #' \dontrun{
    #' analyzer <- OptimizedLDAAnalyzer$new(
    #'   master_df = my_data,
    #'   go_corpus_col = "GO_terms"
    #' )}
    initialize = function(master_df, go_corpus_col) {
      # Check for required columns
      required_cols <- c("identity", "genotype", "Gene", "GO")  # Add required columns
      missing_cols <- setdiff(required_cols, names(master_df))
      if (length(missing_cols) > 0) {
        stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
      }
      
      # Initialize all fields
      self$master_df <- master_df %>%
        mutate(identity = as.character(identity))
      self$go_corpus_col <- go_corpus_col
      private$.valid_go = NULL
      message("Initialization complete.")
    },
    
    validate = function() {
      private$validate_inputs()
      private$validate_go_terms()
      message("Validation complete.")
    },
    
    get_valid_go = function() {
      if (is.null(private$.valid_go)) {
        private$validate_go_terms()  # Auto-populate if missing
      }
      return(private$.valid_go)
    },
    
    set_valid_go = function(value) {
      stop("Use validate_go_terms() to set valid_go")
    },
    
    get_top_terms = function(n = NULL) {
      if (is.null(self$top_terms)) {
        warning("Top terms not extracted yet. Run fit_lda() first.")
        return(NULL)
      }
      if (!is.null(n)) {
        return(self$top_terms %>% group_by(topic) %>% slice_head(n = n))
      }
      return(self$top_terms)
    },
    
    get_formatted_top_terms = function(n = 1, wide_format = TRUE) {
      self$terms <- self$top_terms
      if (wide_format) {
        self$terms %>%
          group_by(topic) %>%
          summarize(
            term = paste(term, collapse = ", "),
            mean_beta = mean(beta)
          )
      } else {
        self$terms
      }
    },
    
    build_go_map = function() {
      if (is.null(self$terms)) {
        message("Formatted top terms is empty, creating tibble...")
        self$get_formatted_top_terms()
      }
      
      all_go_ids <- unique(unlist(strsplit(self$master_df$GO, split = ",\\s*")))
      all_go_ids <- all_go_ids[!is.na(all_go_ids)]
      
      self$go_map <- tryCatch({
        data.frame(lapply(self$terms, gsub, pattern = "go:", replacement = "GO:")) %>%
          mutate(
            TERM = AnnotationDbi::Term(term)
          ) %>% 
          group_by(topic) %>%
          arrange(rank, .by_group = TRUE) %>%
          slice_head(n = 1) %>%
          ungroup() %>%
          mutate(topic = paste0("Topic_", topic))
      }, error = function(e) {
        warning("GO term mapping failed: ", e$message)
        data.frame(GOID = all_go_ids, Term = all_go_ids)
      })
    },
    
    # 1. Optimized DTM Preparation
    #' Prepare Document-Term Matrix
    #' 
    #' Creates a sparse document-term matrix from GO term annotations
    #' 
    #' @param p_sparse Maximum sparsity threshold (0-1)
    #' @return A DocumentTermMatrix object
    #' @examples
    #' \dontrun{
    #' analyzer$prepare_dtm(p_sparse = 0.999)
    #' }
    prepare_dtm = function(p_sparse = 0.999) {
      valid_go <- self$get_valid_go()
      message("Building optimized DTM...")
      self$dtm <- self$master_df %>%
        mutate(across(self$go_corpus_col, list(all = all_of, na = ~ifelse(is.na(.), "", .)
                                              )
                      )
               ) %>%
        unnest_tokens(
          word, 
          !!sym(self$go_corpus_col), 
          token = "regex",
          pattern = ",\\s*"
        ) %>%
        #filter(word %in% valid_go$GO) %>%  # Keep only valid GO terms
        count(Gene, word) %>%
        cast_dtm(Gene, word, n) %>%
        removeSparseTerms(sparse = p_sparse) %>%
        .[rowSums(as.matrix(.)) > 0, ]
      
      if (nrow(self$dtm) == 0) stop("DTM has zero rows after filtering")
      invisible(self)
      return(self$dtm)
    },
    
    # 2. Parallelized Topic Optimization
    #' Find Optimal Topic Number
    #' 
    #' Uses parallel processing to determine the optimal number of topics
    #' 
    #' @param k_range Range of topic numbers to evaluate (e.g., 3:11)
    #' @param metrics Metrics to use for evaluation (default: CaoJuan2009, Griffiths2004)
    #' @return List containing metrics data and optimal k value
    #' @seealso \code{\link[ldatuning]{FindTopicsNumber}}
    find_optimal_k = function(k_range = 3:11, metrics = c("CaoJuan2009", "Griffiths2004")) {
      message("Tuning topics with parallel processing...")
      plan(multisession)
      
      self$k_metrics <- FindTopicsNumber(
        self$dtm,
        topics = k_range,
        metrics = metrics,
        method = "Gibbs",
        control = list(seed = 123, keep = 50),
        mc.cores = availableCores() - 1
      )
      
      self$optimal_k <- self$k_metrics %>%
        pivot_longer(-topics) %>%
        group_by(name) %>%
        mutate(scaled_value = scale(value)) %>%
        group_by(topics) %>%
        summarise(score = mean(scaled_value)) %>%
        filter(score == max(score)) %>%
        pull(topics)
      
      message("Optimal k: ", self$optimal_k)
      invisible(self)
      return(list(k_metrics = self$k_metrics, optimal_k = self$optimal_k))
    },
    
    # 3. Efficient LDA Fitting
    #' Fit a Latent Dirichlet Allocation (LDA) model
    #' 
    #' Uses the optimal k parameter to fit an LDA model to the DTM.
    #' 
    #' @param k A static k parameter if the find_k_optimal() function is not to be utilized. Null if autamated optimal k determination should be used.  (e.g., 3:11)
    #' @param methods The sampling method to use for LDA fitting (default: "Gibbs")
    #' @param iter Number of iterations for the given sampling method (default: 1000)
    #' @param n_top_terms Number of top terms to extract for each topic (default: 10)
    #' @return Returns an LDA model object
    #' @seealso \code{\link[lda_model]{}}
    fit_lda = function(k = NULL, method = "Gibbs", iter = 1000, n_top_terms = 10) {
      k <- ifelse(is.null(k), self$optimal_k, k)
      message("Fitting LDA with k=", k, "...")
      
      self$lda_model <- LDA(
        self$dtm,
        k = k,
        method = method,
        control = list(
          seed = 123,
          keep = 50,
          iter = iter,
          verbose = 1
        )
      )
      
      self$doc_topics <- tidy(self$lda_model, matrix = "gamma") %>%
        mutate(gamma = ifelse(is.nan(gamma), 0, gamma))
      # Store top terms
      private$extract_top_terms(n_top_terms)
      return(self$lda_model)
      invisible(self)
    },
    
    #' Build Genotype Network
    #' 
    #' Constructs a network graph connecting genes to topics based on their genotype associations.
    #' 
    #' @return A tbl_graph object representing the genotype-topic network
    #' @examples
    #' \dontrun{
    #' analyzer$build_genotype_network()
    #' }
    build_genotype_network = function() {
      self$doc_topics <- tidy(self$lda_model, matrix = "gamma")
      
      self$g_edges <- self$doc_topics %>%
        left_join(
          self$master_df %>% select(Gene, genotype, log2FoldChange, identity),
          by = c("document" = "Gene"),
          relationship = "many-to-many"
        )
      print(self$g_edges)
      self$g_nodes <- self$g_edges %>%
        distinct(document, genotype, log2FoldChange) %>%
        mutate(
          importance = abs(log2FoldChange),
          color = cut(
            log2FoldChange,
            breaks = c(-Inf, -1, 1, Inf),
            labels = c("blue", "gray", "red")
          )
        )
      
      self$geno_graph <- tbl_graph(
        nodes = self$g_nodes,
        edges = self$g_edges %>% select(from = document, to = topic, weight = gamma)
      )
      self$geno_graph <- self$geno_graph %>%
        mutate(
          degree = centrality_degree(),
          betweenness = centrality_betweenness(),
          #closeness = centrality_closeness(),
          page_rank = centrality_pagerank(),
          eigen = centrality_eigen(),
          centrality_alpha = centrality_alpha()
        )
      return(self$geno_graph)
    },
    
    #' Build Document-Document Network
    #' 
    #' Creates a network of document (gene) similarities based on their topic distributions
    #' 
    #' @param similarity_threshold Minimum similarity threshold for including edges (default: 0.1)
    #' @return A tbl_graph object representing the document-document network
    #' @examples
    #' \dontrun{
    #' analyzer$build_document_network(similarity_threshold = 0.2)
    #' }
    build_document_network = function(similarity_threshold = 0.1) {
      if (is.null(self$doc_topics)) {
        stop("Run fit_lda() first to generate document-topic distributions")
      }
      
      # Calculate document-document similarities based on topic distributions
      self$doc_doc_edges <- self$doc_topics %>%
        left_join(
          self$master_df %>% distinct(Gene, genotype),
          by = c("document" = "Gene")
        ) %>%
        mutate(full_id = paste(document, genotype),
               type = "doc_doc") %>%
        select(full_id, topic, gamma) %>%
        pairwise_similarity(full_id, topic, gamma, upper = FALSE) %>%
        filter(similarity > similarity_threshold) %>%
        rename(from = item1, to = item2, weight = similarity)
      
      # Create document nodes with metadata
      self$doc_nodes <- self$master_df %>%
        distinct(identity, genotype, log2FoldChange, timepoint, padj, sig_star) %>%
        mutate(full_id = paste(identity, genotype)) %>%
        rename(name = identity, group = genotype)
      
      # Create the document network graph
      self$doc_network <- tbl_graph(
        nodes = self$doc_nodes %>% mutate(type = "document"),
        edges = self$doc_doc_edges %>%
          mutate(type = "document_document")
      ) %>%
        mutate(
          centrality = centrality_degree()
          #community = group_infomap()
        )
      
      message(sprintf(
        "Built document network with:\n- %d documents\n- %d connections",
        igraph::vcount(self$doc_network),
        igraph::ecount(self$doc_network)
      ))
      
      return(self$doc_network)
    },
    
    # Modified: Enhanced topic network with genotype connections
    #' Build Topic Network
    #' 
    #' Constructs a network of topics and genotypes based on similarity measures.
    #' 
    #' @param threshold Minimum similarity threshold for including topic edges (default: 0.001)
    #' @param genotype_edge_weight Scaling factor for genotype edge weights (default: 0.5)
    #' @return A tbl_graph object representing the topic-genotype network
    #' @examples
    #' \dontrun{
    #' analyzer$build_topic_network(threshold = 0.01)
    #' }
    build_topic_network = function(threshold = 0.001, genotype_edge_weight = 0.5) {
      # Original topic-topic network
      self$topic_terms <- tidy(self$lda_model, matrix = "beta")
      self$topic_sim <- self$topic_terms %>%
        pairwise_similarity(topic, term, beta, upper = FALSE) %>%
        filter(similarity > threshold)
      
      # Genotype-genotype edges (using identity column)
      self$genotype_edges <- self$master_df %>%
        group_by(genotype) %>%
        summarise(genes = list(unique(identity))) %>%
        pairwise_count(genotype, genes) %>%
        dplyr::rename(from = item1, to = item2, weight = n) %>%
        mutate(weight = genotype_edge_weight * weight / max(weight))
      
      # Combined graph
      self$topic_network <- tbl_graph(
        nodes = tibble(
          name = c(
            paste0("Topic_", unique(c(self$topic_sim$item1, self$topic_sim$item2))),
            unique(c(self$genotype_edges$from, self$genotype_edges$to))
          ),
          group = c(
            rep("topic", length(unique(c(self$topic_sim$item1, self$topic_sim$item2)))),
            rep("genotype", length(unique(c(self$genotype_edges$from, self$genotype_edges$to))))
          ),
          type = c(
            rep("topic", length(unique(c(self$topic_sim$item1, self$topic_sim$item2)))),
            rep("genotype", length(unique(c(self$genotype_edges$from, self$genotype_edges$to))))
          ),
          funcs = ifelse(
            name %in% self$go_map$topic,
            self$go_map$TERM[match(name, self$go_map$topic)],
            NA_character_
          )
        ),
        edges = bind_rows(
          self$topic_sim %>% 
            dplyr::rename(from = item1, to = item2, weight = similarity) %>%
            mutate(
              from = paste0("Topic_", from),
              to = paste0("Topic_", to),
              type = "topic_topic"
            ),
          self$genotype_edges
        )
      ) %>%
        mutate(degree = centrality_degree(),
               betweenness = centrality_betweenness(),
               centrality_alpha = centrality_alpha())
      
      return(self$topic_network)
    },

    # New: Create topic-document graph
    #' Build Topic-Document Graph
    #' 
    #' Creates a bipartite graph connecting documents (genes) to topics with edges weighted by gamma values.
    #' 
    #' @param gamma_threshold Minimum gamma value for including an edge (default: 0.05)
    #' @return A tbl_graph object representing the topic-document relationships
    #' @examples
    #' \dontrun{
    #' analyzer$build_topic_document_graph(gamma_threshold = 0.1)
    #' }
    build_topic_document_network = function(gamma_threshold = 0.05) {
      if (is.null(self$doc_topics)) {
        stop("Run build_genotype_network() first")
      }
      
      if (is.null(self$go_map)) {
        message("GO map not found, creating...")
        self$build_go_map()
      }
      
      # 1. Create all possible document-topic-genotype combinations
      self$expanded_edges <- self$doc_topics %>%
        left_join(
          self$master_df %>% distinct(Gene, genotype, self$go_corpus_col),
          by = c("document" = "Gene")
        ) %>%
        mutate(
          full_id = paste(document, genotype),  # Create "cyoA WT" style IDs
          topic = paste0("Topic_", topic)
        ) %>%
        filter(gamma > gamma_threshold) %>%
        select(from = full_id, to = topic, weight = gamma)#, top_func = func)
      
      # 2. Prepare document nodes from master_df
      doc_nodes <- self$master_df %>%
        distinct(identity, genotype, log2FoldChange, timepoint, padj, sig_star) %>%
        dplyr::rename(name = identity, group = genotype)
      
      # 3. Build graph (corrected)
      self$topic_doc_network <- tbl_graph(
        nodes = bind_rows(
          doc_nodes %>% mutate(type = "document"),
          tibble(name = unique(self$expanded_edges$to), 
                 group = "topic", 
                 type = "topic",
                 # Add term mapping - this is the key fix
                 funcs = ifelse(
                   name %in% self$go_map$topic,
                   self$go_map$TERM[match(name, self$go_map$topic)],
                   NA_character_
                 )
          )
        ),
        edges = self$expanded_edges %>%
          mutate(type = "topic_document")
      ) %>%
        mutate(
          centrality = centrality_degree(),
          community = group_infomap()
        ) 
      
      message(sprintf(
        "Success! Built graph with:\n- %d documents\n- %d topics\n- %d connections",
        sum(igraph::V(self$topic_doc_network)$type == "document"),
        sum(igraph::V(self$topic_doc_network)$type == "topic"),
        igraph::ecount(self$topic_doc_network)
      ))
      return(self$topic_doc_network)
    },
    
    #' Combine Networks
    #' 
    #' Merges topic-document, document-document, and topic-term networks into a single graph
    #' 
    #' @return A combined tbl_graph object
    #' @examples
    #' \dontrun{
    #' combined <- analyzer$combine_networks()
    #' }
    combine_networks = function() {
      # Ensure all required networks are built
      if (is.null(self$topic_doc_network)) {
        self$build_topic_document_network()
      }
      if (is.null(self$doc_network)) {
        self$build_document_network()
      }
      if (is.null(self$topic_network)) {
        self$build_topic_network()
      }
      
      # Get all graphs
      topic_doc_graph <- self$topic_doc_network
      doc_doc_graph <- self$doc_network
      topic_graph <- self$topic_network
      
      # Rename edges (to or to and from) for clarity.
      topic_doc_graph <- topic_doc_graph %>%
        activate(edges) %>%
        rename(to = to, from = from)
      
      # First combine document networks
      combined_docs <- graph_join(topic_doc_graph, doc_doc_graph)
      
      # Then combine with term network
      self$combined_network <- graph_join(combined_docs, topic_graph)
      
      # Add edge types
      self$combined_network <- self$combined_network %>%
        activate(edges) %>%
        mutate(
          edge_color = case_when(
            type == "topic_document" ~ "blue",
            type == "document_document" ~ "gray",
            type == "topic_topic" ~ "red",
            TRUE ~ "black"
          )
        ) %>%
        activate(nodes) %>%
        mutate(
          node_type = case_when(
            type == "document" ~ "document",
            type == "topic" ~ "topic",
            TRUE ~ "other"
          ),
          node_size = case_when(
            node_type == "document" ~ 2,
            node_type == "topic" ~ 5,
            node_type == "term" ~ 2,
            TRUE ~ 1
          )
        )
      
      message(sprintf(
        "Successfully combined networks with:\n- %d nodes\n- %d edges",
        igraph::vcount(self$combined_network),
        igraph::ecount(self$combined_network)
      ))
      
      return(self$combined_network)
    },
    
    # Visualization for topic-document graph
    #' Plot Topic-Document Graph
    #' 
    #' Visualizes the topic-document network with customizable labeling.
    #' 
    #' @param label_threshold Centrality threshold for node labeling (default: 0.5)
    #' @return A ggplot object showing the topic-document network
    #' @examples
    #' \dontrun{
    #' analyzer$plot_topic_document_graph(label_threshold = 0.7)
    #' }
    plot_topic_document_graph = function(label_threshold = 0.5) {
      if (is.null(self$topic_doc_network)) {
        stop("Run build_topic_document_graph() first")
      }
      
      ggraph(self$topic_doc_network %>% activate(edges) %>% filter(weight > 0), 
             layout = "kk", weights = weight) +
        geom_edge_link(aes(alpha = weight), width = 0.1, color = "gray50") +
        geom_node_point(aes(color = group, size = ifelse(type == "topic", 5, 3))) +
        geom_node_label(
          aes(label = ifelse(type == "topic" | centrality > label_threshold, funcs, "")), 
          repel = TRUE
        ) +
        theme_graph() +
        labs(title = "Topic-Document Network")
    },
    
    # Visualization for enhanced topic network
    #' Plot Enhanced Topic Network
    #' 
    #' Visualizes the combined topic and genotype network.
    #' 
    #' @return A ggplot object showing the enhanced topic network
    #' @examples
    #' \dontrun{
    #' analyzer$plot_enhanced_topic_network()
    #' }
    plot_topic_graph = function() {
      if (is.null(self$topic_network)) {
        stop("Run build_topic_network() first")
      }
      
      ggraph(self$topic_network, layout = "fr") +
        geom_edge_link(aes(alpha = weight, color = ifelse(is.na(similarity), "genotype", "topic"))) +
        geom_node_point(aes(size = degree, color = type)) +
        geom_node_label(aes(label = ifelse(degree > median(degree), funcs, "")), repel = TRUE) +
        scale_color_manual(values = c("topic" = "red", "genotype" = "blue")) +
        theme_graph()
    },
    
    plot_multiplex_graph = function(gamma_threshold = 0.1, similarity_threshold = 0.01) {
      # Prepare both edge types
      self$doc_topic_edges <- self$doc_topics %>%
        filter(gamma > gamma_threshold) %>%
        mutate(from = paste(document, genotype),  # Match your composite IDs
               to = paste0("Topic_", topic),
               type = "doc_topic")
      
      self$term_term_edges <- self$topic_sim %>%
        filter(similarity > similarity_threshold) %>%
        mutate(from = paste0("TERM_", item1),
               to = paste0("TERM_", item2),
               type = "term_term")
      
      self$multiplex_nodes <- bind_rows(
        self$master_df %>%
          distinct(identity, genotype) %>%
          rename(name = identity, group = genotype) %>%
          mutate(type = "document"),
        
        tibble(name = paste0("Topic_", unique(self$doc_topic_edges$topic)), 
               group = "topic", 
               type = "topic"),
        
        tibble(name = unique(c(self$term_term_edges$from,
                               self$term_term_edges$to)),
               group = "term", 
               type = "term")
      )
      
      self$multiplex_edges <- bind_rows(
        self$doc_topic_edges %>% select(from, to, weight = gamma, type),
        self$term_term_edges %>% select(from, to, weight = similarity, type)
      )
      
      # Create unified graph
      self$multiplex_network <- tbl_graph(
        nodes = self$multiplex_nodes,
        edges = self$multiplex_edges,
        directed = FALSE
      )
    
      
      # Plot with distinct edge styles
       # multiplex_graph <- ggraph(self$multiplex_network, layout = "fr") +
       #   geom_edge_link(
       #     aes(filter = type == "doc_topic", 
       #         alpha = weight, 
       #         color = "Document-Topic"),
       #     width = 0.8
       #   ) +
       #   geom_edge_link(
       #     aes(filter = type == "term_term", 
       #         alpha = weight,
       #         color = "Term-Term"),
       #     width = 0.5,
       #     linetype = "dashed"
       #   ) +
       #   geom_node_point(
       #     aes(color = group, size = ifelse(type == "term", 2, 4))
       #   ) +
       #   scale_edge_color_manual(
       #     name = "Edge Type",
       #     values = c("Document-Topic" = "darkblue", "Term-Term" = "firebrick")
       #   ) +
       #   theme_graph() +
       #   labs(title = "Multiplex Network: Documents, Topics, and Terms")
      
      return(multiplex_network)
    }
  ),
    
  private = list(
    .valid_go = NULL,
    
    # Properly named and formatted private methods
    validate_inputs = function() {
      if (!self$go_corpus_col %in% names(self$master_df)) {
        stop("Column '", self$go_corpus_col, "' not found")
      }
      if (all(is.na(self$master_df[[self$go_corpus_col]]))) {
        stop("GO term column contains only NA values")
      }
    },
    
    validate_go_terms = function() {
      all_go_ids <- unique(na.omit(unlist(
        strsplit(self$master_df[[self$go_corpus_col]], ",\\s*")
      )))
      
      if (!requireNamespace("org.EcK12.eg.db", quietly = TRUE)) {
        stop("Required package org.EcK12.eg.db is not installed")
      }
      
      result <- tryCatch({
        res <- AnnotationDbi::select(
          org.EcK12.eg.db,
          keys = all_go_ids,
          columns = c("GO", "TERM"),
          keytype = "GO"
        )
        if (nrow(res) == 0) {
          data.frame(GO = all_go_ids, TERM = all_go_ids)
        } else {
          res %>% distinct(GO, .keep_all = TRUE)
        }
      }, error = function(e) {
        warning("GO validation failed: ", e$message)
        data.frame(GO = all_go_ids, TERM = all_go_ids)
      })
      
      # Ensure we always set .valid_go, even if it's just the raw IDs
      private$.valid_go <- result
      invisible(self)
    },
    
    extract_top_terms = function(n_top_terms) {
      self$top_terms <- tidy(self$lda_model, matrix = "beta")
      
      # Add GO term definitions if available
      val_go <- self$get_valid_go()  # Changed variable name and access method
      if (!is.null(val_go)) {
        self$top_terms <- self$top_terms %>%
          left_join(val_go, by = c("term" = "GO")) %>%
          mutate(term = ifelse(is.na(TERM), term, TERM))
      }
      
      # Rank terms within topics
      self$top_terms <- self$top_terms %>%
        group_by(topic) %>%
        arrange(desc(beta), .by_group = TRUE) %>%
        mutate(rank = row_number()) %>%
        ungroup()
      
      invisible(self)
    }
  )
)
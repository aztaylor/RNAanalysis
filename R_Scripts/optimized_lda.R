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
    doc_topics = NULL,
    dtm = NULL,
    expanded_edges = NULL,
    g_edges = NULL,
    g_nodes = NULL,
    go_corpus_col = NULL,
    go_map = NULL,
    k_metrics = NULL,
    lda_model = NULL,
    master_df = NULL,
    net_graph = NULL,
    optimal_k = NULL,
    terms = NULL,
    tn_graph = NULL,
    top_terms = NULL,
    topic_doc_graph = NULL,
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
      
      self$net_graph <- tbl_graph(
        nodes = self$g_nodes,
        edges = self$g_edges %>% select(from = document, to = topic, weight = gamma)
      )
      self$net_graph <- self$net_graph %>%
        mutate(
          degree = centrality_degree(),
          betweenness = centrality_betweenness(),
          #closeness = centrality_closeness(),
          page_rank = centrality_pagerank(),
          eigen = centrality_eigen(),
          centrality_alpha = centrality_alpha()
        )
      return(self$net_graph)
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
    build_topic_document_graph = function(gamma_threshold = 0.05) {
      if (is.null(self$doc_topics)) {
        stop("Run build_genotype_network() first")
      }
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
      self$topic_doc_graph <- tbl_graph(
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
        edges = self$expanded_edges
      ) %>%
        mutate(
          centrality = centrality_degree(),
          community = group_infomap()
        ) 
      
      message(sprintf(
        "Success! Built graph with:\n- %d documents\n- %d topics\n- %d connections",
        sum(igraph::V(self$topic_doc_graph)$type == "document"),
        sum(igraph::V(self$topic_doc_graph)$type == "topic"),
        igraph::ecount(self$topic_doc_graph)
      ))
      return(self$topic_doc_graph)
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
      genotype_edges <- self$master_df %>%
        group_by(genotype) %>%
        summarise(genes = list(unique(identity))) %>%
        pairwise_count(genotype, genes) %>%
        dplyr::rename(from = item1, to = item2, weight = n) %>%
        mutate(weight = genotype_edge_weight * weight / max(weight))
      
      # Combined graph
      self$tn_graph <- tbl_graph(
        nodes = tibble(
          name = c(
            paste0("Topic_", unique(c(self$topic_sim$item1, self$topic_sim$item2))),
            unique(c(genotype_edges$from, genotype_edges$to))
          ),
          type = c(
            rep("topic", length(unique(c(self$topic_sim$item1, self$topic_sim$item2)))),
            rep("genotype", length(unique(c(genotype_edges$from, genotype_edges$to))))
          )
        ),
        edges = bind_rows(
          self$topic_sim %>% 
            dplyr::rename(from = item1, to = item2, weight = similarity) %>%
            mutate(
              from = paste0("Topic_", from),
              to = paste0("Topic_", to)
            ),
          genotype_edges
        )
      ) %>%
        mutate(degree = centrality_degree(),
          betweenness = centrality_betweenness(),
          centrality_alpha = centrality_alpha())
      
      return(self$tn_graph)
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
      if (is.null(self$topic_doc_graph)) {
        stop("Run build_topic_document_graph() first")
      }
      
      ggraph(self$topic_doc_graph %>% activate(edges) %>% filter(weight > 0), 
             layout = "fr", weights = weight) +
        geom_edge_link(aes(alpha = weight, length = weight), width = 0.1, color = "gray50") +
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
    plot_enhanced_topic_network = function() {
      if (is.null(self$tn_graph)) {
        stop("Run build_topic_network() first")
      }
      
      ggraph(self$tn_graph, layout = "fr") +
        geom_edge_link(aes(alpha = weight, color = ifelse(is.na(similarity), "genotype", "topic"))) +
        geom_node_point(aes(size = degree, color = type)) +
        geom_node_label(aes(label = ifelse(degree > median(degree), name, "")), repel = TRUE) +
        scale_color_manual(values = c("topic" = "red", "genotype" = "blue")) +
        theme_graph()
    },
    
    #' Combine Networks
    #' 
    #' Merges topic-document and topic-term networks into a single graph.
    #' 
    #' @param analyzer Another OptimizedLDAAnalyzer instance to combine with
    #' @return A combined tbl_graph object
    #' @examples
    #' \dontrun{
    #' combined <- analyzer1$combine_networks(analyzer2)
    #' }
    combine_networks = function(analyzer) {
      # Get both graphs
      doc_graph <- analyzer$topic_doc_graph
      term_graph <- analyzer$tn_graph
      
      # Rename term graph nodes to avoid conflicts
      term_graph <- term_graph %>%
        activate(nodes) %>%
        mutate(name = paste0("TERM_", name))  # Prefix term names
      
      # Combine graphs
      combined_graph <- graph_join(doc_graph, term_graph)
      
      # Optional: Add edges between topics in both graphs
      topic_links <- analyzer$doc_topics %>%
        group_by(topic) %>%
        summarise(mean_gamma = mean(gamma)) %>%
        left_join(analyzer$topic_sim, by = c("topic" = "item1"))
      
      return(combined_graph)
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
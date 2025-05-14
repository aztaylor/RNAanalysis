library(R6)
library(tidytext)
library(topicmodels)
library(ldatuning)
library(future)
library(furrr)
library(org.EcK12.eg.db)
library(AnnotationDbi)

OptimizedLDAAnalyzer <- R6Class(
  "OptimizedLDAAnalyzer",
  public = list(
    # Fields
    master_df = NULL,
    go_corpus_col = NULL,
    dtm = NULL,
    lda_model = NULL,
    k_metrics = NULL,
    optimal_k = NULL,
    valid_go = NULL,
    
    # Constructor with input validation
    initialize = function(master_df, go_corpus_col) {
      self$master_df <- master_df
      self$go_corpus_col <- go_corpus_col
      private$validate_inputs()
      private$validate_go_terms()
      print(valid_go)
    },
    
    # 1. Optimized DTM Preparation
    prepare_dtm = function(p_sparse = 0.999) {
      message("Building optimized DTM...")
      self$dtm <- self$master_df %>%
        mutate(across(all_of(self$go_corpus_col), ~ifelse(is.na(.), "", .))) %>%
        unnest_tokens(
          word, 
          !!sym(self$go_corpus_col), 
          token = "regex", 
          pattern = ",\\s*"  # Handles comma-separated GO terms
        ) %>%
        #filter(word %in% self$valid_go$GO) %>%  # Keep only valid GO terms
        count(Gene, word) %>%
        cast_dtm(Gene, word, n) %>%
        removeSparseTerms(sparse = p_sparse) %>%
        .[rowSums(as.matrix(.)) > 0, ]
      
      if (nrow(self$dtm) == 0) stop("DTM has zero rows after filtering")
      invisible(self)
    },
    
    # 2. Parallelized Topic Optimization
    find_optimal_k = function(k_range = 3:11, metrics = c("CaoJuan2009", "Griffiths2004")) {
      message("Tuning topics with parallel processing...")
      plan(multisession)  # Enable parallel
      
      self$k_metrics <- FindTopicsNumber(
        self$dtm,
        topics = k_range,
        metrics = metrics,
        method = "Gibbs",
        control = list(seed = 123, keep = 50),  # keep=50 for Griffiths2004
        mc.cores = availableCores() - 1,

      )
      
      self$optimal_k <- self$k_metrics %>%
        pivot_longer(-topics) %>%
        group_by(topics) %>%
        summarise(score = mean(scale(value))) %>%
        filter(score == max(score)) %>%
        pull(topics)
      
      message("Optimal k: ", self$optimal_k)
      invisible(self)
    },
    
    # 3. Efficient LDA Fitting
    fit_lda = function(k = NULL, method = "Gibbs", iter = 1000) {
      k <- ifelse(is.null(k), self$optimal_k, k)
      message("Fitting LDA with k=", k, "...")
      
      self$lda_model <- LDA(
        self$dtm,
        k = 3,
        method = method,
        control = list(
          seed = 123,
          keep = 50,        # For log-likelihood tracking
          iter = iter,      # Increased iterations
          verbose = 1
        )
      )
      invisible(self)
    },
    
    # 4. Enhanced Visualization
    
    visualize_topics = function(n_topTerms = 10) {
      top_terms <- tidy(self$lda_model, matrix = "beta") %>%
        group_by(topic) %>%
        slice_max(beta, n = n_topTerms, with_ties = FALSE) %>%
        ungroup() %>%
        left_join(self$valid_go, by = c("term" = "GOID")) %>%
        mutate(
          term = ifelse(is.na(.data$term), term, .data$term)  # Safe column reference
        )
      
       ggplot(top_terms, aes(x = reorder(term, beta), y = beta, fill = factor(topic))) +
         geom_col(show.legend = FALSE) +
         facet_wrap(~topic, scales = "free") +
         coord_flip() +
         labs(title = paste("Top", n_topTerms, "Terms per Topic"))
    },
    
    # 5. Network Visualization (Optimized)
    plot_network = function(
    n_top_docs = 15,
    n_top_terms = 10,
    layout = "fr",
    edge_threshold = 0.1
    ) {
      # Extract document-topic relationships
      doc_topics <- tidy(self$lda_model, matrix = "gamma") %>%
        filter(gamma > edge_threshold)  # Threshold for meaningful connections
      
      # Create graph (adapted from your original network code)
      # ... [network creation logic] ...
      graph_data <- doc_topics %>%
        left_join(
          self$master_df %>% select(Gene, genotype, log2FoldChange),
          by = c("document" = "Gene"),
          relationship = "many-to-many"
        ) %>%
        group_by(topic) %>%
        slice_max(gamma, n = 10) %>%
        ungroup()
      
      # Return optimized ggraph object
      ggraph(graph_data, layout = layout) +
        geom_edge_arc(aes(alpha = weight), color = "grey70") +
        geom_node_point(aes(color = type, size = centrality)) +
        geom_node_text(aes(label = ifelse(centrality > 0.1, name, "")), repel = TRUE) +
        theme_graph()
    }
  ),
  
  private = list(
    # Input validation
    validate_inputs = function() {
      if (!self$go_corpus_col %in% names(self$master_df)) {
        stop("Column '", self$go_corpus_col, "' not found. Available: ", 
             paste(names(self$master_df), collapse = ", "))
      }
    },
    
    # GO Term validation (optimized)
    validate_go_terms = function() {
      # 1. Extract GO IDs safely
      all_go_ids <- unique(na.omit(unlist(
        strsplit(self$master_df[[self$go_corpus_col]], ",\\s*")
      )))
      
      # 2. Debug print (optional)
      message("Found ", length(all_go_ids), " unique GO IDs")
      
      # 3. Database lookup with proper error handling
      self$valid_go <- tryCatch({
        result <- AnnotationDbi::select(
          org.EcK12.eg.db,
          keys = all_go_ids,
          columns = c("GO", "TERM"),
          keytype = "GO"
        )
        
        # Ensure consistent column names
        if (nrow(result) > 0) {
          result %>% 
            distinct(GO, .keep_all = TRUE) %>%
            rename(GOID = GO)  # Standardize column name
        } else {
          data.frame(GOID = all_go_ids, TERM = all_go_ids)
        }
      }, error = function(e) {
        warning("GO validation failed: ", e$message)
        data.frame(GOID = all_go_ids, TERM = all_go_ids)
      })
      
      # Return self for chaining
      invisible(self)
    }
  )
)
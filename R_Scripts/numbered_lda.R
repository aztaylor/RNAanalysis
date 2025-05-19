1 library(R6)
2 library(tidyverse)
3 library(topicmodels)
4 library(ldatuning)
5 library(future)
6 library(furrr)
7 library(ggraph)
8 library(igraph)
9 library(tm)
10 library(AnnotationDbi)
11 library(org.EcK12.eg.db)
12 
13 OptimizedLDAAnalyzer <- R6Class(
14   "OptimizedLDAAnalyzer",
15   public = list(
16     # Existing fields
17     master_df = NULL,
18     go_corpus_col = NULL,
19     dtm = NULL,
20     lda_model = NULL,
21     k_metrics = NULL,
22     optimal_k = NULL,
23     net_graph = NULL,
24     doc_topics = NULL,
25     g_edges = NULL,
26     g_nodes = NULL,
27     tn_graph = NULL,
28     top_terms = NULL,
29     topic_terms = NULL,
30     topic_sim = NULL,
31     topic_labels = NULL,
32     topic_gene_matrix = NULL,
33     #topic_quality = NULL,
34     topic_doc_graph = NULL,
35     go_map = NULL,
36     terms = NULL,
37 
38     
39     # Constructor with input validation
40     initialize = function(master_df, go_corpus_col) {
41       # Check for required columns
42       required_cols <- c("identity", "genotype", "Gene", "GO")  # Add required columns
43       missing_cols <- setdiff(required_cols, names(master_df))
44       if (length(missing_cols) > 0) {
45         stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
46       }
47       # Initialize all fields
48       self$master_df <- master_df %>%
49         mutate(identity = as.character(identity))
50       self$go_corpus_col <- go_corpus_col
51       self$terms <- tibble()  # Initialize empty terms tibble
52       self$go_map <- tibble() # Initialize empty GO map
53       private$.valid_go = NULL
54       message("Initialization complete.")
55     },
56     
57     validate = function() {
58       private$validate_inputs()
59       private$validate_go_terms()
60       message("Validation complete.")
61     },
62     
63     get_valid_go = function() {
64       if (is.null(private$.valid_go)) {
65         private$validate_go_terms()  # Auto-populate if missing
66       }
67       return(private$.valid_go)
68     },
69     
70     set_valid_go = function(value) {
71       stop("Use validate_go_terms() to set valid_go")
72     },
73     
74     # 1. Optimized DTM Preparation
75     prepare_dtm = function(p_sparse = 0.999) {
76       valid_go <- self$get_valid_go()
77       message("Building optimized DTM...")
78       self$dtm <- self$master_df %>%
79         mutate(across(self$go_corpus_col, list(all = all_of, na = ~ifelse(is.na(.), "", .)
80                                               )
81                       )
82                ) %>%
83         unnest_tokens(
84           word, 
85           !!sym(self$go_corpus_col), 
86           token = "regex",
87           pattern = ",\\s*"
88         ) %>%
89         #filter(word %in% valid_go$GO) %>%  # Keep only valid GO terms
90         count(Gene, word) %>%
91         cast_dtm(Gene, word, n) %>%
92         removeSparseTerms(sparse = p_sparse) %>%
93         .[rowSums(as.matrix(.)) > 0, ]
94       
95       if (nrow(self$dtm) == 0) stop("DTM has zero rows after filtering")
96       invisible(self)
97       return(self$dtm)
98     },
99     
100     # 2. Parallelized Topic Optimization
101     find_optimal_k = function(k_range = 3:11, metrics = c("CaoJuan2009", "Griffiths2004")) {
102       message("Tuning topics with parallel processing...")
103       plan(multisession)
104       
105       self$k_metrics <- FindTopicsNumber(
106         self$dtm,
107         topics = k_range,
108         metrics = metrics,
109         method = "Gibbs",
110         control = list(seed = 123, keep = 50),
111         mc.cores = availableCores() - 1
112       )
113       
114       self$optimal_k <- self$k_metrics %>%
115         pivot_longer(-topics) %>%
116         group_by(name) %>%
117         mutate(scaled_value = scale(value)) %>%
118         group_by(topics) %>%
119         summarise(score = mean(scaled_value)) %>%
120         filter(score == max(score)) %>%
121         pull(topics)
122       
123       message("Optimal k: ", self$optimal_k)
124       invisible(self)
125       return(list(k_metrics = self$k_metrics, optimal_k = self$optimal_k))
126     },
127     
128     # 3. Efficient LDA Fitting
129     fit_lda = function(k = NULL, method = "Gibbs", iter = 1000, n_top_terms = 10) {
130       k <- ifelse(is.null(k), self$optimal_k, k)
131       message("Fitting LDA with k=", k, "...")
132       
133       self$lda_model <- LDA(
134         self$dtm,
135         k = k,
136         method = method,
137         control = list(
138           seed = 123,
139           keep = 50,
140           iter = iter,
141           verbose = 1
142         )
143       )
144       
145       self$doc_topics <- tidy(self$lda_model, matrix = "gamma") %>%
146         mutate(gamma = ifelse(is.nan(gamma), 0, gamma))
147       # Store top terms
148       private$extract_top_terms(n_top_terms)
149       return(self$lda_model)
150       invisible(self)
151     },
152     
153     get_top_terms = function(n = NULL) {
154       if (is.null(self$top_terms)) {
155         warning("Top terms not extracted yet. Run fit_lda() first.")
156         return(NULL)
157       }
158       if (!is.null(n)) {
159         return(self$top_terms %>% group_by(topic) %>% slice_head(n = n))
160       }
161       return(self$top_terms)
162     },
163     
164     get_formatted_top_terms = function(n = 1, wide_format = TRUE) {
165       self$terms <- self$get_top_terms(n)
166       if (wide_format) {
167         self$terms %>%
168           group_by(topic) %>%
169           summarize(
170             terms = paste(term, collapse = ", "),
171             mean_beta = mean(beta)
172           )
173       } else {
174         self$terms
175       }
176     },
177 
178     build_genotype_network = function() {
179       self$doc_topics <- tidy(self$lda_model, matrix = "gamma")
180       
181       self$g_edges <- self$doc_topics %>%
182         left_join(
183           self$master_df %>% select(Gene, genotype, log2FoldChange, identity),
184           by = c("document" = "Gene"),
185           relationship = "many-to-many"
186         )
187       print(self$g_edges)
188       self$g_nodes <- self$g_edges %>%
189         distinct(document, genotype, log2FoldChange) %>%
190         mutate(
191           importance = abs(log2FoldChange),
192           color = cut(
193             log2FoldChange,
194             breaks = c(-Inf, -1, 1, Inf),
195             labels = c("blue", "gray", "red")
196           )
197         )
198       
199       self$net_graph <- tbl_graph(
200         nodes = self$g_nodes,
201         edges = self$g_edges %>% select(from = document, to = topic, weight = gamma)
202       )
203       self$net_graph <- self$net_graph %>%
204         mutate(
205           degree = centrality_degree(),
206           betweenness = centrality_betweenness(),
207           #closeness = centrality_closeness(),
208           page_rank = centrality_pagerank(),
209           eigen = centrality_eigen(),
210           centrality_alpha = centrality_alpha()
211         )
212       return(self$net_graph)
213     }, 
214     
215     # New: Create topic-document graph
216     build_topic_document_graph = function(gamma_threshold = 0.05) {
217       if (is.null(self$doc_topics)) {
218         stop("Run build_genotype_network() first")
219       }
220       if (is.null(self$terms)) {
221         message("Formatted top terms is empty, creating tibble...")
222         self$get_formatted_top_terms()
223       }
224       
225       all_go_ids <- unique(unlist(strsplit(self$master_df$GO, split = ",\\s*")))
226       all_go_ids <- all_go_ids[!is.na(all_go_ids)]
227       
228       self$go_map <- tryCatch({
229         data.frame(lapply(self$terms, gsub, pattern = "go:", replacement = "GO:")) %>%
230           mutate(
231           TERM = AnnotationDbi::Term(term)
232         )
233       }, error = function(e) {
234         warning("GO term mapping failed: ", e$message)
235         data.frame(GOID = all_go_ids, Term = all_go_ids)
236       })
237       
238       # 1. Create all possible document-topic-genotype combinations
239       expanded_edges <- self$doc_topics %>%
240         left_join(
241           self$master_df %>% distinct(Gene, genotype),
242           by = c("document" = "Gene")
243         ) %>%
244         mutate(
245           full_id = paste(document, genotype),  # Create "cyoA WT" style IDs
246           #topic = self$go_map$TERM[match(gsub("go:", "GO:", topic), self$go_map$term)]#need to fix
247         ) %>%
248         filter(gamma > gamma_threshold) %>%
249         select(from = full_id, to = topic, weight = gamma)
250       
251       # 2. Prepare document nodes from master_df
252       doc_nodes <- self$master_df %>%
253         distinct(identity, genotype, log2FoldChange, timepoint, padj, sig_star) %>%
254         dplyr::rename(name = identity, group = genotype)
255       
256       # 3. Build graph (corrected)
257       self$topic_doc_graph <- tbl_graph(
258         nodes = bind_rows(
259           doc_nodes %>% mutate(type = "document"),
260           tibble(name = unique(expanded_edges$to), 
261                  group = "topic", 
262                  type = "topic")
263         ),
264         edges = expanded_edges
265       ) %>%
266         mutate(
267           centrality = centrality_degree(),
268           community = group_infomap()
269         )
270       
271       message(sprintf(
272         "Success! Built graph with:\n- %d documents\n- %d topics\n- %d connections",
273         sum(igraph::V(self$topic_doc_graph)$type == "document"),
274         sum(igraph::V(self$topic_doc_graph)$type == "topic"),
275         igraph::ecount(self$topic_doc_graph)
276       ))
277       return(self$topic_doc_graph)
278     },
279     
280     # Modified: Enhanced topic network with genotype connections
281     build_topic_network = function(threshold = 0.001, genotype_edge_weight = 0.5) {
282       # Original topic-topic network
283       self$topic_terms <- tidy(self$lda_model, matrix = "beta")
284       self$topic_sim <- self$topic_terms %>%
285         pairwise_similarity(topic, term, beta, upper = FALSE) %>%
286         filter(similarity > threshold)
287       
288       # Genotype-genotype edges (using identity column)
289       genotype_edges <- self$master_df %>%
290         group_by(genotype) %>%
291         summarise(genes = list(unique(identity))) %>%
292         pairwise_count(genotype, genes) %>%
293         dplyr::rename(from = item1, to = item2, weight = n) %>%
294         mutate(weight = genotype_edge_weight * weight / max(weight))
295       
296       # Combined graph
297       self$tn_graph <- tbl_graph(
298         nodes = tibble(
299           name = c(
300             paste0("Topic_", unique(c(self$topic_sim$item1, self$topic_sim$item2))),
301             unique(c(genotype_edges$from, genotype_edges$to))
302           ),
303           type = c(
304             rep("topic", length(unique(c(self$topic_sim$item1, self$topic_sim$item2)))),
305             rep("genotype", length(unique(c(genotype_edges$from, genotype_edges$to))))
306           )
307         ),
308         edges = bind_rows(
309           self$topic_sim %>% 
310             dplyr::rename(from = item1, to = item2, weight = similarity) %>%
311             mutate(
312               from = paste0("Topic_", from),
313               to = paste0("Topic_", to)
314             ),
315           genotype_edges
316         )
317       ) %>%
318         mutate(degree = centrality_degree(),
319           betweenness = centrality_betweenness(),
320           centrality_alpha = centrality_alpha())
321       
322       return(self$tn_graph)
323     },
324 
325     # Visualization for topic-document graph
326     plot_topic_document_graph = function(label_threshold = 0.5) {
327       if (is.null(self$topic_doc_graph)) {
328         stop("Run build_topic_document_graph() first")
329       }
330       
331       ggraph(self$topic_doc_graph %>% activate(edges) %>% filter(weight > 0), 
332              layout = "fr", weights = weight) +
333         geom_edge_link(aes(alpha = weight, length = weight), width = 0.1, color = "gray50") +
334         geom_node_point(aes(color = group, size = ifelse(type == "topic", 5, 3))) +
335         geom_node_label(
336           aes(label = ifelse(type == "topic" | centrality > label_threshold, name, "")), 
337           repel = TRUE
338         ) +
339         theme_graph() +
340         labs(title = "Topic-Document Network")
341     },
342     
343     # Visualization for enhanced topic network
344     plot_enhanced_topic_network = function() {
345       if (is.null(self$tn_graph)) {
346         stop("Run build_topic_network() first")
347       }
348       
349       ggraph(self$tn_graph, layout = "fr") +
350         geom_edge_link(aes(alpha = weight, color = ifelse(is.na(similarity), "genotype", "topic"))) +
351         geom_node_point(aes(size = degree, color = type)) +
352         geom_node_label(aes(label = ifelse(degree > median(degree), name, "")), repel = TRUE) +
353         scale_color_manual(values = c("topic" = "red", "genotype" = "blue")) +
354         theme_graph()
355     },
356     
357     combine_networks <- function(analyzer) {
358       # Get both graphs
359       doc_graph <- analyzer$topic_doc_graph
360       term_graph <- analyzer$tn_graph
361       
362       # Rename term graph nodes to avoid conflicts
363       term_graph <- term_graph %>%
364         activate(nodes) %>%
365         mutate(name = paste0("TERM_", name))  # Prefix term names
366       
367       # Combine graphs
368       combined_graph <- graph_join(doc_graph, term_graph)
369       
370       # Optional: Add edges between topics in both graphs
371       topic_links <- analyzer$doc_topics %>%
372         group_by(topic) %>%
373         summarise(mean_gamma = mean(gamma)) %>%
374         left_join(analyzer$topic_sim, by = c("topic" = "item1"))
375       
376       return(combined_graph)
377     }
378   ),
379     
380   private = list(
381     .valid_go = NULL,
382     
383     # Properly named and formatted private methods
384     validate_inputs = function() {
385       if (!self$go_corpus_col %in% names(self$master_df)) {
386         stop("Column '", self$go_corpus_col, "' not found")
387       }
388       if (all(is.na(self$master_df[[self$go_corpus_col]]))) {
389         stop("GO term column contains only NA values")
390       }
391     },
392     
393     validate_go_terms = function() {
394       all_go_ids <- unique(na.omit(unlist(
395         strsplit(self$master_df[[self$go_corpus_col]], ",\\s*")
396       )))
397       
398       if (!requireNamespace("org.EcK12.eg.db", quietly = TRUE)) {
399         stop("Required package org.EcK12.eg.db is not installed")
400       }
401       
402       result <- tryCatch({
403         res <- AnnotationDbi::select(
404           org.EcK12.eg.db,
405           keys = all_go_ids,
406           columns = c("GO", "TERM"),
407           keytype = "GO"
408         )
409         if (nrow(res) == 0) {
410           data.frame(GO = all_go_ids, TERM = all_go_ids)
411         } else {
412           res %>% distinct(GO, .keep_all = TRUE)
413         }
414       }, error = function(e) {
415         warning("GO validation failed: ", e$message)
416         data.frame(GO = all_go_ids, TERM = all_go_ids)
417       })
418       
419       # Ensure we always set .valid_go, even if it's just the raw IDs
420       private$.valid_go <- result
421       invisible(self)
422     },
423     extract_top_terms = function(n_top_terms) {
424       self$top_terms <- tidy(self$lda_model, matrix = "beta")
425       
426       # Add GO term definitions if available
427       val_go <- self$get_valid_go()  # Changed variable name and access method
428       if (!is.null(val_go)) {
429         self$top_terms <- self$top_terms %>%
430           left_join(val_go, by = c("term" = "GO")) %>%
431           mutate(term = ifelse(is.na(TERM), term, TERM))
432       }
433       
434       # Rank terms within topics
435       self$top_terms <- self$top_terms %>%
436         group_by(topic) %>%
437         arrange(desc(beta), .by_group = TRUE) %>%
438         mutate(rank = row_number()) %>%
439         ungroup()
440       
441       invisible(self)
442     }
443   )
444 )

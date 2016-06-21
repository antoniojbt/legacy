##############
# Generic functions for gene expression analysis from microarryas after limma processing
# Antonio Berlanga-Taylor
# 02 March 2016
##############

##############
# These functions require 
# ggplot2
# illuminaHumanv4
# XGR
##############


##############
## Add annotations to limma's topTable:
get_gene_symbols <- function(topTable_results){
  # Get probe IDs from topTable results:
  illumina_ids  <-  as.character(rownames(topTable_results))
  # Get Entrez IDs:
  probes_by_ENTREZID_RE  <- unlist(mget(illumina_ids, illuminaHumanv4ENTREZREANNOTATED, ifnotfound = NA))
  probes_by_ENTREZID_RE  <- as.data.frame(probes_by_ENTREZID_RE)
  # Get gene symbols:
  probes_by_symbol <- unlist(mget(illumina_ids, illuminaHumanv4SYMBOLREANNOTATED, ifnotfound = NA))
  probes_by_symbol <- as.data.frame(probes_by_symbol)
  # Set IDs as columns for merge:
  topTable_results$probe_ID <- row.names(topTable_results)
  probes_by_symbol$probe_ID <- row.names(probes_by_symbol)
  probes_by_ENTREZID_RE$probe_ID <- row.names(probes_by_ENTREZID_RE)
  # Merge results and annotations:
  topTable_results <- merge(topTable_results, probes_by_symbol)
  topTable_results <- merge(topTable_results, probes_by_ENTREZID_RE)
  # Order lists by adj. p-value in increasing order:
  topTable_results <- topTable_results[order(topTable_results[, 'adj.P.Val'], decreasing = FALSE), ]
  return(topTable_results)
}
##############

##############
# Get ENTREZ IDs for top results from topTable:
get_EntrezIDs_top_genes <- function(topTable_results, adjPvalue){
  # Get topTable of interest at desired p-value:
  top_genes <- which(topTable_results$adj.P.Val < adjPvalue)
  top_genes <- topTable_results[top_genes, ]
  # Set rownames, may cause errors column names are different:
  row.names(top_genes) <- top_genes$probe_ID
  # Get IDs:
  top_genes_ENTREZID  <- unlist(mget(as.character(row.names(top_genes)), illuminaHumanv4ENTREZID, ifnotfound = NA))
  return(top_genes_ENTREZID)
}
##############


##############
# Run XGR as a function. Run separtely as a loop through all ontology terms in the XGR script.
run_XGR_xEnricherGenes <- function(hits_data, background_data, ontology_term, hits_file){
  # Run enrichment analysis:
  enrichment_result <- xEnricherGenes(data = hits_data, background = background_data, ontology = ontology_term, verbose = TRUE)
  # Save enrichment results to file:
  save_enrichment <- xEnrichViewer(enrichment_result, top_num=length(enrichment_result$adjp), sortBy="adjp", details=TRUE)
  top_terms <- save_enrichment[1:20, ]
  output_enrichment <- data.frame(term=rownames(save_enrichment), save_enrichment)
  utils::write.table(output_enrichment, file=sprintf("%s_%s_enrichments.txt", hits_file, ontology_term), 
                     sep="\t", quote = FALSE, row.names=FALSE)
  # Plot:
  png(sprintf("%s_%s_enrichments.png", hits_file, ontology_term), width = 12, height = 12, units = 'in', res = 300)
  plot_gg <- ggplot(top_terms, aes(y=((nOverlap/nAnno)*100), x=name)) +
    geom_point(stat='identity', aes(colour = adjp)) +
    ylab('Percent of genes contained in the set') +
    xlab('Term') +
    coord_flip()
  print(plot_gg) # print() is silent interactively but required here to save to file
  dev.off()
  return(enrichment_result)
}
##############







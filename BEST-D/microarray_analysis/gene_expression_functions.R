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


##############
# Functions for plotting PCs:

# Loop plot for allocation group:
loop_ggplot_PCs <- function(data_file, PCa, PCb, grouping_var) {
  PCa <- sprintf('PC%s', PCa)
  PCb <- sprintf('PC%s', PCb)
  plot_num <- sprintf('%s', PCa)
  plot_num <- ggplot(data=data_file, aes(x=data_file[, c(PCa)], y=data_file[, c(PCb)], colour=grouping_var, shape=grouping_var)) + 
    theme_bw() + geom_point(alpha = 0.8) + 
    scale_colour_manual(values = c("red", "blue", "purple", "darkgreen","black", "orange")) +
    theme(legend.position="bottom", legend.title = element_blank()) +
    labs(x = PCa, y = PCb)
  # file_name <- sprintf('plot_PCA_normalised_filtered_by_allocation_%s_%s.png', PCa, PCb)
  # ggsave(filename = file_name)
  return(plot_num)
}
# PC1 <- loop_ggplot_PCs(pca_by_groups, 1, 2, allocation_group)
# PC1

# Function to loop through:
plot_PCs <- function(i) {
  file_name <- sprintf('p%s', i)
  print(file_name)
  grid_list[[i]] <- file_name
  file_name <- loop_ggplot_PCs(pca_by_groups, i, i+1, allocation_group)
  return(file_name)
}

##############

##############
# Extract the legend as an external object for use with library(gridExtra)
# http://www.sthda.com/english/wiki/ggplot2-easy-way-to-mix-multiple-graphs-on-the-same-page-r-software-and-data-visualization#cowplot-publication-ready-plots

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}
##############


##############
# Multiple plots in one page:
# http://rstudio-pubs-static.s3.amazonaws.com/2852_379274d7c5734f979e106dcf019ec46c.html
# http://www.sthda.com/english/wiki/ggplot2-easy-way-to-mix-multiple-graphs-on-the-same-page-r-software-and-data-visualization#cowplot-publication-ready-plots
# http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
# http://www.sthda.com/english/wiki/ggplot2-easy-way-to-mix-multiple-graphs-on-the-same-page-r-software-and-data-visualization

# multiplot function
# http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/ 
# This is the definition of multiplot. It can take any number of plot objects as arguments, or if it can take a list of plot objects passed to plotlist.
# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
##############

##############
# Write to disk limma's toptable:
write_topTable <- function(coef_var, fit_var) {
  write.table(topTable(fit_var, adjust="BH", coef = coef_var, number = Inf), sep='\t', quote = FALSE, 
              col.names = NA, row.names = TRUE, 
              file=paste0('topTable_', coef_var, '.txt'))  
}
##############

##############

##############


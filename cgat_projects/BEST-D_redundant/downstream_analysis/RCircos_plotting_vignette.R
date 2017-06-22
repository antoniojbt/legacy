#################
# RCircos plots
# Antonio J Berlanga-Taylor
# 18 May 2016
#################

#################
# RCircos facilitates circular plotting and interface with R
# See:
# https://cran.r-project.org/web/packages/RCircos/vignettes/Using_RCircos.pdf
# Check the demo:
# demo("RCircos.Demo.Human")
# Copied to file 'RCircos_demo.R'. The code is the same as the vignette (copied/paraphrased in this file).

# RCircos implements Circos concepts/tools in R, Krzywinski et al 2009:
# http://genome.cshlp.org/content/19/9/1639.full
# http://circos.ca/documentation/

# Other R circular plotting tools are:
# https://bioinformatics.oxfordjournals.org/content/30/19/2811.full
# http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3921174/
# http://bioinformatics.oxfordjournals.org/content/31/17/2912.full.pdf+html
#################

#################
# Load libraries
# source("https://bioconductor.org/biocLite.R")
# biocLite()
# biocLite('RCircos')
# biocLite('cowplot')

library(RCircos)
# library(cowplot) # Complements ggplot2 for publication ready figures (multi-panels, labels, etc.)
# See: https://cran.r-project.org/web/packages/cowplot/vignettes/introduction.html
#################


#################
# RCircos needs data frame as input with chr, start and end 
# Plotting is carried out by:
  # Initialising core components
  # Initialising graphics device
  # Specifying parameters
  # Specifying plot data and locations
  # Adding tracks sequentially
  # Saving to disk
#################

#################
# Examples of RCircos input data and plot types:
# Takes data frame as input with first 3 columns as chr, start, end, then data:
data(RCircos.Histogram.Data)
head(RCircos.Histogram.Data)

# For gene labels and heatmap plots, the gene/probe names must be provided in the fourth column. 
# For other plots, this column could be optional.
data(RCircos.Heatmap.Data)
head(RCircos.Heatmap.Data)


# Plot circle with links between positions:
# Input data for link line plot has only paired genomic position information for each row
# in the order of chromosome name A, chromStart A, chromEnd A, chromosome name B, chromStart B, and chromEnd B.
data(RCircos.Link.Data)
head(RCircos.Link.Data)
#################

#################
# Tracks for RCircos can be of:
# 1.  Connectors
# 2.  Gene labels
# 3.  Heatmap
# 4.  Scatterplot
# 5.  Line plot
# 6.  Histogram
# 7.  Tile plot
# 8.  Link lines
# 9.  Ribbons
# Actual track number to specify would be 10 (9 plus ideogram) 
# Order and number are according to parameters set: if inside or outside and sequence in 
# which they were added:
#################

#################
# Initialize RCircos core components:
# Load the chromosome ideogram data:
data(UCSC.HG19.Human.CytoBandIdeogram)
head(UCSC.HG19.Human.CytoBandIdeogram)
dim(UCSC.HG19.Human.CytoBandIdeogram)
colnames(UCSC.HG19.Human.CytoBandIdeogram)

# Load ideograms from RCircos itself or use a text file manually. Download cytoBandIdeo.txt.gz files from e.g.:
# http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/
# ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/cytoBandIdeo.txt.gz
# e.g.:
hg38 <- read.csv(file = 'cytoBandIdeo.txt', sep = '\t', header = FALSE, 
                 col.names = c("Chromosome", "ChromStart", "ChromEnd", "Band", "Stain"))
head(hg38)
dim(hg38)

# RCircos core components can now be initialized with RCircos.Set.Core.Components() function, which requires:
# cytoinfo = the chromosome ideogram data loaded
# chr.exclude = chromosomes which should be excluded from ideogram, e.g.
# chr.exclude <- c(”chrX”, ”chrY”) # Or set to NULL to include all chrs
# tracks.inside = how many tracks will be plotted inside chromosome ideogram
# tracks.outside = how many data tracks will be plotted outside chromosome ideogram

# An example which initializes core components with all human chromosome ideogram and 
# 10 data track spaces inside of the chromosome ideogram:
chr.exclude <- NULL
cyto.info <- UCSC.HG19.Human.CytoBandIdeogram
tracks.inside <- 10
tracks.outside <- 0
RCircos.Set.Core.Components(cyto.info, chr.exclude, tracks.inside, tracks.outside)

# RCircos core components are saved and can be retrieved with:
rcircos.params <- RCircos.Get.Plot.Parameters(); # Can be user set
rcircos.cyto <- RCircos.Get.Plot.Ideogram(); # Selected by user
rcircos.position <- RCircos.Get.Plot.Positions(); # Modified according to parameters
RCircos.List.Parameters()

# Use RCircos plot parameters to modify plot aspects such as color.map, text.size, etc. To modify use:
rcircos.params <- RCircos.Get.Plot.Parameters();
rcircos.params$base.per.unit <- 30000;
RCircos.Reset.Plot.Parameters(rcircos.params);
RCircos.List.Parameters()
#################

#################
# Worked example (same as demo("RCircos.Demo.Human")):
# So, to plot with RCircos initialise core components and parameters, then add tracks sequentially:
# Initialize Graphics Device:
out.file <- "RCircosDemoHumanGenome.pdf";
pdf(file=out.file, height=8, width=8, compress=TRUE);
RCircos.Set.Plot.Area();

#   Setup RCircos core components:
#   1.  Chromosome ideogram plot information
#   2.  x and y coordinates for a circular line and degrees of the text rotation at each point
#   3.  Plot parameters for image plot control
RCircos.Set.Core.Components(cyto.info=hg19.cyto, chr.exclude=NULL, tracks.inside=10, tracks.outside=0)

# RCircos.Set.Plot.Area() will setup plot area base on total number of
# tracks inside and outside of chromosome ideogram. Modify this with:
par(mai=c(0.25, 0.25, 0.25, 0.25));
plot.new();
plot.window(c(-2.5,2.5), c(-2.5, 2.5));

# Plot Chromosome Ideogram (for example):
RCircos.Chromosome.Ideogram.Plot();
#################

#################
# Continue example:
# Gene Labels and connectors on RCircos Plots
# For best visualization, cex should be >= 0.4 when drawing gene labels

# Connectors are used to mark a genomic position with their names or variant status. 
# Currently, RCircos only provide connector plot between genes and their genomic positions.
data(RCircos.Gene.Label.Data);
head(RCircos.Gene.Label.Data);
name.col <- 4;
side <- "in";
track.num <- 1;
RCircos.Gene.Connector.Plot(RCircos.Gene.Label.Data, track.num, side);
track.num <- 2;
RCircos.Gene.Name.Plot(RCircos.Gene.Label.Data, name.col,track.num, side);
#################

#################
# Continue example:
# Heatmap, Histogram, Line, Scatter, and Tile Plot
# First three columns of input data are chromosome name, start, and end position. 
# Provide plot location (which track and which side of chromosome ideogram)

# Heatmap:
# color.map = Color scales for heatmap plot, default: BlueWhiteRed
data(RCircos.Heatmap.Data);
data.col <- 6; # Data to be plotted
track.num <- 5; # Plot location
side <- "in"; # Determine where track will be (inside or outside, in relation to chr ideogram for example)
RCircos.Heatmap.Plot(RCircos.Heatmap.Data, data.col, track.num, side); # Plot command

# Scatter:
data(RCircos.Scatter.Data);
data.col <- 5;
track.num <- 6;
side <- "in";
by.fold <- 1;
RCircos.Scatter.Plot(RCircos.Scatter.Data, data.col, track.num, side, by.fold);

# Line:
data(RCircos.Line.Data);
data.col <- 5;
track.num <- 7;
side <- "in";
RCircos.Line.Plot(RCircos.Line.Data, data.col, track.num, side);

# Histogram:
data(RCircos.Histogram.Data);
data.col <- 4;
track.num <- 8;
side <- "in";
RCircos.Histogram.Plot(RCircos.Histogram.Data, data.col, track.num, side);

# And tile plot:
data(RCircos.Tile.Data);
track.num <- 9;
side <- "in";
RCircos.Tile.Plot(RCircos.Tile.Data, track.num, side);
#################

#################
# Continue example:
# Links: A Special Plot
# Links presents relationship of two genomic positions and it is always the last
# track inside chromosome ideogram. 
# Input data is a data frame with paired genomic positions in the order of
# chromosome, start, and end position for each one genomic position. 
# RCircos supports two types of link plots: lines and ribbons. 
# Link lines are used for presenting relationship of two small genomic regions 
# ribbons are plotted for bigger genomic regions. 

# chromosome=TRUE (or FALSE) # Modifies colours for same/different chromosome links
data(RCircos.Link.Data);
track.num <- 11;
RCircos.Link.Plot(RCircos.Link.Data, track.num, TRUE);
data(RCircos.Ribbon.Data);
RCircos.Ribbon.Plot(ribbon.data=RCircos.Ribbon.Data, track.num=11, by.chromosome=FALSE, twist=FALSE);

# Close graphics device and end plot:
dev.off();
#################


#################
# The end:
q()
#################
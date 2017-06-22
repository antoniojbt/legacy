# -*- coding: utf-8 -*-
################################################
##1.0 Genome assembly and annotation pipeline ##
################################################

from ruffus import *

import sys, os, csv

import CGAT.IOTools as IOTools

import CGAT.Pipeline as P

import CGAT.Experiment as E

#To do:
#bash compile Cortex
#design calculation for memory, genome complexity, etc.
#readQC
#run case studies
#run available Pseudomonas strains

###################################################
## Pipeline configuration
###################################################

# load options from the config file

P.getParameters(
                ["pipeline.ini" ]
                )

PARAMS = P.PARAMS

PARAMS_ANNOTATIONS = P.peekParameters( PARAMS["annotations_dir"],
                                       "pipeline_annotations.py" )


###################################################
############1.1 Download data


#Create directories and download data

@follows(mkdir('sra-Download.dir'))
@files(None, 'wget-log.sentinel')
def downloadData(infile, outfile): 
    to_cluster = False
    statement = ''' wget -r 
			ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/ERX/ERX003/ERX003921/ERR010476/
			--no-directories
			-o %(outfile)s
			--directory-prefix=sra-Download.dir;	
		''' % locals()
    P.run()


#@follows(downloadData)
#@check_if_uptodate(downloadData)
#def checkDownload():
#    try:
#	with open('wget-log'): pass
#    except IOError:
#    	sys.exit(P.main(sys.argv))
	
#    P.run()

#sra to fastq - split

@follows(downloadData)
@transform("sra-Download.dir/*.sra", regex(r".*/(.*).sra"), r"\1.sra.sentinel")

def sraToFastq(infile, outfile):

    to_cluster = True
    statement = ''' fastq-dump.2.1.7
				--split-3 
                                %(infile)s;
                touch %(outfile)s
                ''' % locals()

    P.run()


#############1.2 Examine quality

#run readQC

#readqc - pipeline

#@follows(sraToFastq, mkdir("readqc-pipeline.dir"))
#def readqcPipeline(infile, outfile):

#    to_cluster = False
#    statement = ''' python 
#                                      %(infile)s;
#               touch %(outfile)s
#                ''' % locals()

#    P.run()


#############1.3 Run genome assemblers

############
#Run Cortex
############

#Decompress files for Cortex

#@follows(readqc pipeline)
#@transform("*.gz", suffix(".fastq.gz"), ".fastq")

#def gunzip(infile, outfile):

#    to_cluster = True
#    statement = ''' gunzip %(infile)s;
#                ''' % locals()

#    P.run()
 
#fastq to vcf - cortex
	#run_calls.pl
		#create index file and sample lists
                #common mistakes so far: missing dots in index file; wrong tabs, structure in index file; different lengths 

@follows(sraToFastq, 'velvetOptimiser')
@files('*.fastq', 'sample-list.tab')
def makeSampleList(infile, outfile):

    outf = IOTools.openFile('outfile2.sentinel', 'w')
    infile = '\n'.join(infile) + '\n'
    outf.write(infile)
    outf.close()
    outf2 = open(outfile, 'w')
    with open('outfile2.sentinel', 'r') as f:
    	for line in f:
		outf2.write(line[:])
    outf2.close()
	

@follows(makeSampleList)
@files('sample-list.tab', 'index.tab')
def makeIndex(infile, outfile):
    destination = open(outfile, 'w')
    source = open(infile,'r')
    for line in source:
        destination.write(line[:-7].strip('\n') + '\t' + 'sample-list.tab' + '\t' + '.' + '\t' + '.' + '\n')
    source.close()
    destination.close()
 

#calculate memory use

#height required for Cortex command
#width required for Cortex command
#Genome size (bp) must be specified by the user
#guess (bp) is an initial number to run, Cortex will stop if this is too low. \
#It serves to allocate the right amount of memory though a few tries may be needed.

#M = memory use calculation (bytes)
#N = total number of kmers
#c = number of samples (colours)
#k = kmer length; has to be specified by user

#@follows(makeIndex)

#height = 0
#width = 0
#genomeSize = 6200000
#guess = genomeSize * 2
#estimate = (2**height) * (width) 
#approximate = guess - estimate
#N = 
#c = 
#k = 
#M = N*((k/32)+(5*c)+1)


#def memoryGuess():
   

#return height, width

#def passParametersToCortex(height, width):



#if using reference build Stampy index
#Change extra_parameters so that multiple file references can be passed
#Run_calls.pl only takes filename.sthash and filename.stidx (ie filename.fasta.st* throws an error)

@follows(makeIndex)
#@transform('*.fasta', regex(r".*/(.*).fasta"), r"\1.stidx")
@transform('*.fasta', suffix('.fasta'), '.stidx', 'PAO1')
def stampyIndex(infile, outfile, extra_parameters1):

    to_cluster = False
    statement = ''' stampy.py -G %(extra_parameters1)s %(infile)s;
                ''' % locals()

    P.run()

# if using reference build Stampy hash

@follows(stampyIndex)
#@transform('*.fasta', regex(r".*/(.*).fasta"), r"\1.sthash")
@transform('*.fasta', suffix('.fasta'), '.sthash', 'PAO1')
def stampyHash(infile, outfile, extra_parameters1):

    to_cluster = False
    statement = ''' stampy.py -g %(extra_parameters1)s -H %(extra_parameters1)s;
                ''' % locals()

    P.run()


#To do: pass options for cortex executable, build if missing, pass kmer_size option, re-run if two kmer lengths;
#
#Create reference list for Cortex, at the moment only for single end sequenced files, needs .pelist:    

@follows(stampyHash)
@transform('*.fasta', suffix('.fasta'), '.selist')
def referenceList(infile, outfile):
    
    outf = IOTools.openFile(outfile, 'w')
    infile = ''.join(infile) + '\n'
    outf.write(infile)
    outf.close()

#build Cortex binary at desired kmer lengths

@follows(referenceList)
@transform('*.selist', suffix('.selist'), '.k31.ctx')
def buildCortexKmer_1(infile, outfile):

    to_cluster = False
    statement = ''' cortex_var_31_c1 
                    --kmer_size 31 
                    --mem_height 17 
                    --mem_width 100 
                    --se_list %(infile)s  
                    --dump_binary %(outfile)s 
                    --sample_id REF
                ''' % locals()

    P.run()


@follows(buildCortexKmer_1)
@transform('*.selist', suffix('.selist'), '.k61.ctx')
def buildCortexKmer_2(infile, outfile):

    to_cluster = False
    statement = ''' cortex_var_63_c1 
                    --kmer_size 61 
                    --mem_height 17 
                    --mem_width 100 
                    --se_list %(infile)s
                    --dump_binary %(outfile)s 
                    --sample_id REF
                ''' % locals()

    P.run()
    
#Run run_calls.pl from Cortex; specify options, refer to manual
#Delete run_calls.dir and log file to re-run
#Make so that run_calls can be re-run without having to generate binaries again

@follows(buildCortexKmer_2)
@files('index.tab', 'run_calls.dir', 'PAO1-ref-Klockgether.', 'log.run_calls.sentinel')
def doRunCalls(infile, outfile, extra_parameters1, extra_parameters2):

    to_cluster = False
    statement = ''' run_calls.pl 
                    --first_kmer 31 
                    --last_kmer 61 
                    --kmer_step 30 
                    --fastaq_index %(infile)s 
                    --auto_cleaning yes 
                    --bc yes 
                    --pd no 
                    --outdir %(outfile)s 
                    --outvcf %(extra_parameters1)s 
                    --ploidy 1 
                    --stampy_hash PAO1 
                    --stampy_bin /ifs/apps/bio/stampy-1.0.17/stampy.py 
                    --list_ref_fasta PAO1.selist 
                    --refbindir ./ 
                    --genome_size 6300000 
                    --qthresh 5 
                    --mem_height 18 
                    --mem_width 100 
                    --vcftools_dir /ifs/apps/src/vcftools_0.1.11/vcftools_0.1.11/ 
                    --do_union yes 
                    --logfile %(extra_parameters2)s 
                    --workflow joint 
                    --ref CoordinatesAndInCalling 
                    --squeeze_mem 
                    --apply_pop_classifier 
                    ''' % locals()

    P.run()

#Convert binaries to fasta files from cortex

#for sample lists; leaving a log file with the same name when re-running; \
#error in memory calculation (note where the error is located to re-run automatically);  

#load raw and decomp Cortex vcfs to database

#summary of vcf?

#Remove sites from VCFs that have coverage on both alleles (most likely to be repeats)


#######################
#Run Velvet
#######################

#Velvet web page: http://www.ebi.ac.uk/~zerbino/velvet/Manual.pdf
#To run with VelvetOptimiser: http://www.vicbioinformatics.com/velvetoptimiser.manual.txt
#Usage: VelvetOptimiser.pl [options] -f 'velveth input line'
#Only required parameter is -f; -o passes additional velvetg parameters
#VelvetOptimiser needs velveth -f option without hash or directory name as:
#{[-file_format][-read_type] filename} repeated for as many read files as you have.
#Turn this into a passable parameter so as to run for as many read files as present in a directory

#When setting parameters, velveth needs:
#the file format (fasta, fastq, fastq.gz, etc.), read category (short -Illumina-, paired, etc.)
#-separate indicates read pairs are in separate files

#Not yet: Currently set for a low starting hash value (-s) and a high end hash value (-e) with 


#velvetg options can include: -clean yes; 
 
#Example: VelvetOptimiser.pl -s 33 -e 41 -f '-fastq.gz -shortPaired -separate SRR292770_1.fastq.gz SRR292770_2.fastq.gz' \
#                            -o ' -min_contig_lgth 200' -p SRR292770

@follows(sraToFastq)
@transform('*.fastq', suffix('.fastq'), '.VelvetOptimiser.dir/')

def velvetOptimiser(infile, outfile): 
    to_cluster = True
    statement = ''' VelvetOptimiser.pl 
                    -s 33 
                    -e 41 
                    -f '-fastq -short %(infile)s' 
                    -o '-min_contig_lgth 200 -clean yes'
                    -d %(outfile)s
                ''' % locals()
    P.run()


#@follows(sraToFastq, mkdir("VelvetOptimiser.dir))
#@transform('*.fastq', regex('(\S+).fastq'), r'VelvetOptimiser.dir/\1.contigs.fa')
#
#def velvetOptimiser(infile, outfile): 
#    to_cluster = True
#    
#    tempdir = P.getTemporaryDirectory()
#    statement = ''' VelvetOptimiser.pl 
#                    -s 33 
#                    -e 41 
#                    -f '-fastq -short %(infile)s' 
#                    -o '-min_contig_lgth 200 -clean yes'
#                    -d %(tempdir)s/contig.fa
#                    ; checkpoint; mv %(tempdir)s/contigs.fa VelvetOptimiser.dir/%(outfile)s 
#                ; rm -rf %(tempdir)s
#                ''' % locals()
#    P.run()

#Run velveth if not using VelvetOptimiser, together with KmerGenie: http://kmergenie.bx.psu.edu/


#Example: velveth out_data_35 35 -fastq.gz -shortPaired -separate SRR292770_1.fastq.gz SRR292770_2.fastq.gz
#To do: run function to check file input (fastq, fastq.gz, etc.), type of reads (short, short paired, long, etc.), if input files are paired and  \
#separate (pair_1.fastq, pair_2.fastq, etc.)
#Set k-mer length as parameter (=extra parameter) using the results from KmerGenie


@follows(sraToFastq, mkdir('velveth.dir'))
@transform('*.fastq', suffix('.fastq'), '/velveth.dir/*', '35')

def runVelvetH(infile, outfile, extra_parameters1): 
    to_cluster = True
    statement = ''' velveth %(outfile)s %(extra_parameters1)s 
                    -fastq
                    -short %(infile)s
		''' % locals()
    P.run()


#Run velvetg if not using VelvetOptimiser:

#The next Velvet step to run is velvetg to build the graph.
#Example: velvetg out_data_35 -clean yes -exp_cov 21 -cov_cutoff 2.81 -min_contig_lgth 200

#@follows(mkdir('xxxxx.dir'))
@follows(runVelvetH)
@transform('velveth.dir', suffix(''), 'velveth.dir/*.fa')

def runVelvetG(infile, outfile): 
    to_cluster = True
    statement = ''' velvetg %(infile)s
                    -clean yes 
                    -min_contig_lgth 200
		''' % locals()
    P.run()

#To do: check the Velvet Columbus module in order to assemble with a reference that assists the process, http://www.ebi.ac.uk/~zerbino/velvet/Columbus_manual.pdf

#Copy the contigs file from the Velvet output folder and rename it:
#Reason? Example: cp out_data_35/contigs.fa SRR292770_unordered.fasta .
#Need to run regex for multiple files.

#@follows(runVelvetG, mkdir('xxxxx.dir'))
#@files(None, 'xxxx.sentinel')
#
#def copyVelvetFiles(infile, outfile): 
#    to_cluster = True
#    statement = ''' cp %(infile)s %(outfile)s .
#		''' % locals()
#    P.run()



#############1.4 Order contigs against a reference

#After assembling into contigs order them against a reference genome for further downstream analysis and comparison.
#‘Move Contigs’ from Mauve can do this:
#Reference:
#Darling, A. E., Mau, B. and Perna, N. T., “progressiveMauve: multiple genome alignment with gene gain, loss and rearrangement”. \
#PLoS One,2010 5(6):e11147.
#Alternative is Abacas (MUMmer wrapper). Abacas:
#http://abacas.sourceforge.net/Manual.html#7._Default_output_files

#The command for ordering against a reference genome is: 
#Example: abacas.1.3.1.pl –r 'reference_file'.fasta -q 'query_file'_unordered.fasta –p ‘nucmer’ -c –m –b –o 'output_file'.fasta
#Example: abacas -r *.fasta -q VelvetOptimiser.dir/*.fa -p nucmer -b -m -o abacas
#nucmer files will be rewritten at this point; need to use Pipeline.py utilitites for renaming and moving to temp files
#Default output files for abacas:

#Running ABACAS with default options will generate the following files:
#    1. Ordered and orientated sequence file (reference_query.fasta or prefix.fasta)
#    2. Feature file (reference_query.tab or prefix.tab)
#    3. Bin file that contains contigs that are not used in ordering (reference_query.bin or prefix.bin)
#    4. Comparison file (reference_query.crunch or prefix.crunch)
#    5. Gap information (reference_query.gaps, prefix.gaps)
#    6. Information on contigs that have a mapping information but could not be used in the ordering (unused_contigs.out)
#    7. Feature file to view contigs with ambiguous mapping (reference.notMapped.contigs.tab). This file should be uploaded on the reference side of ACT view.
#    8. A file that shows how repetitive the reference genome is (reference.Repeats.plot).
#Files 7 & 8 should be uploaded on the reference side of ACT view. 

@follows(mkdir('abacas.dir'))
#@posttask(touch_file('abacas.sentinel'))
@transform('*301.dir/contigs.fa', suffix('.fa'), '.abacas', 'PAO1.fasta')

def runAbacas(infile, outfile, extra_parameters1):
    to_cluster = True
    statement = ''' abacas
                    -r %(extra_parameters1)s 
                    -q %(infile)s
                    -p nucmer
                    -b -m 
                    -o %(outfile)s
                    ; checkpoint
                    ; mv %(outfile)s abacas.dir/
                    ; mv nucmer* abacas.dir/
                    ; mv unused_contigs.out abacas.dir/
                    ; touch %(outfile)s
                    ; checkpoint
                ''' % locals()
    P.run()

#############1.4.1 Viewing the ordered contigs (Mauve), a java based window tool
#Mauve: Genome alignments can identify evolutionary changes in the DNA by aligning homologous regions of sequence. 
#Mauve is a software package that attempts to align orthologous and xenologous regions among two or more genome sequences that have undergone \
#both local and large-scale changes.
#For CLI and usage check: http://asap.ahabs.wisc.edu/mauve-aligner/mauve-user-guide/using-progressivemauve-from-the-command-line.html


@follows(runAbacas, mkdir('progressiveMauve.dir'))
#@posttask(touch_file('progressiveMauve.sentinel'))
@transform('*.VelvetOptimiser.dir/contigs.fa', regex('(\S+).VelvetOptimiser.dir/contigs.fa'), r'\1.progressiveMauve')

def viewMauve(infile, outfile): 
    to_cluster = True
    statement = ''' progressiveMauve
                    --output=%(outfile)s.xmfa
                    --output-guide-tree=%(outfile)s.tree --backbone-output=%(outfile)s.backbone
                    %(infile)s
                    ; checkpoint
                    ; mv %(outfile)s* progressiveMauve.dir/
                    ; touch %(outfile)s
                    ; checkpoint 
                ''' % locals()
    P.run()


#############1.4.2 Viewing the ordered contigs with Artemis, Artemis Comparison Tool and DNA plotter (Artemis suite)
#Artemis runs with art, ACT with act, and dnaplotter as such. ACT runs on Artemis.
#Check http://www.sanger.ac.uk/resources/software/act/
#ACT is only the viewer, comparison files need to be run externally.
#Use BLAST first to identify if an unknown sequence already exists in a public database:

#    Megablast is intended for comparing a query to closely related sequences and works best if the target percent identity is 95% or more but is very fast.
#    Discontiguous megablast uses an initial seed that ignores some bases (allowing mismatches) and is intended for cross-species comparisons.
#    BlastN is slow, but allows a word-size down to seven bases.

#Check:
#Web based entry:    http://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&BLAST_PROGRAMS=megaBlast&PAGE_TYPE=BlastSearch&SHOW_DEFAULTS=on&LINK_LOC=blasthome
#Manual:    http://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=ProgSelectionGuide#tab31

#blastn requires a configuration file ('.ncbirc') which can be located in the working directory. Specify -db path there and other variables.
#BLAST searches can take time so create a database first (may need to specify path in a .ncbirc file in your woring directory)
#Create a database for BLAST using makeblastdb after concatenating all genomes, example:
#makeblastdb -in ../all_bacteria_multifasta.fna -dbtype nucl -title all_bacteria -input_type fasta

#There's something going on when running makeblastdb and directory structure when using .ncbirc file and -db option for database searching (blastn looks one directory above \
#the one specified when searching for the database files). 
#Can run with -db and full path though, for example:
#blastn -query ERR010476.VelvetOptimiser.dir/contigs.fa \
#-db /ifs/mirror/ncbi/bacteria/blastdb/all_bacteria_multifasta.fna -task megablast \
#-evalue 1 -out test4.megablast -outfmt "6 qseqid qstart qend sseqid sstart send evalue bitscore pident score qseq sseq" -num_threads 8

#If running against a downloaded ncbi genome i.e. 
#/ifs/mirror/ncbi/bacteria/Pseudomonas_aeruginosa_DK2_uid168996/NC_018080.fna
#use -subject option instead of -db

#Also run against NCBI protein (nucleotide to protein) databases, other organisms, and check UCSB and Ensembl for additional genomes (?)
 

@follows(mkdir('ACT.dir'))
@transform('*301.dir/contigs.fa', suffix('.fa'), 'ACT.dir/megablast', 'megablast_vs_nt_bacteriaNCBI.db', '/ifs/mirror/ncbi/bacteria/all_bacteria_multifasta.fna')

def runMegaBlast(infile, outfile, extra_parameters1, extra_parameters2): 
    
    to_cluster = True
    statement = ''' python %(scriptsdir)s/farm.py --split-at-regex="^>(\S+)" --chunksize=10 
                    "
                    blastn 
                    -query %(infile)s
                    -db %(extra_parameters2)s
                    -task megablast
                    -evalue 1 
                    -out %(outfile)s 
                    -outfmt "6 qseqid qstart qend sseqid sstart send evalue bitscore pident score qseq sseq"
                    "
                    |
                    python %(scriptsdir)s/blast2table.py > %(extra_parameters1)s
                '''
    P.run()

# % locals() is substituted in by Pipeline.py

@follows(runMegaBlast)
@transform('ACT.dir/*.megablast', suffix('.megablast'), 'ACT.dir/.act')

def viewACT(infile, outfile): 
    to_cluster = True
    statement = ''' act
                    -quiet
                    -Dblack_belt_mode=false
                    %(infile)s 
                    %(outfile)s
                ''' % locals()
    P.run()


#@follows(viewACT)
#@transform('VelvetOptimiser.dir/*.fa', suffix('.fa'), '.dnaplotter')
#
#def viewDNAPlotter(infile, outfile): 
#    to_cluster = True
#    statement = ''' dnaplotter
#                    %(infile)s
#                    %(outfile)s
#                ''' % locals()
#    P.run()

#############1.5 Mauve Assembly Metrics –Statistical View of the Contigs
#http://code.google.com/p/ngopt/wiki/How_To_Score_Genome_Assemblies_with_Mauve

#1. Score each assembly
#2. Generate summary plots comparing the assemblies
#3. Visualize PDFs of accuracy plots. A large number of PDFs will be generated in the current directory. \
#These plots can be used to decide on the assembly strategy that best meets your project's objectives. \
#For example, you may prefer the assembly with the highest fraction of coding regions reconstructed, \
#or the assembly with the fewest structural errors, or the assembly that captures the true genomic content most precisely.
#4. Examine the table of assembly characteristics In addition to the visual output, a tab-delimited table containing assembly metrics for each \
#assembly is reported in the file summaries.txt. This file could be used in conjunction with an automated method to select assembly \
#parameters, for example. 

#Also check the Assemblathon metrics: http://genome.cshlp.org/content/21/12/2224/T7.expansion.html

#Mauve assembly metrics command example: 

# run the scoring:
#[java -cp Mauve.jar org.gel.mauve.assembly.ScoreAssembly] -reference Haloferax_volcanii_DS2.gbk \
#                -assembly assembly1.fasta -reorder assembly1_scores -outputDir assembly1_scores

# Generate plots of metrics
#./mauveAssemblyMetrics.pl ./

@follows(velvetOptimiser, mkdir('Mauve_Assembly_Metrics.dir'))
@transform('VelvetOptimiser.dir/*.fa', suffix('.fa'), 'Mauve_Assembly_Metrics.dir/*.txt')

def mauveAssemblyMetrics(infile, outfile): 
    to_cluster = True
    statement = ''' mauveAssemblyMetrics
                    -reference 
                    -assembly %(infile)s
                    -reorder
                    -outputDir Mauve_Assembly_Metrics.dir
                    ; checkpoint
                    ; mauveAssemblyMetrics .
                ''' % locals()
    P.run()


##############1.6 Annotation using RAST -- web based server

##############1.6.1 Alternatives to RAST -- for the command line, Prokka or BG7:

#A number of command line tools are available for annotation on a local machine.
#For fast de novo annotation try Prokka
#http://www.vicbioinformatics.com/software.prokka.shtml

@follows(viewACT)
#@transform('VelvetOptimiser.dir/*.fa', suffix('.fa'), 'Mauve_Assembly_Metrics.dir/*.txt')

def annotateGenomesProkka(infile, outfile): 
    to_cluster = True
    statement = ''' Prokka
                    %(infile)s
                    %(outfile)s
                ''' % locals()
#    P.run()


#For comparative annotation try BG7 
#http://bg7.ohnosequences.com/

@follows(annotateGenomesProkka)
#@transform('VelvetOptimiser.dir/*.fa', suffix('.fa'), 'Mauve_Assembly_Metrics.dir/*.txt')

def compareGenomesBG7(infile, outfile): 
    to_cluster = True
    statement = ''' BG7
                    %(infile)s
                    %(outfile)s
                ''' % locals()
#    P.run()



##########################################
#2.Comparative genome analysis
##########################################
#Most of these tools are graphics based, check CLI options 

##############2.1 Identify and download genome sequences to compare to

@follows(compareGenomesBG7)
#@transform('VelvetOptimiser.dir/*.fa', suffix('.fa'), 'Mauve_Assembly_Metrics.dir/*.txt')

def downloadGenomesForComparison(infile, outfile): 
    to_cluster = True
    statement = ''' wget 
                    %(infile)s
                    %(outfile)s
                ''' % locals()
#    P.run()


##############2.2 Use Mauve for multiple genome alignment

@follows(downloadGenomesForComparison)
#@transform('VelvetOptimiser.dir/*.fa', suffix('.fa'), 'Mauve_Assembly_Metrics.dir/*.txt')

def multipleGenomeAlignmentMauve(infile, outfile): 
    to_cluster = True
    statement = ''' Mauve
                    %(infile)s
                    %(outfile)s
                ''' % locals()
#    P.run()


##############2.3 Use ACT for detailed pairwise comparisons

@follows(multipleGenomeAlignmentMauve)
#@transform('VelvetOptimiser.dir/*.fa', suffix('.fa'), 'Mauve_Assembly_Metrics.dir/*.txt')

def pairwiseComparisonsACT(infile, outfile): 
    to_cluster = True
    statement = ''' ACT
                    %(infile)s
                    %(outfile)s
                ''' % locals()
#    P.run()


##############2.4 Use BRIG for to visualise comparisons against reference for multiple sequences

@follows(pairwiseComparisonsACT)
#@transform('VelvetOptimiser.dir/*.fa', suffix('.fa'), 'Mauve_Assembly_Metrics.dir/*.txt')

def visualiseSequencesAgainstReferenceBRIG(infile, outfile): 
    to_cluster = True
    statement = ''' BRIG
                    %(infile)s
                    %(outfile)s
                ''' % locals()
#    P.run()


##########################################
#3.Typing and specialist tools
##########################################
## Many of these tools are web based servers or visualising tools. 

##############3.1 Use PHAST to identify phage sequences

## For web based services check if they can be made local or sue their URLAPI:
# PHAST, for example, has:
# http://phast.wishartlab.com/how_to_use.html
# PHAST's URLAPI can be intergrated into a user's local program:
#1) For a GenBank ACCESSION number or a GI number use:
#http://phast.wishartlab.com/cgi-bin/phage_command_line.cgi?acc=<NC_number or GI_number>  (The '<' and '>' is not needed)
#or
#2) If you wish to submit a fasta sequence then use:
#http://phast.wishartlab.com/cgi-bin/phage_command_line.cgi?fasta_seq=<Fasta_sequence>
#
#You can use any html request tool, such as wget:
#Regex or multiple files?

@follows(doRunCalls, velvetOptimiser, mkdir('PHAST.dir'))
@transform('*.fa', suffix('.fa'), '/PHAST.dir/*.phast')

def sendPHAST(infile, outfile): 
    to_cluster = True
    statement = ''' wget 
                        http://phast.wishartlab.com/cgi-bin/phage_command_line.cgi?acc=%(infile)s 
			-O %(outfile)s
		''' % locals()
    P.run()

#Our URLAPI can only handle single submissions at a time. 
#It will not work for multiple raw DNA sequences in one input file, multiple ACCESSION numbers or GI numbers \
#in one input file, or mutilple GenBank format files in a single input file.


##############3.2 Use ResFinder to identify resistance genes


##############3.3 Perform multilocus sequence typing for database searching


##############3.4 Use PATRIC to make comparisons: http://patricbrc.org/portal/portal/patric/Home




if __name__=="__main__":

    sys.exit(P.main(sys.argv))



##########################################
#.Fasta statistics
##########################################
#Check fasta and other possibly useful tools in CGAT scripts:
#http://www.cgat.org/~andreas/documentation/cgat/CGATReference.html

#Check:
#index_fasta.py - Index fasta formatted files; Build an index for a fasta file. Pre-requisite for many CGAT tools.
#fasta2bed.py - compute G+C profile for a DNA sequence
#fasta2counts.py - basic stats from collection of sequences
#fasta2gaps.py - output assembly gaps from fasta file
#fasta2spliced.py - extract splice junctions from fasta file
#fasta2table.py - analyze sequences for codon bias and other sequence properties
#fasta2variants.py - create sequence variants from a set of sequences
#fastq2table.py - compute stats on fastq files
#vcf2vcf.py - manipulate vcf files
#diff_fasta.py - compare contents of two fasta files
#concatenate_sequences.py - concatenate sequences from multiple fasta files

#Check also blast scripts, system scripts, trees and orthology scripts
#Check also CGAT modules for phylogeny and other tools: http://www.cgat.org/~andreas/documentation/cgat/modules.html#modules


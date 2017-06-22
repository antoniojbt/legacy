#################################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id: pipeline_snps.py 2870 2010-03-03 10:20:29Z andreas $
#
#   Copyright (C) 2009 Andreas Heger
#
#   This program is free software; you can redistribute it and/or
#   modify it under the terms of the GNU General Public License
#   as published by the Free Software Foundation; either version 2
#   of the License, or (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
#################################################################################
"""
===========================
Pipeline genome_assembly
===========================

:Author: Antonio J Berlanga-Taylor
:Release: $Id$
:Date: |today|
:Tags: Python, genome assembly



The genome assembly pipeline takes reads from single genomes from one or more experiments. It takes fastq files as input, runs quality assessment for assembly purporses, builds an assembly (outputs fasta files), measures various assembly metrics and runs assembly improvement tools (aligns, breaks misassemblies, gap fills, error corrects and annotates de novo or reference based).

Overview
========

The pipeline assumes data derives from NGS experiments with purified (single) genomes which derive from multiple experiments (:term:'experiment') with one or more biological and/or technical replicates (:term:'replicate'). A :term:'replicate' within each :term:'experiment' is a :term:'track'. 

Usage
=====

See :ref:`PipelineSettingUp` and :ref:`PipelineRunning` on general information how to use CGAT pipelines.


Configuration
-------------

The pipeline requires a configured :file:`pipeline.ini` file. 

The sphinxreport report requires a :file:`conf.py` and :file:`sphinxreport.ini` file 
(see :ref:`PipelineDocumenation`). To start with, use the files supplied with the
:ref:`Example` data.

Input
-----

Optional inputs
+++++++++++++++

Requirements
------------

The pipeline requires the results from :doc:`pipeline_annotations`. Set the configuration variable 
:py:data:`annotations_database` and :py:data:`annotations_dir`.

On top of the default CGAT setup, the pipeline requires the following software to be in the 
path:

+--------------------+-------------------+------------------------------------------------+
|*Program*           |*Version*          |*Purpose*                                       |
+--------------------+-------------------+------------------------------------------------+
| Flash              |                   |                                                |
+--------------------+-------------------+------------------------------------------------+
| Velvet             |                   |                                                |
+--------------------+-------------------+------------------------------------------------+
| VelvetOptimiser    |                   |                                                |
+--------------------+-------------------+------------------------------------------------+
| Cortex             |                   |                                                |
+--------------------+-------------------+------------------------------------------------+
| Mauve Assembly Metrics |               |                                                |
+--------------------+-------------------+------------------------------------------------+
| Stampy             |                   |                                                |
+--------------------+-------------------+------------------------------------------------+
|Lots, I need to update this             |                                                |
+--------------------+-------------------+------------------------------------------------+



Pipeline output
===============

The major output are fasta files of assembled genomes. Other files are produced (annotations and quality metrics).
is in the database file :file:`csvdb`.



Example
=======

Example data is available at http://www.cgat.org/~andreas/sample_data/pipeline_genome_assemblyxxxxx .
To run the example, simply unpack and untar::

   wget http://www.cgat.org/~andreas/sample_data/pipeline_genome_assembly.tgz
   tar -xvzf pipeline_genome_assembly.tgz
   cd pipeline_genome_assembly
   python <srcdir>/pipeline_genome_assembly.py make full

.. note:: 
   For the pipeline to run, install the :doc:`pipeline_annotations` as well.

Glossary
========

.. glossary::


Code
====

"""

# Import ruffus
from ruffus import *

# Import Python modules
import sys, os, csv, glob, optparse, shutil, gzip, itertools, re, math, types, collections, time 

# Import external modules 
import sqlite3, numpy

import logging as L
from rpy2.robjects import r as R
import rpy2.robjects as ro
import rpy2.robjects.vectors as rovectors
from rpy2.rinterface import RRuntimeError

# Import CGAT tools
import CGAT.IOTools as IOTools
import CGATPipelines.Pipeline as P
import CGAT.Experiment as E
import CGATPipelines.PipelineTracks as PipelineTracks
import CGAT.CSV2DB as CSV2DB
import CGAT.Database as Database
#import PipelineMetagenomeAssembly
import CGAT.CSV as CSV
import CGATPipelines.PipelineMapping as PipelineMapping
import CGATPipelines.PipelineMappingQC as PipelineMappingQC
import CGAT.Stats as Stats
import CGAT.GTF as GTF
import CGAT.IndexedFasta as IndexedFasta


#To do:
#Create PipelineGenomeAssembly utilities script, modify pipeline.ini file, add assemblers, options for each, executables, etc. 
#Continue from mapping pipeline line 1086

#bash compile Cortex
#design calculation for memory, genome complexity, etc.
#readQC
#run case studies
#run available Pseudomonas strains
#Assembly metrics**
#Stitching with flash

# Mapping and variant calling example from Harris et al Nat Biotech 2013	
#Box 1: Methods
#From Read and assembly metrics inconsequential for clinical utility of whole-genome sequencing in mapping outbreaks

#Reads were mapped against the chromosome of an EMRSA-15 reference, HO 50960412 (accession number HE681097), using SMALT version 0.6.4 
#(http://www.sanger.ac.uk/resources/software/smalt/). Reads that aligned equally well to two or more locations in the reference sequence were not mapped. 
#For HiSeq and MiSeq data, paired reads were mapped with an expected insert size of between 50 and 1,000 bp, whereas Ion Torrent data 
#reads were mapped as single ended. In all cases, reads containing indels were realigned using the Genome Analysis Toolkit10 
#RealignerTargetCreator and IndelRealigner modules. Variants were called using a combination of samtools11 mpileup and bcftools. 
#For HiSeq and MiSeq data, default options were used, whereas for Ion Torrent data, settings were chosen as used in the Ion Torrent Variant Caller 
#Plugin: for samtools mpileup, a minimum base quality of 7 was used instead of the default of 13, the coefficient of homopolymer errors was reduced 
#to 50 from the default of 100, the minimum number of reads required for an indel candidate was increased to 4 from a default of 1 and 
#the phred-scaled gap opening and extension probabilities were reduced to 10 and 17 from the defaults of 40 and 20, respectively. 
#Pseudosequences for each isolate were then created by filtering all bases using the following filters: 
#(i) base must be covered by at least four reads, of which at least two must be on each strand; and 
#(ii) bases must have a quality score >50 and a mapping quality >30, and the 
#majority base must be present in at least 75% of reads on each strand. 
#Phylogenetic trees of variable sites were reconstructed with RAxML12, and SNPs reconstructed onto the trees using parsimony. 
#Raw reads for each platform are available from the European Nucleotide Archive under project accession number ERP000985. Accessions for individual samples are shown in Supplementary Table 1.




################################################
##1.0 Genome assembly pipeline ##
################################################


###################################################
###################################################
###################################################
## Pipeline configuration
###################################################

# load options from the config file

P.getParameters(
    ["%s/pipeline.ini" % os.path.splitext(__file__)[0],
     "../pipeline.ini",
     "pipeline.ini" ],
    defaults = {
        'annotations_dir' : "",
        'paired_end' : False } )

PARAMS = P.PARAMS

PARAMS_ANNOTATIONS = P.peekParameters( PARAMS["annotations_dir"],
                                       "pipeline_annotations.py" )


# define input files

#INPUT_FORMATS = ("*.fastq.1.gz", "*.fastq.gz", "*.sra", "*.csfasta.gz")
#REGEX_FORMATS = regex( r"(\S+).(fastq.1.gz|fastq.gz|sra|csfasta.gz)")


###################################################################
###################################################################
## Helper functions mapping tracks to conditions, etc
###################################################################

# collect sra, fasta and and fastq compressed or not tracks
TRACKS = PipelineTracks.Tracks( PipelineTracks.Sample ).loadFromDirectory(
    glob.glob( "*.sra" ), "(\S+).sra" ) +\
    PipelineTracks.Tracks( PipelineTracks.Sample ).loadFromDirectory(
        glob.glob( "*.fastq.gz" ), "(\S+).fastq.gz" ) +\
        PipelineTracks.Tracks( PipelineTracks.Sample ).loadFromDirectory(
            glob.glob( "*.fastq.1.gz" ), "(\S+).fastq.1.gz" ) +\
            PipelineTracks.Tracks( PipelineTracks.Sample ).loadFromDirectory(
                glob.glob( "*.csfasta.gz" ), "(\S+).csfasta.gz" ) +\
                PipelineTracks.Tracks( PipelineTracks.Sample ).loadFromDirectory(
                    glob.glob( "*.fasta.gz" ), "(\S+).fasta.gz" ) +\
                    PipelineTracks.Tracks( PipelineTracks.Sample ).loadFromDirectory(
                        glob.glob( "*.fasta" ), "(\S+).fasta" ) +\
                        PipelineTracks.Tracks( PipelineTracks.Sample ).loadFromDirectory(
                            glob.glob( "*.fastq.gz" ), "(\S+).fastq.gz" ) +\
                            PipelineTracks.Tracks( PipelineTracks.Sample ).loadFromDirectory(
                                glob.glob( "*.fastq" ), "(\S+).fastq" )


# define some tracks if needed
#TRACKS = PipelineTracks.Tracks( PipelineTracks.Sample ).loadFromDirectory( 
#    glob.glob("*.ini" ), "(\S+).ini" )

###################################################################
## Global flags
###################################################################
ASSEMBLERS = P.asList( PARAMS["assemblers"] )

###################################################################
###################################################################
###################################################################
def connect():
    '''connect to database.

    Use this method to connect to additional databases.

    Returns a database connection.
    '''

    dbh = sqlite3.connect( PARAMS["database"] )
    statement = '''ATTACH DATABASE '%s' as annotations''' % (PARAMS["annotations_database"])
    cc = dbh.cursor()
    cc.execute( statement )
    cc.close()

    return dbh

###################################################################
###################################################################
###################################################################
##
###################################################################
if os.path.exists("pipeline_conf.py"):
    L.info( "reading additional configuration from pipeline_conf.py" )
    execfile("pipeline_conf.py")


###################################################################
###################################################################
###################################################################
## Examples from pipeline_quickstart.py for Ruffus tasks
###################################################################

@transform( ["%s.ini" % x.asFile() for x in TRACKS], suffix( ".ini"), ".counts")
def dummyTask( infile, outfile ):
    '''a dummy task - counts number of words in ini files'''

    # can run on cluster
    to_cluster = True

    # the statement
    statement = '''awk 'BEGIN { printf("word\\tfreq\\n"); } 
                        {for (i = 1; i <= NF; i++) freq[$i]++}
                        END { for (word in freq) printf "%%s\\t%%d\\n", word, freq[word] }'
                < %(infile)s > %(outfile)s'''
    P.run()

@transform( dummyTask, suffix(".counts"), "_counts.load" )
def loadDummyTask( infile, outfile ):
    '''load results of word counting into database.'''
    P.load( infile, outfile, "--index=word" )


########################################################################
#########################################################################
#########################################################################
## File types
#########################################################################
SEQUENCEFILES=("*.fastq.1.gz",
               "*.fastq.gz",
               "*.sra",
               "*.export.txt.gz",
               "*.csfasta.gz",
               "*.csfasta.F3.gz",
	       "*.fastq",
	       "*.fasta.gz",
	       "*.fasta",
               )
SEQUENCEFILES_REGEX=regex( r"(\S+).(fastq.1.gz|fastq.gz|sra|csfasta.gz|csfasta.F3.gz|export.txt.gz|fastq|fasta|fasta.gz)")

###################################################################
###################################################################
###################################################################
## Load number of reads
###################################################################
@transform( SEQUENCEFILES,
            SEQUENCEFILES_REGEX,
            r"\1.nreads" )
def countReads( infile, outfile ):
    '''count number of reads in input files.'''
    to_cluster = True
    m = PipelineMapping.Counter()
    statement = m.build( (infile,), outfile )
    P.run()

###################################################
##############1.1 Merge paired reads if necessary (from Nick's pipeline?)
###################################################

@follows( mkdir("merged_reads_Flash.dir" ) )
@transform( SEQUENCEFILES,
            SEQUENCEFILES_REGEX,
            add_inputs( os.path.join( PARAMS["genome"] + ".fa") ),
            r"merged_reads_Flash.dir/\1.flash" )
def flashMergePairedReads(infile, outfile):

    '''# Run Flash to merge paired reads designed with overlapping insert sizes:
# http://ccb.jhu.edu/software/FLASH/MANUAL
# http://ccb.jhu.edu/software/FLASH/FLASH-reprint.pdf
    '''

    #job_options= "-pe dedicated %i -R y" 
    to_cluster = True
    m = PipelineGenomeAssembly.Flash( executable = P.substituteParameters( **locals() )["merge_reads_executable"] )
    infile, reffile = infiles
    #IMS remove reporting options to the ini
    #bowtie_options = "%s --best --strata -a" % PARAMS["bowtie_options"] 
    statement = m.build( (infile,), outfile )
    P.run()

###################################################
##############1.1 Examine quality

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


#############1.2 Run genome assemblers

############
#Run Cortex
############
 
#fastq to vcf - cortex
	#run_calls.pl
		#create index file and sample lists
                #common mistakes so far: missing dots in index file; wrong tabs, structure in index file; different lengths 


@files('*.fastq.1.gz', 'sample-list.tab')
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

    #to_cluster = False
    statement = ''' stampy.py -G %(extra_parameters1)s %(infile)s;
                ''' % locals()

    P.run()

# if using reference build Stampy hash

@follows(stampyIndex)
#@transform('*.fasta', regex(r".*/(.*).fasta"), r"\1.sthash")
@transform('*.fasta', suffix('.fasta'), '.sthash', 'PAO1')
def stampyHash(infile, outfile, extra_parameters1):

    #to_cluster = False
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

    #to_cluster = False
    job_options = '-l h=cgatgpu1'
    statement = ''' cortex_var_31_c96 
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

    #to_cluster = False
    job_options = '-l h=cgatgpu1'
    statement = ''' cortex_var_63_c96 
                    --kmer_size 61 
                    --mem_height 17 
                    --mem_width 100 
                    --se_list %(infile)s
                    --dump_binary %(outfile)s 
                    --sample_id REF
                ''' % locals()

    P.run()
 
@follows(buildCortexKmer_2)
@transform('*.selist', suffix('.selist'), '.k127.ctx')
def buildCortexKmer_3(infile, outfile):

    #to_cluster = False
    job_options = '-l h=cgatgpu1'
    statement = ''' cortex_var_127_c96
                    --kmer_size 127 
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

@follows(buildCortexKmer_3)
@files('index.tab', 'run_calls.dir', 'PAO1_phase1.', 'log_QT_10.run_calls.sentinel')
def CortexRunCalls(infile, outfile, extra_parameters1, extra_parameters2):

    #job_options ='-l mem_free=130G'
    job_options = '-l h=cgatgpu1'
    #to_cluster = False
    statement = ''' run_calls.pl 
                    --first_kmer 31 
                    --last_kmer 61 
                    --kmer_step 30 
                    --fastaq_index %(infile)s 
                    --auto_cleaning yes 
                    --bc yes 
                    --pd yes 
                    --outdir %(outfile)s 
                    --outvcf %(extra_parameters1)s 
                    --ploidy 1 
                    --stampy_hash PAO1 
                    --stampy_bin /ifs/apps/bio/stampy-1.0.17/stampy.py 
                    --list_ref_fasta PAO1_nolines.selist 
                    --refbindir ./ 
                    --genome_size 6300000 
                    --mem_height 22 
                    --mem_width 200 
                    --vcftools_dir /ifs/apps/src/vcftools_0.1.11/vcftools_0.1.11/ 
                    --do_union yes 
                    --logfile %(extra_parameters2)s 
                    --workflow independent 
		    --ref CoordinatesAndInCalling 
                    --squeeze_mem 
                    --apply_pop_classifier
		    --qthresh 10
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


#Merge separate paired end files into one using velveth tool

#@transform('*.fastq.gz', suffix('.fastq.gz'), '.VelvetOptimiser.dir/')

#def shuffleFastq(infile, outfile): 
#    to_cluster = True
#    statement = ''' /ifs/apps/bio/velvet-1.2.08/bin/contrib/shuffleSequences_fasta/shuffleSequences_fastq.pl
#                    %(infile)s' 
#                    %(outfile)s
#                ''' % locals()
#    P.run()



#VelvetOptimiser.pl -s 155 -e 301 -f '-fastq -long S1-1-1.flash.extendedFrags.fastq' -o '-min_contig_lgth 200 -clean yes -read_trkg yes -amos_file yes -exp_cov auto -cov_cutoff auto' -d S1-1-1.flash.extendedFrags.fastq.VelvetOptimiser.kmer155to301.dir/ -t 10

@transform('*.fastq.gz', suffix('.gz'), '.VelvetOptimiser.dir/')

def velvetOptimiser(infile, outfile): 
    
    job_options ='-l mem_free=30G' 
    to_cluster = False 
    statement = ''' VelvetOptimiser.pl 
			-s 155 
			-e 301 
			-f '-fastq -long %(infile)s' 
			-o '-min_contig_lgth 200 -clean yes -read_trkg yes -amos_file yes -exp_cov auto -cov_cutoff auto' 
			-d %(outfile)s
			-t 12   
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


@follows(mkdir('velveth.dir'))
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



#######################
#Run SOAPdenovo
#######################

#http://soap.genomics.org.cn/soapdenovo.html
#${bin} all -s config_file -K 63 -R -o graph_prefix 1>ass.log 2>ass.err

@follows(mkdir('SOAPdenovo.dir'))
@transform('SOAPdenovo_config_file.txt', suffix('.txt'), '.SOAPdenovo')

def runSOAPdenovo(infile, outfile): 
    
    job_options ='-l mem_free=130G' 
    to_cluster = True
    statement = ''' SOAPdenovo-63mer all
		    -s %(infile)s
		    -K 63
		    -R
		    -o %(outfile)s
		    1>ass.log 2>ass.err
		    -p 10
		    -V 
		    -N 6300000 
                ''' % locals()
    P.run()


#######################
#Run Ray
#######################


#@active_if("ray" in ASSEMBLERS)
#@follows(mkdir("ray.dir"))
#@transform( SEQUENCE_TARGETS[PARAMS["pool_reads"]],
#            SEQUENCEFILES_REGEX
#            , r"ray.dir/\1.contigs.fa")
#def runRay(infile, outfile):
#    '''
#    run Ray on each track
#    '''
#    to_cluster = False
##    job_options=" -l h=!cgatsmp1,h=!cgatgpu1,h=!andromeda,h=!gandalf,h=!saruman -pe mpi 10 -q all.q "
#    job_options = "-l mem_free=30G"
#    statement = PipelineMetagenomeAssembly.Ray().build(infile)
#    P.run()


#Gunzip first as Ray requires LIBZ compiled at installation to run fastq.gz files

@transform('*.gz', suffix('.gz'), '.fastq')

def gunzip(infile, outfile): 
    
    to_cluster = True
    statement = ''' gunzip %(infile)s
                ''' % locals()
    P.run()


#@follows(gunzip, mkdir('Ray.dir'))
@transform('*.fastq', suffix('.fastq'), '.Ray.dir/')

def runRay(infile, outfile): 
    
   # job_options ='-l h=!andromeda,h=!cgatgpu1,h=!cgatsmp1,h=!gandalf,h=!saruman -l mem_free=130G -pe mpi 10 -q all.q' 
    job_options = '-pe mpi 10 -q mpi.q'
    to_cluster = True
    statement = ''' Ray
		    -k 63
		    -s %(infile)s
		    -o %(outfile)s
		    --minimum-contig-length 1000
                ''' % locals()
    P.run()


#######################
#Run SGA
#######################


@transform('*.fastq.gz', suffix('.fastq.gz'), '.SGA.preprocess.fastq')
def runSGApreprocess(infile, outfile): 
# Preprocess the data to remove ambiguous basecalls
# Example: sga preprocess -m -q -o S1-1-1.flash.extendedFrags.SGA.preprocess S1-1-1.flash.extendedFrags.fastq.gz 
    
    to_cluster = True
    statement = ''' sga preprocess
		    -m
		    -q
		    -o %(outfile)s
		    %(infile)s
                ''' % locals()
    P.run()


@follows(runSGApreprocess)
@transform('*.preprocess.fastq', suffix('.preprocess.fastq'), '.preprocess.fastq.index.sentinel')
def runSGAindex(infile, outfile): 
# FM-index the preprocessing data
# Example: sga index -v -a sais -t 2 -c S1-1-1.flash.extendedFrags.SGA.preprocess
# -v option adds information to the JSON file that later causes an error (ie no longer a JSON object)     
    to_cluster = True
    statement = ''' sga index
		    -a sais
		    -t 4
		    -c
		    %(infile)s
		    ;
		    touch %(outfile)s
                ''' % locals()
    P.run()


@follows(runSGAindex)
@transform('*.preprocess.fastq', suffix('.preprocess.fastq'), '.preqc', 'PAO1_nolines.fa')
def runSGApreqc(infile, outfile, extra_parameters1): 
# Run quality control checks on preprocessed data to fine tune assembling
# Example: sga preqc --reference=PAO1.fasta -t 2 S1-1-1.flash.extendedFrags.fastq.sample.SGA.preprocess > S1-1-1.flash.extendedFrags.fastq.sample.SGA.no-v.preqc
   
    to_cluster = True
    statement = ''' sga preqc 
		    --reference=%(extra_parameters1)s
		    -t 4
		    %(infile)s
		    > %(outfile)s
                ''' % locals()
    P.run()


@follows(runSGApreqc)
@transform('*.preqc', suffix('.preqc'), '.preqc-report.pdf')
def runSGApreqcReport(infile, outfile): 
# Generate PDF report for QC. To do: pull out for SphinxReport
# Example: sga-preqc-report.py S1-1-1.flash.extendedFrags.fastq.sample.SGA.copy.preqc 

    to_cluster = True
    statement = ''' sga-preqc-report.py
		    %(infile)s
		    -o %(outfile)s
                ''' % locals()
    P.run()

@follows(runSGApreqcReport)
@transform('*.preprocess.fastq', regex(r'(.*).preprocess.fastq'), r'\1.corrected.fastq', r'\1.error-correction-output')
def runSGAcorrect(infile, outfile, extra_parameters1): 
# Error correct and build the index that will be used for error correction. Error correction is k-mer based. --learn option allows k-mer cutoff parameter to be learned automatically.    
# Example: sga correct -a kmer -t 4 --discard --metrics=S1-1-1.flash.extendedFrags.error-correction-output --learn -o S1-1-1.flash.extendedFrags.SGA.corrected S1-1-1.flash.extendedFrags.SGA.preprocess &
# Need to summarise and plot error correction output, save stderr for reads discarded, corrected, etc. and use error-correction-output files for histograms


    to_cluster = True
    statement = ''' sga correct 
		    -a kmer 
		    -t 4
		    --discard
		    --learn
		    --metrics=%(extra_parameters1)s 
		    -o %(outfile)s
		    %(infile)s
		''' % locals()
    P.run()

@follows(runSGAcorrect)
@transform('*.corrected.fastq', suffix('.corrected.fastq'), '.corrected.fastq.index.sentinel')
def runSGAindexCorrected(infile, outfile): 
# Run FM-index on corrected basecalls fastq files for contig assembly
    
    to_cluster = True
    statement = ''' sga index
		    -a sais
		    -t 4
		    -c %(infile)s
		    ;
		    touch %(outfile)s
                ''' % locals()

    P.run()
  


@follows(runSGAindexCorrected)
@transform('*.corrected.fastq', suffix('.corrected.fastq'), '.filter.pass.fasta')
def runSGAfilter(infile, outfile): 
# Remove exact-match duplicates and reads with low-frequency k-mers
# Example: sga filter -t 4 --homopolymer-check --low-complexity-check S1-1-1.flash.extendedFrags.SGA.corrected
   
    to_cluster = True
    statement = ''' sga filter
		    --homopolymer-check
		    --low-complexity-check
		    -t 4
		    -o %(outfile)s
		    %(infile)s
		 ''' % locals()

    P.run()


@follows(runSGAfilter)
@transform('*.filter.pass.fasta', suffix('.fasta'), '.fasta.overlap.sentinel')
def runSGAoverlap(infile, outfile): 
# Compute the structure of the string graph
# Example: sga overlap -t 4 reads.ec.filter.pass.fa
   
    to_cluster = True
    statement = ''' sga overlap
		    -t 4
		    %(infile)s
		    ;
		    touch %(outfile)s
		 ''' % locals()

    P.run()

@follows(runSGAoverlap)
@transform('*.filter.pass.asqg.gz', regex(r'(.*).asqg.gz'), r'\1.assemble.sentinel', r'\1.assemble')
def runSGAassemble(infile, outfile, extra_parameters1): 
# Perform the contig assembly
# Example: sga assemble -m 111 --min-branch-length 400 -o primary reads.ec.filter.pass.asqg.gz

    to_cluster = True
    statement = ''' sga assemble
		    -o %(extra_parameters1)s 
		    %(infile)s
		    ;
		    touch %(outfile)s
		 ''' % locals()

    P.run()

@follows(runSGAassemble)
@transform('*contigs.fa', suffix('contigs.fa'), 'contigs.fa.bwa-index.sentinel')
def runBWAindex(infile, outfile):
# Scaffolding: align the reads to the contigs
# Use bwa and samtools to index and align: bwa index, bwa aln, bwa sampe/samse | samtools view > *.bam 

    to_cluster = True
    statement = ''' bwa index %(infile)s
		    ;
		    touch %(outfile)s
		''' % locals()

    P.run()

@follows(runBWAindex)
@transform('*.SGA.filter.pass.assemble-contigs.fa', regex(r'(.*).SGA.filter.pass.assemble-contigs.fa'), r'\1.SGA.filter.pass.assemble.bwa-aln.sai', r'\1.fastq.gz')
def runBWAalign(infiles, outfile, extra_parameters1): 
   
    to_cluster = True
    statement = ''' bwa aln -t 4 %(infiles)s %(extra_parameters1)s  
		    >
		    %(outfile)s 
		''' % locals()

    P.run()

@follows(runBWAalign)
@transform('*.SGA.filter.pass.assemble-contigs.fa', regex(r'(.*).SGA.filter.pass.assemble-contigs.fa'), r'\1.SGA.filter.pass.assemble.bwa-aln.bam', r'\1.SGA.filter.pass.assemble.bwa-aln.sai', r'\1.fastq.gz')
def generateBAMs(infile, outfile, extra_parameters1, extra_parameters2): 
    
    to_cluster = True
    statement = ''' bwa samse %(infile)s %(extra_parameters1)s %(extra_parameters2)s  
		    | samtools view -Sb - > %(outfile)s 
		''' % locals()

    P.run()

@follows(generateBAMs)
@transform('*.SGA.filter.pass.assemble.bwa-aln.bam', regex(r'(.*).SGA.filter.pass.assemble.bwa-aln.bam'), r'\1.SGA.filter.pass.assemble.bwa-aln.contig-distances.sentinel', r'\1.SGA.filter.pass.assemble.bwa-aln-bam')
def estimateContigDistances(infile, outfile, extra_parameters1): 
# Convert BAM files into a set of contig-contig distance estimates
# sga-bam2de.pl -n 10 -m 200 --prefix libPE libPE.bam

   
    to_cluster = True
    statement = ''' sga-bam2de.pl
		    -n 10
		    -m 200
		    -t 4
		    --prefix %(extra_parameters1)s 
		    %(infile)s
		    ;
		    touch %(outfile)s 
		''' % locals()

    P.run()

@follows(estimateContigDistances)
@transform('*.SGA.filter.pass.assemble.bwa-aln.bam', regex(r'(.*).SGA.filter.pass.assemble.bwa-aln.bam'), r'\1.SGA.filter.pass.assemble.bwa-aln.contig-copies.astat')
def estimateContigCopies(infile, outfile): 
# Compute copy number estimates of the contigs
# sga-astat -m 200 libPE.bam > libPE.astat
   
    to_cluster = True
    statement = ''' sga-astat.py
		    -m 200
		    %(infile)s
		    >
		    %(outfile)s 
		''' % locals()

    P.run()

@follows(estimateContigCopies)
@transform('*.SGA.filter.pass.assemble-contig.fa', regex(r'(.*).SGA.filter.pass.assemble.bwa-aln.bam'), r'\1.SGA.filter.pass.assemble.scaffold.sentinel', r'\1.SGA.filter.pass.assemble.bwa-aln.contig-copies.astat', r'\1.SGA.scaffold.scaf')
def buildScaffolds(infile, outfile, extra_parameters1, extra_parameters2): 
# Build the scaffolds
# sga scaffold -m 200 -a libPE.astat -o scaffolds.scaf --pe libPE.de $PRIMARY_CONTIGS
# Use --strict to build very conservative scaffolds that should be highly accurate. 
   
    to_cluster = True
    statement = ''' sga scaffold
		    -m 200
		    -a %(extra_parameters1)s
		    --remove-conflicting
		    -o %(extra_parameters2)s
		    ;
		    touch %(outfile)s
		''' % locals()

    P.run()


@follows(buildScaffolds)
@transform('*.SGA.filter.pass.assemble-contig.fa', regex(r'(.*).SGA.filter.pass.assemble.bwa-aln.bam'), r'\1.SGA.filter.pass.assemble.scaffold-fasta.sentinel', r'\1.SGA.filter.pass.asqg.gz', r'\1.SGA.scaffold.fa')
def runScaffoldToFasta(infile, outfile, extra_parameters1, extra_parameters2): 
# Convert the scaffolds to FASTA format
# sga scaffold2fasta --use-overlap --write-unplaced -m 200 -a $PRIMARY_GRAPH -o sga-scaffolds.fa scaffolds.scaf
  
    to_cluster = True
    statement = ''' sga scaffold2fasta
		    -m 200
		    --use-overlap
		    --write-unplaced
		    -a %(extra_parameters1)s
		    -o %(extra_parameters2)s
		    ;
		    touch %(outfile)s
		''' % locals()

    P.run()


#######################
#Run SPAdes
#######################

#--careful option requires paired end files to run, outputs with error but generates files (not contigs or scaffolds though)
#RE-RUN with --careful option? 

@mkdir('SPAdes.dir')
@transform('*.1.fastq.gz', regex(r'(.*).(.).fastq.gz'), r'\1.SPAdes', r'\1.2.fastq.gz')
def runSPAdes(infile, outfile, extra_parameters1): 
    
    job_options ='-l mem_free=30G' 
    to_cluster = True
    statement = ''' cd SPAdes.dir ;
		    spades.py 
		    -o %(outfile)s 
		    -1 ../%(infile)s
		    -2 ../%(extra_parameters1)s
		    --careful ;
		    mv %(outfile)s/contigs.fasta %(outfile)s/%(outfile)s.fasta ;
		    cd .. ;
		    touch %(outfile)s 
                ''' % locals()
    P.run()



#############1.5 Assess quality of assemblies

#Run QUAST and other tools to determine various quality metrics of each assembly

#######################
#Run QUAST
#######################
#quast.py -R PAO1_nolines.fa -L --gene-finding -t 5 --strict-NA S1-1-1.flash.extendedFrags.SPAdes.dir/contigs.fasta


@follows(runSPAdes)
@files(None, 'QUAST_post_Assembly', 'PAO1_nolines.fa')
def runQUASTpostAssembly(infile, outfile, extra_parameters1): 
    to_cluster = True
    statement = ''' ln -s SPAdes.dir/trimmomatic_S*/*SPAdes.fasta . ;
		    rename trimmomatic 2 *SPAdes.fasta ;
		    quast.py
                    -R %(extra_parameters1)s
                    -t 5
		    --strict-NA
		    --gene-finding
		    --gage
		    -o %(outfile)s *SPAdes.fasta  
                ''' % locals()
    P.run()



#######################
#Run Mauve Assembly Metrics
#######################


#Mauve Assembly Metrics: Statistical View of the Contigs:

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

#@follows(velvetOptimiser, mkdir('Mauve_Assembly_Metrics.dir'))
@follows(mkdir('Mauve_Assembly_Metrics.dir'))
@transform('*.dir/*.fasta', suffix('.fasta'), '.Mauve_Assembly_Metrics.dir/')

def runMauveAssemblyMetrics(infile, outfile): 
    to_cluster = True
    statement = ''' mauveAssemblyMetrics.sh
                    -reference PAO1.fasta 
                    -assembly %(infile)s
                    -reorder
                    -outputDir %(outfile)s
                    ; checkpoint
                ''' % locals()
    P.run()

#                    ; mauveAssemblyMetrics .


#######################
#Run contigs2stats 
#######################


@transform('*.dir/contigs.fasta', suffix('.fasta'), '.contigs2stats')

def runContigs2Stats(infile, outfile):
    to_Cluster = True
    statement = ''' python /ifs/devel/antoniob/cgat/scripts/contigs2stats.py
		    -n 50 
		    -I %(infile)s
		    -S %(outfile)s
		''' % locals()

    P.run()




#############1.6 Post-assembly improvement
#Quality assess assemblies produced (QUAST, done above), detect and break misassemblies (REAPR), align, gap fill, error correct and annotate with PAGIT tools. 
 
#Run PAGIT tools to align the genome (ABACAS), fill gaps (IMAGE), correct errors (ICORN) and transfer annotations from a reference (RATT).
#Repeat assembly quality assessment after each step (run QUAST). 
#Transfer notes on Mac (/CGAT/Pseudomonas.../PAGIT...)


#######################
#Run REAPR for detecting and breaking misassemblies
#######################

#Manual:ftp://ftp.sanger.ac.uk/pub/resources/software/reapr/Reapr_1.0.17.manual.pdf
#REAPR uses fastq reads to error correct the assemblies. Ideally long insert size reads together with high quality short insert size reads are provided but it can be run with one library (typically short insert Illumina data). 
#Output files can be used for Abacas (contig ordering and alignment for post-assembly improvement) and for plotting. 

#Get fasta files from best assembler:
#ln -s ../genome_assembly_tests_flash_SPAdes-3.0_1.dir/S1-1-1.flash.extendedFrags.SPAdes.dir/contigs.fasta .

#Get paired fastq data from original files:
#ln -s /ifs/projects/proj018/runs/readQC_Flash.dir/*fastq .

#Run facheck if necessary (if already providing bam files of paired reads mapped to assembly), use *.info file to rename with seqrename command afterwards if there are errors thrown during facheck:
#reapr facheck contigs.fasta S1-facheck-output

#Rename sequences in fasta file if necessary with:
#reapr seqrename

#Run perfectmap (optional) to help score each base of the genome. Requires short insert Illumina reads and assumes that all reads are of the same length:
#reapr perfectmap assembly.fa reads_1.fq reads_2.fq <insert_size> perfect
#reapr perfectmap S1-1-1.SPAdes.contigs.fasta S1-1-1.1.fastq S1-1-1.2.fastq S1-1-1.REAPR.smaltmap.bam

#If not providing bam files, map paired reads with smaltmap command *better as maps, sorts and removes duplicates):
#reapr smaltmap S1-1-1.SPAdes.contigs.fasta S1-1-1.1.fastq S1-1-1.2.fastq S1-1-1.REAPR.smaltmap.bam

#Run perfectmap if useful, this requires both short and long paired reads to map back to assembly, otherwise go to pipeline command after smaltmap:
#reapr perfectmap S1-1-1.SPAdes.contigs.fasta S1-1-1.1.fastq S1-1-1.2.fastq 385 S1-1-1.SPAdes.REAPR.perfectmap.

#Run reapr pipeline to break and correct assembly:
#reapr pipeline S1-1-1.SPAdes.contigs.fasta S1-1-1.REAPR.smaltmap.bam S1-REAPR-pipeline.dir

#Tests for REAPR run in /ifs/projects/proj018/runs/PAGIT_runs.dir/ for S1 samples with SGA and SPAdes assemblies.
#Ran with reapr smaltmap, reapr perfectmap and reapr pipeline.
#Results are similar with the error "Low fragment coverage within a contig: 2766" called and most bases converted to 'N'. For this case run:
	
	#reapr pipeline break -b=1 to ignore fragments with low coverage without a gap

#TO DO: run facheck on fasta files to make sure sequence names won't break the pipeline.

#Output and report check for perfectbam (per base coverage):
	#out.perfect_cov.gz (chromosomes with zero coverage are not present in this file)
	#out.hist (anything over 100 is reported as 100 in this file)

#Sanity checks for REAPR: check 00.Sample directories and look at GC coverage and insert size plots: 
	#gcvscov.lowess.pdf ; 
	#the fragment size distribution: insert.in.pdf , insert.out.pdf , insert.same.pdf . The in plot should be your expected insert size distribution and the other two noise.
	#Check insert.stats.txt 
	#In the wd check 05.summary.report.txt and 03.score.errors.gff.gz for summary stats, errors and warnings.

#Output files for ACT, Abacas and downstream are:
	#04.break.broken_assembly.fa (result of breaking contigs at gaps and turning indels to N's)
	#04.break.broken_assembly_bin.fa (contains the actual sequences of what was turned into an N)
	#Run reaper plot to get Artemis files for specific contigs, these aren't run by default
	#

#reapr perfectmap S1-1-1.SPAdes.contigs.fasta S1-1-1.1.fastq S1-1-1.2.fastq 385 output_prefix

#Rename fastqs for IMAGE later on (it adds an underscore to the fastq prefix provided and messes the files called giving silent errors)
def renameFASTQsForIMAGE():
    
    to_cluster = False
    statement = '''rename 1.fastq paired_1.fastq *fastq ;
                   rename 2.fastq paired_2.fastq *fastq 
      		'''
    P.run()

#@follows(renameFASTQsForIMAGE)
@transform('*.fasta', regex('(\w+[\-]\w+[\-]\w+).(?:.*).fasta'), r'\1.REAPR.perfectmap', r'\1.paired_1.fastq', r'\1.paired_2.fastq', '385')
def runPerfectMap(infile, outfile, extra_parameters1, extra_parameters2, extra_parameters3): 
    
    job_options ='-l mem_free=30G' 
    to_cluster = True
    statement = ''' reapr perfectmap %(infile)s %(extra_parameters1)s %(extra_parameters2)s %(extra_parameters3)s %(outfile)s ; touch %(outfile)s
                ''' % locals()
    P.run()


@follows(runPerfectMap)
@transform('*.fasta', regex('(\w+[\-]\w+[\-]\w+).(?:.*).fasta'), r'\1.REAPR.smaltmap.bam', r'\1.paired_1.fastq', r'\1.paired_2.fastq')
def runSmaltMap(infile, outfile, extra_parameters1, extra_parameters2):
    
    job_options ='-l mem_free=30G' 
    to_cluster = True
    statement = ''' reapr smaltmap %(infile)s %(extra_parameters1)s %(extra_parameters2)s %(outfile)s
                ''' % locals()

    P.run()

#If running reapr pipeline with perfectmap, add the perfectmap output_prefix at the end of the command
@follows(runSmaltMap)
@transform('*.fasta', regex('(\w+[\-]\w+[\-]\w+).(?:.*).fasta'), r'\1.REAPR.dir', r'\1.REAPR.smaltmap.bam', r'\1.REAPR', r'\1.REAPR.perfectmap')
def runREAPRpipeline(infile, outdir, extra_parameters1, extra_parameters2, extra_parameters3): 
    
    job_options ='-l mem_free=30G' 
    to_cluster = True
    statement = ''' reapr pipeline -break b=1 %(infile)s %(extra_parameters1)s %(outdir)s %(extra_parameters3)s > %(outdir)s.pipeline.log ;
		    mv %(extra_parameters2)s.dir.* %(extra_parameters2)s.perfectmap.* %(extra_parameters1)s %(extra_parameters1)s.bai %(outdir)s ; touch %(extra_parameters1)s 
                ''' % locals()
    P.run()


@follows(runREAPRpipeline)
#@posttask(touch_file('QUASTpostREAPR_completed.sentinel'))
@files(None, 'QUAST_post_REAPR', 'PAO1_nolines.fa')
def runQUASTpostREAPR(infile, outfile, extra_parameters1): 
    to_cluster = True
    statement = ''' quast.py
                    -R %(extra_parameters1)s
                    -t 5
		    --strict-NA
		    --gene-finding
		    -L
		    --gage
		    -o %(outfile)s
                    *.REAPR.dir/04.break.broken_assembly.fa 
                ''' % locals()
    P.run()

#######################
#Setup PAGIT
#######################
#Setup directories for ABACAS, IMAGE, ICORN, RATT
#Place reference embl annotation file and .fa and .fa.fai files in working directory
#Source PAGIT environment

@follows(runQUASTpostREAPR)
def sourcePAGIT():
    to_cluster = False
    statement = '''source /ifs/apps/bio/PAGIT-v1/PAGIT/sourceme.pagit
		'''
    P.run()


#######################
#Run ABACAS
#######################
#abacas.pl 

#perl $PAGIT_HOME/ABACAS/abacas.pl -r Refsequence.fasta -q assembly.fasta -p nucmer -m -b -t -o myPrefix
#mkdir %(extra_parameters2)s ; 
#cd %(extra_parameters2) ;
#-m	print ordered contigs to file in multifasta format 
#-b	print contigs in bin to file 
#-l	int	minimum contig length [default 100] Contigs below this are not used.
#-t 	run tblastx on contigs that are not mapped (much slower)
#abacas manual is at http://abacas.sourceforge.net/Manual.html#5._Options

@follows(sourcePAGIT)
@mkdir('ABACAS.dir')
@mkdir('*.paired_*', formatter('.+/(?P<basename>.*).paired_(.).fastq'), ['{path[0]}/ABACAS.dir/{basename[0]}.ABACAS'])
#Create blast format database for reference genome in order to blast with ABACAS
@originate('formatdb.completed')
def formatDB(outfile):
    statement = ''' formatdb -p F -i PAO1_nolines.fa ; touch %(outfile)s
		''' % locals()
    P.run()

'''
@follows(formatDB)
@mkdir('*.REAPR.dir', formatter('.+/(?P<basename>.*).REAPR.dir'), ['{path[0]}/ABACAS.dir/{basename[0]}.ABACAS'])
@transform('/ifs/projects/proj018/runs/PAGIT_runs.dir/PAGIT_on_Ruffus_3.dir/*.REAPR.dir/04.break.broken_assembly.fa', regex('(.*).REAPR.dir/04.break.broken_assembly.fa'), r'\1.ABACAS', '../../PAO1_nolines.fa', r'ABACAS.dir/\1.ABACAS/')
def runABACAS2(infile, outfile, extra_parameters1, extra_parameters2): 

   # job_options ='-l mem_free=30G' 
    to_cluster = True
    statement = '' cd %(extra_parameters2)s ;
		    abacas 
		    -r %(extra_parameters1)s 
		    -q %(infile)s 
		    -p nucmer
		    -m -b -d -l 500 
		    -o %(outfile)s > %(outfile)s.log ; 
		    cat %(outfile)s.fasta %(outfile)s.contigsInbin.fas > %(outfile)s.mappedAndUnmapped.fasta ; 
		    touch %(outfile)s
		''
#mv blast.out %(extra_parameters2)s if using -t option
    P.run()
'''

@follows(formatDB)
@jobs_limit(1)
#@transform('*.SPAdes.fasta', suffix('.SPAdes.fasta'), r'ABACAS.dir/\1.ABACAS/\1', 'PAO1_nolines.fa', r'ABACAS.dir/\1.ABACAS/')
@transform('*.REAPR.dir/04.break.broken_assembly.fa', suffix('.REAPR.dir/04.break.broken_assembly.fa'), r'ABACAS.dir/\1.ABACAS/\1', 'PAO1_nolines.fa', r'ABACAS.dir/\1.ABACAS/')
def runABACAS(infile, outfile, extra_parameters1, extra_parameters2): 

   # job_options ='-l mem_free=30G' 
    to_cluster = True
    statement = ''' abacas 
		    -r %(extra_parameters1)s 
		    -q %(infile)s 
		    -p nucmer
		    -m -b -d -l 500 
		    -o %(outfile)s > %(outfile)s.log ; 
		    mv nucmer* unused_contigs.out %(extra_parameters2)s ;
		    cat %(outfile)s.fasta %(outfile)s.contigsInbin.fas > %(outfile)s.mappedAndUnmapped.fasta ; 
		    touch %(outfile)s
		'''
#mv blast.out %(extra_parameters2)s if using -t option
    P.run()

#######################
#Run IMAGE
#######################

#No restart, quick option is (as most gaps get closed in the first 2 or 3 iterations): perl $PAGIT_HOME/IMAGE/image.pl -scaffolds inputScaffolds.fasta -prefix pairedReadsPart -iteration 1 -all_iteration 9 -dir_prefix ite -kmer 55
#nohup perl $PAGIT_HOME/IMAGE/image.pl -scaffolds S1-1-1.SGA.ABACAS.mappedAndUnmapped.fasta -prefix ../../S1-1-1.paired -iteration 1 -all_iteration 5 -dir_prefix ite -kmer 55 &

#**IMAGE adds an underscore just before the pair number for the -prefix option (eg _1.fastq), if it can't find the file it gives a silent (in stdout but doesn t stop iterations and provides output files) 
#error when calling smalt during its pipeline.
#ie paired fastq files need to be renamed to eg *_1.fastq: rename .fastq.1 .paired_1.fastq *1 rename .fastq.2 .paired_2.fastq *2

#Use restartIMAGE.pl to optimise gap closing, eg restartIMAGE.pl last_iteration_folder kmer number_of_iterations Illumina_lane_prefix 
#files have to be at the top directory
#nohup image.pl -scaffolds S1-1-1.mappedAndUnmapped.fasta -prefix S1-1-1.paired -iteration 1 -all_iteration 5 -dir_prefix S1-1-1.ite & works with everything in one directory
#nohup restartIMAGE.pl S1-1-1.ite5 111 2 partitioned & works in the same directory as above (/ifs/projects/proj018/runs/PAGIT_runs.dir/PAGIT_on_Ruffus_single_file_tests.dir/IMAGE_tests.dir)

#Summary files are: 
	#new.fa = updated contigs created during gap-closing
	#new.read.placed = maps of contigs to saffolds for the current iteration
	#walk2.summary = short description of gap-closing results after each iteration, re-run if still extending or closing gaps
	#'paritioned' = fastq files with subset of reads that only contain reads involved in spanning gaps (reduces execution time)
	#iteration summary stats from restartIMAGE.pl script = summary of gap-closing results from all iterations

#IMAGE generates several files in the base directory (delete if re-using image.pl script in the same base directory):
	#contigs.fa@
	#contigs.fa.new
	#contigs.fa.new.list
	#contigs.fa.original@
	#insert_size
	#partitioned_1.fastq@
	#partitioned_2.fastq@
	#read.placed@
	#read.placed.new
	#read.placed.original@
	#submit.S1-1-1.ite.log
	#tmp.fa

#image.pl also changes the abacas mappedAndUnmapped files to:
	#S1-1-1.mappedAndUnmapped.fasta.bak@
	#S1-1-1.mappedAndUnmapped.fasta.bak.above300.fasta
	#S1-1-1.mappedAndUnmapped.fasta.bak.above300.read.placed

@follows(runABACAS)
@mkdir('IMAGE.dir')
#@files(None, 'linFASTQs.sentinel') 
#def linkFASTQs(infile, outfile):
    
#    to_cluster = False
#    statement = '''cd IMAGE.dir/ ;
#		   ln -s /ifs/projects/proj018/backup/gzip_temp/*fastq* . ;
#		   rename .fastq.1 .paired_1.fastq *1 ;
#                   rename .fastq.2 .paired_2.fastq *2 ; 
#		   touch %(outfile)s  
#      		'''
#    P.run()


#@follows(renameFastqsForIMAGE)
#def IMAGEabsPaths():
#    IMAGEinput = ''.join(os.path.abspath(os.getcwd())+'/ABACAS.dir/S*/*.mappedAndUnmapped.fasta')
    #IMAGEregex = ''.join(os.path.abspath(os.getcwd())+'/ABACAS.dir/(.*).ABACAS/(.*).mappedAndUnmapped.fasta')
    #IMAGEextra1 = ''.join(IMAGEregex+r'\1.paired')

#@follows(IMAGEabsPaths)
#@follows(linkFASTQs)
@mkdir('*.paired_*', formatter('.+/(?P<basename>.*).paired_(.).fastq'), ['{path[0]}/IMAGE.dir/{basename[0]}.IMAGE'])
#@transform([''.join(os.path.abspath(os.getcwd())+'/ABACAS.dir/S*/*.mappedAndUnmapped.fasta')], regex(r"''.join(os.path.abspath(os.getcwd())+'/ABACAS.dir/(.*).ABACAS/(.*).mappedAndUnmapped.fasta')"), r'\1.ite', r"''.join(os.path.abspath(os.getcwd())+'\1.paired')", r'IMAGE.dir/\1.IMAGE')
#@transform([os.path.join((os.getcwd()),'/ABACAS.dir/S*/*.mappedAndUnmapped.fasta')], regex(os.path.join(os.getcwd(),r'/ABACAS.dir/(.*).ABACAS/(.*).mappedAndUnmapped.fasta')), r'\1.ite', r'os.path.join(os.getcwd(), r'\1.paired')', r'IMAGE.dir/\1.IMAGE')
#@transform(IMAGEabsPaths, regex(r'(.*).mappedAndUnmapped.fasta'), r'\1.ite', os.path.join(os.getcwd(), r'\1.paired'), r'IMAGE.dir/\1.IMAGE')
#@transform('/ifs/projects/proj018/runs/PAGIT_runs.dir/PAGIT_on_Ruffus_single_file_tests.dir/tests_2.dir/ABACAS.dir/S*/*.mappedAndUnmapped.fasta', regex(r'/ifs/projects/proj018/runs/PAGIT_runs.dir/PAGIT_on_Ruffus_single_file_tests.dir/tests_2.dir/ABACAS.dir/(.*).ABACAS/(.*).mappedAndUnmapped.fasta'), r'\1.ite', r'/ifs/projects/proj018/runs/PAGIT_runs.dir/PAGIT_on_Ruffus_single_file_tests.dir/tests_2.dir/\1.paired', r'IMAGE.dir/\1.IMAGE')
@jobs_limit(1)
@transform('/ifs/projects/proj018/runs/PAGIT_runs.dir/PAGIT_on_Ruffus_all_3.dir/ABACAS.dir/S*/*.mappedAndUnmapped.fasta', regex(r'/ifs/projects/proj018/runs/PAGIT_runs.dir/PAGIT_on_Ruffus_all_3.dir/ABACAS.dir/(.*).ABACAS/(.*).mappedAndUnmapped.fasta'), r'\1.ite', r'\1.paired', r'IMAGE.dir/\1.IMAGE', r'\1')
def runIMAGE(infile, outdir, extra_parameters1, extra_parameters2, extra_parameters3): 

    #job_options ='-l mem_free=30G' 
    to_cluster = True
    statement = '''cd %(extra_parameters2)s ;
		   ln -s /ifs/projects/proj018/backup/gzip_temp/%(extra_parameters3)s.fastq* . ;
		   rename .fastq.1 .paired_1.fastq *1 ;
		   rename .fastq.2 .paired_2.fastq *2 ;
		   image.pl 
		   -scaffolds %(infile)s 
		   -prefix %(extra_parameters1)s 
		   -iteration 1
		   -all_iteration 3 
		   -dir_prefix %(outdir)s
		   -kmer 31 > %(outdir)s.log ;
		   restartIMAGE.pl %(outdir)s3 61 3 partitioned > %(outdir)s.6.log ;
		   restartIMAGE.pl %(outdir)s6 101 2 partitioned > %(outdir)s.8.log ;
		   restartIMAGE.pl %(outdir)s8 151 2 partitioned > %(outdir)s.10.log ;
		   image_run_summary.pl %(outdir)s > %(outdir)s.iteration_summary_stats ;
		   cd %(outdir)s10 ;
		   contigs2scaffolds.pl new.fa new.read.placed 300 0 %(outdir)s.final.IMAGE ;
		   touch %(outdir)s ;
                ''' % locals()
    P.run()

#######################
#Run ICORN
#######################
#ln -s ../../S1-1-1.paired_* .
#ln -s ../runIMAGE.dir/S1-SPAdes.dir/ite5/S1-SPAdes.IMAGE.scaffolds.ite5.fa .

#For usage:
#icorn.start.sh
#Only run one ICORN instance per directory!

#icorn.start.sh <Sequence to correct> <iteration start> <iteration end> <fastq_1> <fastq_2> <Insert range i.e. 100,400> <mean insert size>
#nohup icorn.start.sh S1-SPAdes.IMAGE.scaffolds.ite5.fa 1 6 S1-1-1.paired_1.fastq S1-1-1.paired_2.fastq 50,490 385 &
# IMAGE.dir/S1-1-1.IMAGE/S1-1-1.ite18/S1-1-1.ite.final.IMAGE.fa
#ICORN outputs 
	#output.perfectcoverage.txt
	#PerfectCoverageplots/PerfectMapping.Chromosome.plot.gz
	#PerfectCoverageStats.txt
	#Stats.Correction.csv
	#Stats.Mapping.csv
	#Stats.2.Correction.csv
	#ICORN.overview.txt
#Problems:
#cat ICORN.dir/S1-1-1.ICORN/ICORN.overview.txt 
#The file PerfectCoverageStats.txt doesn't exist. Please be sure to run the getPerfectcoverage.2lanes.sh script.

@follows(runIMAGE)
@mkdir('ICORN.dir')
@mkdir('*.paired_*', formatter('.+/(?P<basename>.*).paired_(.).fastq'), ['{path[0]}/ICORN.dir/{basename[0]}.ICORN'])
@jobs_limit(1)
@transform('/ifs/projects/proj018/runs/PAGIT_runs.dir/PAGIT_on_Ruffus_all_3.dir/IMAGE.dir/S*/*.ite10/*.ite.final.IMAGE.fa', regex(r'/ifs/projects/proj018/runs/PAGIT_runs.dir/PAGIT_on_Ruffus_all_3.dir/IMAGE.dir/(.*).IMAGE/(.*).ite10/(.*).ite.final.IMAGE.fa'), r'\1.ICORN', r'\1.paired_1.fastq', r'\1.paired_2.fastq', r'ICORN.dir/\1.ICORN', r'\1')
#@transform([''.join(os.path.abspath(os.getcwd())+'/ICORN.dir/S*/*.ite18/*.ite.final.IMAGE.fa')], regex(r"''.join(os.path.abspath(os.getcwd())+'/IMAGE.dir/(.*).IMAGE/(.*).ite18/(.*).ite.final.IMAGE.fa')"), r'\1.ICORN', r'../../\1.fastq.1', r'../../\1.fastq.2',  r'ICORN.dir/\1.ICORN')
def runICORN(infile, outdir, extra_parameters1, extra_parameters2, extra_parameters3, extra_parameters4): 
    
   # job_options ='-l mem_free=30G' 
    to_cluster = True
    statement = '''cd %(extra_parameters3)s  ;
                   ln -s /ifs/projects/proj018/backup/gzip_temp/%(extra_parameters4)s.fastq* . ;
                   rename .fastq.1 .paired_1.fastq *1 ;
                   rename .fastq.2 .paired_2.fastq *2 ;
		   icorn.start.sh %(infile)s 1 3 %(extra_parameters1)s %(extra_parameters2)s 50,490 385 > %(outdir)s.log ;
		   icorn.start.sh %(infile)s 4 6 %(extra_parameters1)s %(extra_parameters2)s 50,490 385 > %(outdir)s.2.log ; 
		   mv Final.icorn.fasta %(outdir)s.final.fasta ;
		   cd ../.. ;
		   touch %(outdir)s
                ''' % locals()
    P.run()


#######################
#Run RATT
#######################

#start.ratt.sh
#http://ratt.sourceforge.net/documentation.html#Installation

@follows(runICORN)
@mkdir('RATT.dir')
@mkdir('*.paired_*', formatter('.+/(?P<basename>.*).paired_(.).fastq'), ['{path[0]}/RATT.dir/{basename[0]}.RATT'])
@mkdir('RATT.dir/EMBL_annotations.dir')
def EMBLlinks():
    to_cluster = False
    statement = '''cd RATT.dir/EMBL_annotations.dir/ ;
                   ln -s /ifs/projects/proj018/runs/P.aeruginosa_other_references.dir/PAO1_annotation_NC_002516.embl . ;
		'''
    P.run()

#nohup start.ratt.sh ./EMBL_annotations.dir/ S1-SPAdes.IMAGE.scaffolds.ite5.fa S1-SPAdes.from_IMAGE.multiple_RATT Multiple > ratt.output.txt &
#Final RATT file is 'userprefix.fastaHeader'.final.embl

@follows(EMBLlinks)
@transform('/ifs/projects/proj018/runs/PAGIT_runs.dir/PAGIT_on_Ruffus_all_3.dir/ICORN.dir/S*/*.final.fasta', regex('/ifs/projects/proj018/runs/PAGIT_runs.dir/PAGIT_on_Ruffus_all_3.dir/ICORN.dir/S(.*)/(.*).ICORN.final.fasta'), r'\2', r'RATT.dir/\2.RATT/')
def runRATT(infile, outfile, extra_parameters1): 
    
    job_options = '-l mem_free=30G'
    to_cluster = True
    statement = '''cd %(extra_parameters1)s ; 
		   start.ratt.sh ../EMBL_annotations.dir %(infile)s %(outfile)s.RATT Strain > %(outfile)s.RATT.log 
                ''' % locals()
    P.run()


@follows(runRATT)
@files(None, 'QUAST_post_RATT', 'PAO1_nolines.fa')
def runQUASTpostRATT(infile, outfile, extra_parameters1): 
    to_cluster = True
    statement = ''' quast.py
                    -R %(extra_parameters1)s
                    -t 5
		    --strict-NA
		    --gene-finding
		    --gage
		    -o %(outfile)s 
                    ICORN.dir/S*/*.final.fasta  
                ''' % locals()
    P.run()




##############Comparative genomics
#TO DO: run mugsy, convert maf to phy, and others (Sibelia, C-Sibelia)


#Run jModelTest2 to identify best model
#This cmd checks all available models and best fits tests:

#jModelTest -d mugsy.refined.phy -f -i -g 4 -s 11 -AIC -BIC -AICc -DT -p -a -w -v

#-f        include models with unequals base frecuencies (e.g., -f) (default is false)
#-i         include models with a proportion invariable sites (e.g., -i) (default is false)
#-g numberOfCategories          include models with rate variation among sites and number of categories (e.g., -g 8) (default is false & 4 categories)
#-o outputFileName
#-p         calculate parameter importances
#-s numberOfSubstitutionSchemes         number of substitution schemes (e.g., -s 11) (it has to be 3,5,7,11,203; default is 3)
#-v         do model averaging and parameter importances (e.g., -v) (default is false)
#-w         write PAUP block (e.g., -w) (default is false)
#-a         estimate model-averaged phylogeny for each active criterion (e.g., -a) (default is false)


#@follows(runRATT)
#@transform('*.phy', suffix('.phy'), '.')
@files('mugsy.refined.phy','jModelTest.output')
def runjModelTest2(infile, outfile): 
    to_cluster = True
    job_options = '-l mem_free=30G'    
    statement = ''' jModelTest
		    -d %(infile)s 
		    -f -i -g 4 -s 11 
		    -AIC -BIC -AICc -DT
		    -p -a -w -v
		    -tr 12
		    -o %(outfile)s
                ''' % locals()
    P.run()






##############1.7 Variant calling with samtools and GATK

## Seqanswers http://seqanswers.com/forums/showthread.php?t=21116&highlight=haploid+variant+calling: 
##We found that the following workflow works well for variant calling in haploids using samtools. To generate a consensus sequence from the realigned BAM file in FASTQ format, we used:

##samtools mpileup -uf ref.fasta realigned.bam | bcftools view -cg -s sample.txt - | perl vcfutils.pl vcf2fq > consensus.fq

##where sample.txt is a plain text file with the sample name in the first column and the ploidy level in the second column. For example,

##sample_name 1

##The columns should be tab-delimited. 

#Second example http://ged.msu.edu/angus/tutorials-2013/snp_tutorial.html:
#samtools mpileup -uD -r 2L:100,000-150,000 -f /data/drosophila/dmel-all-chromosome-r5.37.fasta \
#/data/snp_calling/RAL357_full_bwa.sorted.bam /data/snp_calling/RAL391_full_bwa.sorted.bam | \
#bcftools view -bvcg - > RAL_samtools.raw.bcf

#-u 	output into an uncompressed bcf file (rather than compressed)
#-D  	keep read depth for each sample
#-r  	specify which chromosome region to call SNPs for (you can omit this if you want to do the whole genome, but in the interest of speed, we picked a 50kb region)
#-f 	reference genome file
#Pipe output to bcftools, which does our SNP calling based on those likelihoods.
#-b 	output to BCF format (rather than VCF); 
#-c 	do SNP calling
#-v 	only output potential variant sites (i.e., exclude monomorphic ones); 
#-g 	call genotypes for each sample in addition to just calling SNPs. 
#Then run:
#bcftools view RAL_samtools.raw.bcf | vcfutils.pl varFilter -D100
#   > RAL_samtools.vcf

#Convert the BCF file into a VCF file , 
#and pipe to vcfutils.pl with the varFilter -D100 option to filters SNPs that had read depth higher than 100 
#(high coverage SNPs may represent variation between variable copy number repeats, 

#Third example http://samtools.sourceforge.net/mpileup.shtml:
#samtools mpileup -uf ref.fa aln1.bam aln2.bam | bcftools view -bvcg - > var.raw.bcf  
#bcftools view var.raw.bcf | vcfutils.pl varFilter -D100 > var.flt.vcf  

#######################
#Run Samtools
#######################

@transform('*.fasta', suffix('.fasta'), '.fasta.fai')
def createIndexSamtools(infile, outfile):
    to_cluster = True
    statement = ''' samtools faidx %(infile)s ''' % locals()

    P.run()

#@follows(createIndexSamtools)
@transform('*bowtie.bam', suffix('.bam'), '.position_sorted.bam')
def sortBAMSamtools(infile, outfile):
    to_cluster = True
    outfile = P.snip(outfile, '.bam')
    statement = ''' samtools sort %(infile)s %(outfile)s ''' % locals()
   
    P.run()


#@follows(sortBAMSamtools)
@transform('*sorted.bam', regex(r'(.*).bam'), r'\1.deduplicated.bam', r'\1.dedup_metrics.picard')
def deduplicatePicard(infile, outfile, extra_parameters1):
    to_cluster = True
    job_options = '-l mem_free=30G' 
    statement = ''' MarkDuplicates I=%(infile)s O=%(outfile)s REMOVE_DUPLICATES=true METRICS_FILE=%(extra_parameters1)s
		    ; samtools index %(outfile)s
		''' % locals()
    
    P.run()


#Check minimum quality scores for filtering
  
@follows(deduplicatePicard, mkdir('Samtools.dir'))
@transform('*deduplicated.bam', suffix('.bam'), '.bcf', '*.fa')
#@files(None, 'all_baseline_samples.bcf', '*.fa')
def callSamtoolsMpileupBCFtoolsView(infiles, outfiles): 
    to_cluster = True
    statement = ''' samtools mpileup 
    		    -P ILLUMINA -up -t DP 
		    -f %(extra_parameters1)s %(infiles)s | bcftools call -cv -O u -o %(outfiles)s
		    ; checkpoint
 		''' % locals()
    P.run()


#@follows(deduplicatePicard, mkdir('Samtools.dir'))
#@transform('*deduplicated.bam', suffix('.bam'), '.bcf')
#@files(None, 'all_baseline_samples.bcf', '*.fa')
#def callSamtoolsMpileupBCFtoolsView(none, outfile): 
#    to_cluster = True
#        statement = ''' samtools mpileup 
#                        -P ILLUMINA -up -t DP 
#                        -f %(extra_parameters1)s S1-1-1.bowtie.position_sorted.deduplicated.bam S24_1-1-1.bowtie.position_sorted.deduplicated.bam S30-1-1.bowtie.position_sorted.deduplicated.bam S40-1-1.bowtie.position_sorted.deduplicated.ba
#	                ; checkpoint
#                    ''' % locals()
#   P.run()



#Should run varFilter with  -D100 (max depth)  -d10 (min depth)
#@follows(callSamtoolsMpileupBCFtoolsView)
@transform('*.bcf', suffix('.bcf'), '.vcf')
def callVCFutilsVarFilter(infile, outfile): 
    to_cluster = True
    statement = ''' bcftools view %(infile)s 
 		    | vcfutils.pl varFilter -D100 -d10 -p
		    > %(outfile)s 
		    ; checkpoint
                ''' % locals()
    P.run()

###########################
#Samtools has changed: update and cleanup needed:

#@transform('*.bam', suffix('.bam'), '.sorted.bam', '*.fa')
@transform('*.bam', suffix('.bam'), '.sorted.bam')
def samtoolsSort(infile, outfile, extra_parameters1):

    job_options ='-l mem_free=50G'
    to_cluster = True
    statement = ''' 
                    samtools index %(infile)s ; checkpoint ;
                    samtools sort %(infile)s -O bam -T samtools.temp > %(outfile)s ; checkpoint ; 
                    samtools index %(outfile)s ; checkpoint ;
                    bedtools genomecov -ibam %(outfile)s | grep ^genome > %(outfile)s.coverage.txt ; checkpoint 
                '''
    P.run()


#@follows(deduplicatePicard)
@follows(samtoolsSort)
@transform('*.sorted.bam', suffix('.sorted.bam'), '', '*.fa')
def callVariants(infile, outfile, extra_parameters1):

    job_options ='-l mem_free=50G'
    to_cluster = True
    statement = ''' 
                    samtools mpileup -uf %(extra_parameters1)s %(infile)s | bcftools view -bvcg - > %(outfile)s_vars_raw.bcf ; checkpoint ;
                    bcftools view %(outfile)s_vars_raw.bcf | vcfutils.pl varFilter -d4 -D 200 > %(outfile)s_vars_filt.vcf ; checkpoint ;
               '''
    P.run()


###########################


#@follows(callVCFutilsVarFilter)
@follows(callVariants)
@transform('*.vcf', suffix('.vcf'), '.SNPs_stats')
def countSNPs(infile, outfile): 
    to_cluster = True
    statement = ''' cat %(infile)s | grep -v '#' | wc -l > %(outfile)s 
                ''' % locals()
    P.run()


###########################
###########################
#Samtools update 11/05/2015:
## Call variants with samtools after mapping with bowtie2 with local realignment (--local option):
# Check tutorial and workflow with bwa mapping -> GATK raw gapped alignment -> samtools
# http://www.htslib.org/workflow/#mapping_to_variant

# TO DO: Clean up this pipeline: separate assembly, variant calling, etc. 
# TO DO: Clean up samtools calls: old samtools and an updated samtools (done quickly) and latest workflow below.

@mkdir('tmp')
@transform('*.bam', suffix('.bam'), '.sorted.bam')
def fixAndSort(infile, outfile):

    to_cluster = True
    job_options = '-l mem_free=50G'
    statement = '''samtools sort -n %(infile)s -O bam -T tmp/%(infile)s.tmp | samtools fixmate -O bam - %(infile)s.fixmate.bam ; checkpoint ;
    		   samtools sort %(infile)s.fixmate.bam -O bam -T tmp/%(infile)s.fixmate.bam.tmp > %(outfile)s ; rm -f *fixmate.bam ; checkpoint 
		''' % locals()
    P.run()

@follows(fixAndSort)
@transform('*.sorted.bam', suffix('.sorted.bam'), '.readGroups.bam', '../PAO1_nolines.fa')
def GATKReadGroups(infile, outfile, genome,
                   library="unknown", platform="Illumina",
                   platform_unit="1", track="unknown"):
    '''Reorders BAM according to reference fasta and adds read groups'''

    if track == 'unknown':
        track = P.snip(os.path.basename(infile), "sorted.bam")
    tmpdir_gatk = P.getTempDir('.')
    job_options = '-l mem_free=1.4G -l picard=1'
    statement = '''ReorderSam
                    INPUT=%(infile)s
                    OUTPUT=%(tmpdir_gatk)s/%(track)s.reordered.bam
                    REFERENCE=%(genome)s
                    ALLOW_INCOMPLETE_DICT_CONCORDANCE=true
                    VALIDATION_STRINGENCY=SILENT ; checkpoint ;''' % locals()

    statement += '''samtools index %(tmpdir_gatk)s/%(track)s.reordered.bam ;
                    checkpoint ;''' % locals()

    statement += '''AddOrReplaceReadGroups
                    INPUT=%(tmpdir_gatk)s/%(track)s.reordered.bam
                    OUTPUT=%(outfile)s
                    RGLB=%(library)s
                    RGPL=%(platform)s
                    RGPU=%(platform_unit)s
                    RGSM=%(track)s
                    VALIDATION_STRINGENCY=SILENT ; checkpoint ;''' % locals()

    statement += '''samtools index %(outfile)s ;
                    checkpoint ;''' % locals()
    statement += '''rm -rf %(tmpdir_gatk)s ;''' % locals()

    P.run()

## Quality and INDEL alignment improvements:
# Sort bam files into coordinate order, quality score adjustment with GATK (analyse covariates and compensate for those 
#covariates by adjusting quality scores):

@follows(GATKReadGroups)
@transform('*.readGroups.bam', suffix('.readGroups.bam'), '.realigned_1.bam', '../PAO1_nolines.fa')
def GATKrealign(infile, outfile, extra_parameters1):

    to_cluster = True
    job_options = '-l mem_free=1.4G -l picard=1'
    statement = '''GenomeAnalysisTK -T RealignerTargetCreator -R %(extra_parameters1)s -I %(infile)s -o %(infile)s.intervals ; checkpoint ;
    		   GenomeAnalysisTK -T IndelRealigner -R %(extra_parameters1)s -I %(infile)s -targetIntervals %(infile)s.intervals -o %(outfile)s
		''' % locals()
    P.run()


#Requires a known variants vcf file: 
@follows(GATKrealign)
@transform('*.realigned_1.bam', suffix('.realigned_1.bam'), '.recal.bam', '../*.fa')
def GATKRecalibrate(infile, outfile, extra_parameters1):

    to_cluster = True
    job_options = '-l mem_free=1.4G -l picard=1'
    statement = '''GenomeAnalysisTK -T BaseRecalibrator -R %(extra_parameters1)s -I %(infile)s -o %(infile)s.recal.table ; checkpoint; 
		   GenomeAnalysisTK -T PrintReads -R %(extra_parameters1)s -I %(infile)s --BSQR %(infile)s.recal.table -o %(outfile)s
		''' % locals()
    P.run()

# Mark PCR and optical duplicates with Picard Markduplicates, calculate read depth, compile into one bam, merge, sort and index:

#@follows(GATKRecalibrate)
@follows(GATKrealign)
#@transform('*.recal.bam', regex(r'(.*).recal.bam'), r'\1.deduplicated.bam', r'\1.dedup_metrics.picard')
@transform('*.realigned_1.bam', regex(r'(.*).realigned_1.bam'), r'\1.deduplicated.bam', r'\1.dedup_metrics.picard')
def deduplicateAndCoverage(infile, outfile, extra_parameters1):
    
    to_cluster = True
    job_options = '-l mem_free=1.4G -l picard=1'
    statement = ''' MarkDuplicates I=%(infile)s O=%(outfile)s REMOVE_DUPLICATES=true METRICS_FILE=%(extra_parameters1)s ; checkpoint ;
		    bedtools genomecov -ibam %(outfile)s | grep ^genome > %(outfile)s.coverage.txt ; checkpoint
		''' % locals()
    
    P.run()


@follows(deduplicateAndCoverage)
@files(None, 'all-samples-merged.bam')
def mergeSamtools(infile, outfile):
    
    to_cluster = True
    statement = '''ls -1 *deduplicated.bam > bam_file_list.txt ;
    		   samtools merge -b bam_file_list.txt %(outfile)s ; checkpoint ;
		   samtools sort %(outfile)s -O bam -T .tmp  > %(outfile)s.sorted.bam ; checkpoint ;
		   samtools index %(outfile)s.sorted.bam ; checkpoint ;
		   touch %(outfile)s
		''' % locals()
    
    P.run()


@follows(mergeSamtools)
#@transform('*.sorted.bam', suffix('.sorted.bam'), '.realigned_2.bam', '*.fa')
@files('all-samples-merged.bam.sorted.bam', 'all-samples-merged.realigned_2.bam', '../PAO1_nolines.fa')
def GATKrealign2(infile, outfile, extra_parameters1):

    to_cluster = True
    job_options = '-l mem_free=1.4G -l picard=1'
    statement = '''GenomeAnalysisTK -T RealignerTargetCreator -R %(extra_parameters1)s -I %(infile)s -o %(infile)s.intervals ; checkpoint ;
    		   GenomeAnalysisTK -T IndelRealigner -R %(extra_parameters1)s -I %(infile)s -targetIntervals %(infile)s.intervals -o %(outfile)s ; checkpoint ;
		   samtools index %(outfile)s 
		''' % locals()
    P.run()


#@follows(deduplicateAndCoverage)
@transform('*deduplicated.bam', suffix('deduplicated.bam'), 'sorted.bam')
def sortSamtools(infile, outfile):
    
    to_cluster = True
    statement = '''samtools sort %(infile)s -O bam -T %(infile)s.tmp > %(outfile)s ; checkpoint ;
		   samtools index %(outfile)s ; checkpoint 
		''' % locals()
    
    P.run()



## Variant calling:
# Produce bcf file and call variants, index bcf, produce stats from variant calls and plot results:
# Run either merging all samples and calling on single bam or calling samples individually.

#@follows(mergeSamtools, mkdir('vcf_stats_plots.dir'))
@follows(GATKrealign2, mkdir('vcf_stats_plots.dir'))
@files('all-samples-merged.realigned_2.bam', 'all-samples-merged.vcf.gz', '../PAO1_nolines.fa', 'vcf_stats_plots.dir')
def callVariantsMerged(infile, outfile, extra_parameters1, extra_parameters2):

    job_options = '-l h=cgatgpu1'
    to_cluster = True
    statement = '''samtools mpileup -d 100 -P Illumina -ugpf %(extra_parameters1)s %(infile)s | bcftools call -vmO z -o %(outfile)s ; checkpoint ;
		   tabix -p vcf %(outfile)s ; checkpoint ; 
		   bcftools stats -F %(extra_parameters1)s -s - %(outfile)s > %(outfile)s.stats ; checkpoint ;
		   plot-vcfstats -s -p %(extra_parameters2)s %(outfile)s.stats
                ''' % locals()
    P.run()


@follows(sortSamtools)
@transform('*.sorted.bam', suffix('.sorted.bam'), '.vcf.gz', '../PAO1_nolines.fa')
def callVariantsSingleFiles(infile, outfile, extra_parameters1):

    job_options ='-l mem_free=25G'
    to_cluster = True
    statement = '''samtools mpileup -d 100 -P Illumina -ugpf %(extra_parameters1)s %(infile)s | bcftools call -vmO z -o %(outfile)s ; checkpoint ;
		   tabix -p vcf %(outfile)s ; 
		   bcftools stats -F %(extra_parameters1)s -s - %(outfile)s > %(outfile)s.stats 
                ''' % locals()
    P.run()


@follows(callVariantsMerged)
#@follows(callVariantsSingleFiles)
#@mkdir('*.sorted.bam', suffix('.sorted.bam'), r'\1.vcf_plots.dir')
@transform('*.vcf.gz.stats', suffix('.vcf.gz.stats'), r'\1.vcf_plots.dir/')
def vcfStats(infile, outdir):

    job_options ='-l h=cgatgpu1'
    to_cluster = True
    statement = '''plot-vcfstats -s -p %(outdir)s %(infile)s
                ''' % locals()
    P.run()


## Benchmark pipeline results:
#http://www.bioplanet.com/gcat/reports/119/variant-calls/ion-torrent-225bp-se-exome-30x/bowtie2-gatk-ug/compare-183-169/group-read-depth


#######################
#Run GATK
#######################



#######################
#Check Genome coverage
#######################

@follows(mkdir('Genome_coverage.dir'))
@transform('bowtie.dir/*.bam', suffix('.bam'), '.genome_coverage')
def getGenomeCoverage(infile, outfile): 
    to_cluster = True
    statement = ''' genomeCoverageBed -ibam %(infile)s -g PAO1_nolines.fa | grep genome > %(outfile)s 
		    ; checkpoint
		    ; mv %(outfile)s Genome_coverage.dir/
                ''' % locals()
    P.run()


#######################
#Phylogenetics from SNPs
#######################

@transform('*.vcf.gz', suffix('.vcf.gz'), '.SNPs_in_fasta')
def vcfToFasta(infile, outfile): 
    to_cluster = True
    statement = '''zcat %(infile)s | vcf-to-tab > %(outfile)s ; checkpoint ; 
                ''' % locals()
    P.run()





#######################
#Strain specific mapping
#######################



#######################
#Check Genome coverage
#######################


## Commands to set up strain specific mapping:

# Create subdirectories:
#mkdir P001.dir # there are 16 patients, I've created these manually as groupings don't have common prefixes.

# Soft links to all fastqs (will use the processed samples, have done this manually for now for each patient directory and each sample):
#ln -s /ifs/projects/proj018/runs_second_round.dir/readQC_processed_2.dir/processed.dir/*.gz .
#ln -s /ifs/projects/proj018/runs_second_round.dir/readQC_processed.dir/processed.dir/trimmomatic-S03_TP_D7_001_TP_D5_002-1-1.fastq* P001.dir/

# Soft link of strain specific assembly:
#ln -s /ifs/projects/proj018/runs_second_round.dir/assembly.dir/SPAdes.dir/trimmomatic_S*.SPAdes/trimmomatic_S*.SPAdes.fasta P001.dir/


# Trimmomatic processed files require renaming:
#rename - _ *.gz

# Rename / create link for commands downstream (regex otherwise) of strain assembly: 
#ln -s *SPAdes strain_specific_assembly.fa

# Config pipeline:
#python /ifs/devel/antoniob/CGATPipelines/CGATPipelines/pipeline_mapping.py config

# TO DO: Copy ini file over. This has bowtie2 setup, ideally I need to create an annotations file and database for each (using PAO1 for now, TO DO: check this one is actually OK, then create specific annotations if possible [from predictions, maybe not useful though]):
#cp ../../mapping_processed_2.dir/pipeline.ini .

# Run pipeline:
#nohup python /ifs/devel/antoniob/CGATPipelines/CGATPipelines/pipeline_mapping.py -v 5 -p 20 --force make full &

# Build report:
#nohup python /ifs/devel/antoniob/CGATPipelines/CGATPipelines/pipeline_mapping.py -v 5 -p 20 --force make build_report &

# Get genome coverage:
#nohup python /ifs/devel/antoniob/projects/Pseudomonas-018/pipeline_genome_assembly.dir/src/pipeline_genome_assembly.py -v 5 -p 20 --force make getGenomeCoverage &

# TO DO: Plot coverage, doesn't work:
#cgatreport-test -t coverage --path ./Genome_coverage.dir -r r-ggplot -o statement='aes(coverage,fraction,color=Slice)+geom_point()+xlim(0,150)+ylab("fraction of genome")+xlab("coverage(no.reads)")'

# Also check:
#cgatreport-test -t coverage --path ./Genome_coverage.dir -r r-ggplot -o statement='aes(coverage,fraction,color=Slice)+geom_point()+xlim(0,150)+ylab("fraction of genome")+xlab("coverage(no.reads)")' -m cumulative


## Strain specific variant calling:
# Create a subdirectory within each strain's directory (P0xx.dir):
#mkdir samtools.dir

# Rename / create link for commands downstream (regex otherwise) of strain assembly: 
#ln -s ../*SPAdes.fasta strain_specific_assembly.fa

# Create links to mapped files and indexes:
#ln -s ../bowtie.dir/*bai .
#ln -s ../bowtie.dir/*bam .

# Create bowtie2 indexes for strain assembly:
#bowtie2-build strain_specific_assembly.fa strain_specific_assembly &> strain_specific_assembly.stout

# TO DO: Run SNP and INDEL calling with samtools (needs clean-up, see commands below):
#nohup python /ifs/devel/antoniob/projects/Pseudomonas-018/pipeline_genome_assembly.dir/src/pipeline_genome_assembly.py -v 3 -p 20 --force make countSNPs &




#########################################################################
# Load statistics: readQC pipeline example:

#@jobs_limit( 1, "db" )
#@transform( runFastqc, suffix(".fastqc"), "_fastqc.load" )
#def loadFastqc( infile, outfile ):
#    '''load FASTQC stats.'''
#    
#    track = P.snip( infile, ".fastqc" )
#
#    def section_iterator( infile ):
#
#        data = []
#        for line in infile:
#            if line.startswith( ">>END_MODULE" ): 
#                yield name, status, header, data
#            elif line.startswith(">>"):
#                name, status = line[2:-1].split("\t")
#                data = []
#            elif line.startswith("#"):
#                header = "\t".join([ x for x in line[1:-1].split("\t") if x != ""] )
#            else:
#                data.append( "\t".join([ x for x in line[:-1].split("\t") if x != ""] ) )



###################################################################
###################################################################
###################################################################
## primary targets
###################################################################
@follows( loadDummyTask )
def full(): pass

###################################################################
###################################################################
###################################################################
## primary targets
###################################################################

@follows( mkdir( "report" ) )
def build_report():
    '''build report from scratch.'''

    E.info( "starting report build process from scratch" )
    P.run_report( clean = True )

@follows( mkdir( "report" ) )
def update_report():
    '''update report.'''

    E.info( "updating report" )
    P.run_report( clean = False )

@follows( update_report )
def publish_report():
    '''publish report.'''

    E.info( "publishing report" )
    P.publish_report()

if __name__== "__main__":

    # P.checkFiles( ("genome.fasta", "genome.idx" ) )
    sys.exit( P.main(sys.argv) )


################################################################################
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
Pipeline MinION
===========================

:Author: Antonio Berlanga-Taylor
:Release: $Id$
:Date: 10 Nov 2014
:Tags: Python, MinION, nanopore

A pipeline template.

Overview
========

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
|                    |                   |                                                |
+--------------------+-------------------+------------------------------------------------+

Pipeline output
===============

The major output is in the database file :file:`csvdb`.

Example
=======

Example data is available at http://www.cgat.org/~andreas/sample_data/pipeline_MinION.tgz.
To run the example, simply unpack and untar::

   wget http://www.cgat.org/~andreas/sample_data/pipeline_MinION.tgz
   tar -xvzf pipeline_MinION.tgz
   cd pipeline_MinION
   python <srcdir>/pipeline_MinION.py make full

.. note:: 
   For the pipeline to run, install the :doc:`pipeline_annotations` as well.

Glossary
========

.. glossary::


Code
====

"""
from ruffus import *

import sys, glob, gzip, os, itertools, re, math, types, collections, time, csv
import optparse, shutil
import sqlite3
import rpy2.robjects as robjects
from rpy2.robjects import r as R
#import pandas as pandas
#from pandas import *
from numpy import *
from matplotlib import *


import CGAT.Experiment as E
import CGAT.IOTools as IOTools
import CGAT.Database as Database

###################################################
##pandas settings
###################################################

#pandas.set_option('display.mpl_style', 'default') # Better looking graphs
#figsize(20, 10)
#matplotlib.rc('xtick', labelsize=12) 
#matplotlib.rc('ytick', labelsize=12)
#font = {'family' : 'normal', 'weight' : 'bold', 'size' : 12}
#matplotlib.rc('font', **font)

###################################################
###################################################
###################################################
## Pipeline configuration
###################################################

# load options from the config file
import CGAT.Pipeline as P
P.getParameters( 
    ["%s.ini" % __file__[:-len(".py")],
     "../pipeline.ini",
     "pipeline.ini" ] )

PARAMS = P.PARAMS
PARAMS_ANNOTATIONS = P.peekParameters( PARAMS["annotations_dir"],
                                       "pipeline_annotations.py" )

###################################################################
###################################################################
## Helper functions mapping tracks to conditions, etc
###################################################################
import CGATPipelines.PipelineTracks as PipelineTracks

# define some tracks if needed
TRACKS = PipelineTracks.Tracks( PipelineTracks.Sample ).loadFromDirectory( 
    glob.glob("*.ini" ), "(\S+).ini" )


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
## worker tasks
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



###################################################################
###################################################################

#MinION results processing after basecalling from Metrichor

#Check example: https://nanoporetech.atlassian.net/wiki/display/~JLord/2014/10/30/Cruchaga+Lab+Lambda+Burn+In+Run

#poretools stats --type fwd /downloads_072514/ > burnin_fwd_072514
#poretools stats --type rev /downloads_072514/ > burnin_rev_072514

#Process, plot and inspect fast5 results:
#poretools combine -o second_flowcell_run_1_all_files.fast5.tar.gz /ifs/projects/antoniob/Nanopore/burn_in/MAP_Antonio/Data/reads/downloads/the_rest/FGU097_2490_55727_13_nov_2014_run_1*


#ln -s /ifs/projects/antoniob/Nanopore/burn_in/MAP_Antonio/Data/reads/downloads/metrichor*/*.fast5 .

@collate('*.fast5', regex(r'FGU097_(.+)_RUN_(\d)_(.+)_lib_(\d)_(.+).fast5'), r'FGU097_\1_RUN_\2_\3_lib_\4.fast5.tar.gz')
def mergeFiles(infile, outfile):

    job_options ='-l mem_free=30G'
    to_cluster = True
    list1 = str(infile).replace("'", " ")
    list2 = str(list1).replace("(", "")
    list3 = str(list2).replace(")", "")
    list4 = str(list3).replace(",", " ")
    statement = ''' poretools combine -o %(outfile)s %(list4)s
		'''
    P.run()

#rm -rf *.fast5 links (thousands of files)


#@follows(mergeFiles)
#poretools squiggle generates a plot for each file, here generating for the longest read -- TO DO: poretools squiggle --theme-bw --saveas png winner_%(outfile)s.fasta ; checkpoint ; needs the fast5 format, not fasta
@jobs_limit(1)
@transform('*.tar.gz', regex('(.*).tar.gz'), r'\1')
def runPoretools(infile, outfile):
   
    job_options ='-l mem_free=30G' 
    to_cluster = True
    statement = ''' poretools stats --type=all --full-tsv %(infile)s > %(outfile)s.stats_all ; checkpoint ;
		    poretools stats --type=2D %(infile)s > %(outfile)s.stats_2D ; checkpoint ;
		    poretools fastq --type=all %(infile)s > %(outfile)s.fastq ; checkpoint ;
		    poretools hist --min-length=100 --theme-bw --saveas hist_%(outfile)s.png %(infile)s ; checkpoint ;
		    poretools tabular %(infile)s > tabular_%(outfile)s.tab ; checkpoint ;
		    poretools nucdist %(infile)s > nuc_dist_%(outfile)s.nuc_dist ; checkpoint ;
		    poretools qualdist %(infile)s > qual_dist_%(outfile)s.qualdist ; checkpoint ;
		    poretools winner %(infile)s > winner_%(outfile)s.fasta ; checkpoint ;
		    poretools times %(infile)s > times_%(outfile)s.tab ; checkpoint ;
		    poretools yield_plot --theme-bw --saveas yield_reads_%(outfile)s.png %(infile)s ; checkpoint ;
		    poretools yield_plot --theme-bw --plot-type=basepairs --saveas yield_bp_%(infile)s.png %(infile)s ; checkpoint ;
		    poretools occupancy --saveas occupancy_%(outfile)s.png %(infile)s ; checkpoint ;
		    touch %(outfile)s
		'''
    P.run()

#mkdir poretools.dir
#mv nuc_dist* *png *stats* winner* times* tabular* qual* poretools.dir/


#For alignment/variant calling:
#lastal options = -T0 local aligment ; -Q1 use sequence quality scores fastq sanger format ; -f1 output MAF ; -s2 use both strands
#maf-convert is from last 

#ln -s /ifs/mirror/ncbi/Enterobacteria_phage_lambda.fasta .
#lastdb Enterobacteria_phage_lambda.last_index Enterobacteria_phage_lambda.fasta
#samtools faidx ref.fasta

@follows(runPoretools)
@transform('*.fastq', regex('(.*).fast5.fastq'), r'\1', 'Enterobacteria_phage_lambda.last_index', 'Enterobacteria_phage_lambda.fasta')
def alignMinION(infile, outfile, extra_parameters1, extra_parameters2):

    job_options ='-l mem_free=30G'
    to_cluster = True
    statement = ''' lastal -Q1 -f1 -s2 -T0 %(extra_parameters1)s %(infile)s > %(outfile)s.maf ; checkpoint ;
		    maf-convert -d sam %(outfile)s.maf > %(outfile)s.sam ; checkpoint ;
		    samtools view -T %(extra_parameters2)s -bS %(outfile)s.sam > %(outfile)s.bam ; checkpoint ; 
		    samtools index %(outfile)s.bam ; checkpoint ;
		    samtools sort %(outfile)s.bam -O bam -T samtools.temp > %(outfile)s.sorted.bam ; checkpoint ; 
		    samtools index %(outfile)s.sorted.bam ; checkpoint ;
		    bedtools genomecov -ibam %(outfile)s.sorted.bam | grep ^genome > %(outfile)s.coverage.txt ; checkpoint ;
		    touch %(outfile)s
		'''
    P.run()

#Get coverage:
#@follows(alignMinION)
#@transform('*.coverage.txt', suffix('*.txt'), '.png')
#def plotCoverage(infile, outfile):

    #job_options ='-l mem_free=30G'
#    to_cluster = True
    #rtable = robjects.r['read.table']
    #rtable('%(infile)s') 
    #rplot =  robjects.r['plot']
    #rplot(cov[1:51,2], cov[1:51,5], type='h', col='darkgreen', lwd=3, xlab='Depth', ylab='Fraction of genome at depth',)
    # %(outfile)s
#    readTable = pandas.Dataframe(%(infiles)s)
#    plot()        
                
#    P.run()


#Call variants:
#@follows(plotCoverage)
@transform('*.sorted.bam', regex('(.*).fastq'), r'\1', '*.fasta')
def callVariants(infile, outfile, extra_parameters1):

    job_options ='-l mem_free=30G'
    to_cluster = True
    statement = ''' samtools mpileup -uf %(extra_parameters1)s %(infile)s | bcftools view -eg - | vcfutils.pl vcf2fq > %(outfile)s_consensus.fastq ; checkpoint ;
		    samtools mpileup -uf %(extra_parameters1)s %(infile)s | bcftools view -bvcg - > %(outfile)s_vars_raw.bcf ; checkpoint ;
		    bcftools view %(outfile)s_vars_raw.bcf | vcfutils.pl varFilter -d4 -D 200 > %(outfile)s_vars_filt.vcf ; checkpoint ;
		    blastn -query %(outfile)s_consensus.fasta -subject %(extra_parameters1)s -out %(outfile)s_blastn.txt ; checkpoint
               '''
    P.run()


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
def publish():
    '''publish report and data.'''

    E.info( "publishing report" )
    P.publish_report()

if __name__== "__main__":

    # P.checkFiles( ("genome.fasta", "genome.idx" ) )
    sys.exit( P.main(sys.argv) )

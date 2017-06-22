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
Pipeline fastq_qc_blast
===========================

:Author: Antonio J Berlanga-Taylor
:Release: $Id$
:Date: |today|
:Tags: Python


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

Example data is available at http://www.cgat.org/~andreas/sample_data/pipeline_fastq_qc_blast.tgz.
To run the example, simply unpack and untar::

   wget http://www.cgat.org/~andreas/sample_data/pipeline_fastq_qc_blast.tgz
   tar -xvzf pipeline_fastq_qc_blast.tgz
   cd pipeline_fastq_qc_blast
   python <srcdir>/pipeline_fastq_qc_blast.py make full

.. note:: 
   For the pipeline to run, install the :doc:`pipeline_annotations` as well.

Glossary
========

.. glossary::


Code
====

"""
from ruffus import *

import sys, glob, gzip, os, itertools, re, math, types, collections, time
import optparse, shutil
import sqlite3

import CGAT.Experiment as E
import CGAT.IOTools as IOTools
import CGAT.Database as Database

###################################################
###################################################
###################################################
## Pipeline configuration
###################################################

# load options from the config file
import CGAT.Pipeline as P

P.getParameters(
    ["%s/pipeline.ini" % os.path.splitext(__file__)[0],
     "../pipeline.ini",
     "pipeline.ini"],
    only_import=__name__ != "__main__")

PARAMS = P.PARAMS

PARAMS.update(P.peekParameters(
    PARAMS["annotations_dir"],
    "pipeline_annotations.py",
    on_error_raise=__name__ == "__main__",
    prefix="annotations_",
    update_interface=True))

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
    P.load( infile, outfile, "--add-index=word" )


###################################################################
@transform('*.gz', suffix( ".gz"), ".counts.gz")
def downSample( infile, outfile ):
    '''down sample fastq reads to avoid long runs downstream'''

    to_cluster = True

    statement = ''' zcat %(infile)s | python /ifs/devel/antoniob/cgat/scripts/fastq2fastq.py --method=sample --sample-size=0.1 > %(outfile)s '''
    P.run()



#Diamond manual:
#http://ab.inf.uni-tuebingen.de/data/software/diamond/download/public/manual.pdf
#Run once: diamond makedb -t 4 -v 3 --in nr_2015 -d nr_2015_db (creates a diamond binary database .dmnd
#/ifs/mirror/ncbi/ncbi_nr/nr_2015.dir/nr_2015_db.dmnd

#For diamond 0.6.13 diamond blastx --threads 2 --db /ifs/mirror/ncbi/ncbi_nr/nr_2015.dir/nr_2015_0.6.13 --query CRCZ_0_TP_005-1-1.fastq.1.copy -o CRCZ_0_TP_005-1-1.fastq.1.copy.dmnd --tmpdir scratch.dir

#diamond blastx --threads 5 --db /ifs/mirror/ncbi/ncbi_nr/nr_2015.dir/nr_2015_0.5.2 --query CRCZ_0_TP_005-1-1.fastq.1.copy -o CRCZ_0_TP_005-1-1.fastq.1.copy.dmnd --tmpdir scratch.dir

#For diamond 0.6.13: diamond blastx -d %(extra_parameters1)s -q %(infile)s --sam %(outfile)s -t scratch.dir --compress 1


#@follows(downSample)
@mkdir('scratch.dir')
@transform('*.gz', suffix('.gz'), '.dmnd' , 'nr_2015_0.5.2')
def runDiamond(infile, outfile, extra_parameters1):
    '''Blast reads vs database using Diamond'''

    to_cluster = True
    
    statement = ''' zcat %(infile)s > %(infile)s.out ;
		    diamond blastx -d %(extra_parameters1)s -q %(infile)s.out -o %(outfile)s --tmpdir scratch.dir
    		'''
    P.run()


@follows(runDiamond)
@transform( ["%s.ini" % x.asFile() for x in TRACKS], suffix( ".ini"), ".counts")
def runMEGAN( infile, outfile ):
    '''Take Diamond output and classify using MEGAN'''

    to_cluster = True

    statement = ''' %(infile)s > %(outfile)s'''
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

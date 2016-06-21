# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <rawcell>

# To do:
#     Integrate exo papers methods
#     Detail pseudocode
#     Code!
#     Test!
#     Push back to pipeline_peakcalling.py?

# WTCHG LCL exo vs Chip seq VDR, RFX5, CTCF. File locations:
# Code (this file):
# /ifs/devel/antoniob/projects/chIP-Exo

# Data:
# Giuseppe: /ifs/projects/antoniob/data/exoVDR_Ram/
# WTCHG LCLs: /ifs/projects/antoniob/data/exoWTCHG/
# VDR ChIP-seq: /ifs/projects/antoniob/data/chIPseqVDR/

# Runs:
# /ifs/projects/antoniob/chIP/  


# input paramaters for:
#     GEM
#     ENCODE quality metrics
#     GeneTrack
#     ?
#     Giuseppe, code quality metrics
#     
# in pipeline modules
# 
# 
# check 
# !gedit /ifs/devel/antoniob/src/pipeline_template.py
# !gedit /ifs/devel/antoniob/src/Pipeline.py
# !gedit /ifs/devel/antoniob/src/pipeline_quickstart.py

# <codecell>

#pwd

# <codecell>

#!python /ifs/devel/antoniob/src/pipeline_quickstart.py --name=exo

# <rawcell>

# Once the directory is set-up start coding within pipeline_exo.py:
# 
# ie    ChIP-Exo pipeline example (pseudocode)

# <rawcell>

# 
# First step
# Document pipeline:
# Date, author, licence, objective, overview, usage, configuration, required tools, file inputs, file naming, outputs, directory structure, glossary, example, example data, etc.

# <rawcell>

# ################################################################################
# #
# #   MRC FGU Computational Genomics Group
# #
# #   $Id: xxxxxxxxxxxxxxxxxxxxxxxxxx $
# #
# #   Copyright (C) 2009 Andreas Heger
# #
# #   This program is free software; you can redistribute it and/or
# #   modify it under the terms of the GNU General Public License
# #   as published by the Free Software Foundation; either version 2
# #   of the License, or (at your option) any later version.
# #
# #   This program is distributed in the hope that it will be useful,
# #   but WITHOUT ANY WARRANTY; without even the implied warranty of
# #   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# #   GNU General Public License for more details.
# #
# #   You should have received a copy of the GNU General Public License
# #   along with this program; if not, write to the Free Software
# #   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
# #################################################################################
# 
# 
# """
# ===========================
# ChIP-Exo Pipeline
# ===========================
# 
# :Author: Antonio J Berlanga-Taylor
# :Release: $Id$
# :Date: 01/04/2013
# :Tags: Python
# 
# A pipeline to analyse ChIP-exo data using both paired and single end reads.
# This pipeline is based on pipeline_peakcalling.py and several other tools/code from CGAT/FGU/ENCODE.
# 
# ChIP-Exo original reference: Rhee and Pugh Cell 2011
# 
# General analytical pipeline for exo data:
# 
# QC --> Mapping --> Peak calling and QC ENCODE metrics --> Motif finding --> Annotation and integration --> Study specific questions
# 
# In pipelines:
# 
# pipeline_fastqc 
# --> pipeline_mapping.py with BWA and STAMPY and QC with SAMTools and Picard (mark duplicates, multiple metrics, estimate library complexity, uniquely mapped reads)
# --> (this pipeline, peaks) uses MACS2, GEM, GeneTrack
# --> (this pipeline, QC) ENCODE metrics (sequencing saturation curve, # peaks, # peaks with motifs, # reads in peaks, FRiP, IDR, saturation, etc.) and separately visualise with IGV 
# --> pipeline_intervals.py;
# --> motif searching: PipelineMotifs.py, Motifs.py; check specific tools in src (run grep again on /ifs/devel/antoniob/src/*)
# --> Study specific pipeline (allele specific events)
# 
# Questions:
# - Relate peak locations to known genomic features (TSS, introns, exons, etc.; genetic variation); - Motif discovery and degeneracy?; - Allele specific binding; - Expression levels of assocaited genes; - Motif cluster analysis; - Evolution analysis 
# 
# 
# Overview
# ========
# 
# Usage
# =====
# 
# The peak calling pipeline applies various peak calling algorithms on mapped reads. It is designed for ChIP-Exo data, preferably paired end, and uses narrow peak callers.
# The ChIP-Exo pipeline imports reads from one or more ChIP-Exo experiments and performs the following tasks (depending on the peak caller used):
# 
#    * calls peaks 
#    * describe de-novo motifs
#    * find motifs within intervals
# 
# 
# pipeline_peakcalling takes as input reads aligned to genomic sequence as :term:`bam` formatted files and calls peaks. The pipeline implements several peak callers:
# 
# macs_ (excluded from this pipeline as it does not handle paired data or single end exo data adequately)
#    Model-based Analysis of ChIP-Seq (MACS), for identifying transcript factor binding sites. MACS captures the influence of genome complexity to evaluate the significance of enriched ChIP regions, and MACS improves the spatial resolution of binding sites through combining the information of both sequencing tag position and orientation. MACS can be easily used for ChIP-Seq data alone, or with control sample with the increase of specificity.
# 
# macs2_ (used as a substitute for macs1.4 and seems to run well on paired and on exo data)
#    MACS 2 is the new release of the MACS peak caller. Among other improvements it adds support for handling paired end reads.
# 
# GEMS_
# 
# 
# GeneTrack_
# 
# Excluded:
# spp_ (may be sueful, particularly as ENCODE metrics packages can run on this but haven't tested properly on exo data)
#    SPP is a R package especially designed for the analysis of Chip-Seq data from Illummina platform. The package was developed by Peter Park's group from Harvard Medical School.
# 
# Don't seem to work (can't handle lack of controls, paired data or exo data or a combination of these):
# zinba_
# sicer_narrow_
# sicer_broad_
# peakranger_ranger_
# peakranger_ccat_
# 
# Check pipeline_peakcalling.py for more on these peakcallers.
# 
# Peak callers have different strengths and weaknesses. Some might work well on broad peaks such as some histone marks, others work better for narrow, sharp peaks. Many callers these days attempt to call both types of peaks.
# 
# See :ref:`PipelineSettingUp` and :ref:`PipelineRunning` on general information how to use CGAT pipelines.
# 
# 
# Configuration
# -------------
# 
# The pipeline requires a configured :file:`pipeline.ini` file. 
# 
# The sphinxreport report requires a :file:`conf.py` and :file:`sphinxreport.ini` file  (see :ref:`PipelineDocumenation`). To start with, use the files supplied with the :ref:`Example` data.
# 
# 
# Input
# -----
# 
# Mapped reads
# ++++++++++++
# 
# The principal input of this pipeline is a collection of reads mapped to a reference genome. 
# Mapped reads are imported by placing files are linking to files in the :term:`working directory`.
# 
# The default file format assumes the following convention:
# 
#    <sample>-<condition>-<replicate>.genome.bam
# 
# ``sample`` and ``condition`` make up an :term:`experiment`, while ``replicate`` denotes the :term:`replicate` within an :term:`experiment`. 
# 
# Please not the suffix ``genome.bam`` which is required to distinguish the input :term:`bam` formatted files from those that are created in the pipeline after duplication removal (``.prep.bam`` and ``.call.bam``)
# 
# 
# Optional inputs
# +++++++++++++++
# 
# Requirements
# ------------
# 
# The pipeline requires the results from :doc:`pipeline_annotations`. Set the configuration variable 
# :py:data:`annotations_database` and :py:data:`annotations_dir`.
# 
# On top of the default CGAT setup, the pipeline requires the following software to be in the path:
# 
# +--------------------+-------------------+------------------------------------------------+
# |*Program*           |*Version*          |*Purpose*                                       |
# +--------------------+-------------------+------------------------------------------------+
# GEM (already setup for CGAT)
# 	Peak caller and motif finder, check http://cgs.csail.mit.edu/gem/ |
# 
# MACS2 (already setup for CGAT)
# 	version 1.4 documentation: http://liulab.dfci.harvard.edu/MACS/ 
# 
# GeneTrack:
# 	singleton filtering (peak on one strand only), peak pair distance calculation (maximum distance used to determine peak pairs), peak pairing (match peaks on opposing strands), tag count criteria (minimum number of reads to call a peak), minimum number of replicates to use to call a peak. Check Rhee and Pugh Cell 2011 original method paper for criteria in supplementary methods.
# 
# 	https://code.google.com/p/genetrack/
# 	http://bioinformatics.oxfordjournals.org/content/24/10/1305.full
# 
# 
# Pipline Output
# ==============
# 
# The results of the computation are all stored in an sqlite relational database :file:`csvdb`.
# 
# The major output is in the database file :file:`csvdb`. For each peak caller there are tables called:
# 
# <track>_<caller>_regions
# <track>_<caller>_summits
# 
# Each of these tables contains the following columns:
# 
# +------------------+--------------------------------------------------------------+
# |*Column*          |*Content*                                                     |
# +------------------+--------------------------------------------------------------+
# |avgval            |Average read depth in interval                                |
# +------------------+--------------------------------------------------------------+
# |contig            |Contig                                                        |
# +------------------+--------------------------------------------------------------+
# |control_avgval    |Average read depth in control within inter val                |
# |                  |                                                              |
# +------------------+--------------------------------------------------------------+
# |control_length    |Interval length                                               |
# +------------------+--------------------------------------------------------------+
# |control_npeaks    |Number of peaks in control                                    |
# +------------------+--------------------------------------------------------------+
# |control_nreads    |Number of control reads in interval                           |
# +------------------+--------------------------------------------------------------+
# |control_peakcenter|Peak center of control                                        |
# +------------------+--------------------------------------------------------------+
# |control_peakval   |Number of reads at peak in control                            |
# +------------------+--------------------------------------------------------------+
# |end               |End coordinate of interval                                    |
# +------------------+--------------------------------------------------------------+
# |length            |Length of interval                                            |
# +------------------+--------------------------------------------------------------+
# |npeaks            |Number of peaks in interval                                   |
# +------------------+--------------------------------------------------------------+
# |nreads            |Number of reads in interval                                   |
# +------------------+--------------------------------------------------------------+
# |peakcenter        |Peak center in interval                                       |
# +------------------+--------------------------------------------------------------+
# |peakval           |Number of reads at peak                                       |
# +------------------+--------------------------------------------------------------+
# |start             |246251                                                        |
# +------------------+--------------------------------------------------------------+
# 
# The unprocessed output files created by the peak callers are in individual subdirectories for each caller (:file:`macs2.dir`, etc.).
# 
# IDR analysis
# ------------
# 
# The output of the IDR analysis is in the :file:`idr.dir` directory.
# 
# Example
# =======
# 
# Example data is available at http://www.cgat.org/~andreas/sample_data/pipeline_mapping.tgz.
# To run the example, simply unpack and untar::
# 
#    wget http://www.cgat.org/~andreas/sample_data/pipeline_mapping.tgz
#    tar -xvzf pipeline_mapping.tgz
#    cd pipeline_mapping
#    python <srcdir>/pipeline_mapping.py make full
# 
# .. note:: 
#    For the pipeline to run, install the :doc:`pipeline_annotations` as well.
# 
# .. _macs: http://liulab.dfci.harvard.edu/MACS/00README.html
# .. _spp: http://compbio.med.harvard.edu/Supplements/ChIP-seq/tutorial.html
# 
# 
# Glossary
# ========
# 
# .. glossary::
# 
# This pipeline implements the following nomenclature for peaks:
# 
#    region
#       A broad region defined by a peak caller. Regions are usually created in the first step of
#       peak calling, before peak refinement or detection of subpeaks is performed.
# 
#    summit
#       A narrow region defined by a caller. These are the output of any peak refinement
#       or subpeak detection steps
# 
#    interval
#       Generic term for an interval in which the read density is higher then expected. Both
#       a :term:`region` or a :term:`summit` are an :term:`interval`
#       
#    peak
#       Within a :term:`region` or :term:`summit` the position with the highest base coverage.
# 
# The pipeline computes some basic measures to validate peak calling. In order to fully annotate peaks, use :doc:`pipeline_intervals`.
# 
# .. note::
#    The pipeline currently expects that multi-mapping reads (reads mapping to multiple locations) have been removed.
# 
# 
# 
# QC
# ---
# 
# The pipeline implements the following QC measures. See :pmid:`22955991`.
# 
# NSC
#    normalized strand coefficient.
# 
# RSC
#   relative strand correlacion.
# 
# The pipeline will also do and IDR analysis (see https://sites.google.com/site/anshulkundaje/projects/idr) for spp called peaks. 

# <rawcell>

# Code
# ====
# 
# """

# <rawcell>

# Second step
# import necessary python packages, modules and pipelines 
# ruffus, numpy, sqlite, etc.
# PipelinePeakcalling,PipelineMotifs, etc.

# <codecell>

# load modules
from ruffus import *

import Experiment as E
import logging as L
import Database, CSV

import sys, os, re, shutil, itertools, math, glob, time, gzip, collections, random

import numpy, sqlite3
import GFF, GTF, Bed 
import IOTools
import IndexedFasta

import PipelinePeakcalling
import PipelineMotifs
import PipelineTracks
import PipelineMappingQC

# <rawcell>

# 
# First task:
# Pipeline configuration --> load options from the pipeline configuration file (*pipeline.ini*) file
# 
# Function is: P.getParameters()
# 
# Adjust paired end default when necessary
# 

# <codecell>

###################################################
###################################################
###################################################
## Pipeline configuration
###################################################

# load options from the config file
import Pipeline as P
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

# <rawcell>

# @follows
# Second task:
# Pre-process and set-up directories:
#     collect input files and identify experiments, conditions and replicates
#         Input files should be of the type:  <sample>-<condition>-<replicate>.genome.bam
#         Functions to get appropriate tracks are of the type: 
#             required functions are: getUnstimulated(track), getUnsubstracted(track), getSubstracted(track), getBamFiles(infile, suffix)
#             Others: getControl(track)
# 
#     prepare databases:
#         Function is: connect()
# 
#     make directories: mkdir("directory-name")
#     make connections to reporting and publishing directories and locations (read *sphinxreport.ini* and *conf.py* files)
#     

# <codecell>

###################################################################
###################################################################
## Helper functions mapping tracks to conditions, etc
###################################################################
import PipelineTracks

Sample = PipelineTracks.Sample3

# collect tracks, exclude any control tracks
TRACKS = PipelineTracks.Tracks( Sample ).loadFromDirectory( 
    [ x for x in glob.glob( "*.genome.bam" ) if PARAMS["tracks_control"] not in x ],
    "(\S+).genome.bam" )

ALL = Sample()
EXPERIMENTS = PipelineTracks.Aggregate( TRACKS, labels = ("condition", "tissue" ) )
CONDITIONS = PipelineTracks.Aggregate( TRACKS, labels = ("condition", ) )
TISSUES = PipelineTracks.Aggregate( TRACKS, labels = ("tissue", ) )

###################################################################
def getControl( track ):
    '''return appropriate control for a track
    '''
    n = track.clone()
    n.condition = PARAMS["tracks_control"]
    n.replicate = "R1"
    return n

def getUnstimulated( track ):
    '''return unstimulated condition for a track
    '''
    n = track.clone()
    n.condition = PARAMS["tracks_unstimulated"]
    return n

def getSubtracted( track ):
    '''return "subtracted" condition for a track
    '''
    n = track.clone()
    n.condition += PARAMS["tracks_subtract"]
    return n

def getUnsubtracted( track ):
    '''return "unsubtracted" condition for a track
    '''
    n = track.clone()
    if PARAMS["tracks_subtract"]:
        if n.condition.endswith( PARAMS["tracks_subtract"] ):
            n.condition = n.condition[:-len(PARAMS["tracks_subtract"])]
    return n

def getBamFiles( infile, suffix ):
    '''return associated bamfiles with infiles.'''
    track = P.snip( os.path.basename(infile), suffix)
    bamfile = P.snip( os.path.basename(infile), suffix ) + ".call.bam"
    assert os.path.exists( bamfile ), "bamfile %s does not exist" % bamfile

    controlfile = "%s.call.bam" % getControl(Sample(track)).asFile()
    if not os.path.exists( controlfile ):
        L.warn( "no controlfile '%s' for track '%s' not found " % (controlfile, track ) )
        controlfile = None

    return bamfile, controlfile

# <codecell>

###################################################################
###################################################################
###################################################################
def connect():
    '''connect to database.

    This method also attaches to helper databases.
    '''

    dbh = sqlite3.connect( PARAMS["database"] )
    statement = '''ATTACH DATABASE '%s' as annotations''' % (PARAMS["annotations_database"])
    cc = dbh.cursor()
    cc.execute( statement )
    cc.close()

    return dbh

# <codecell>

###################################################################
###################################################################
###################################################################
##
###################################################################
if os.path.exists("pipeline_conf.py"): 
    L.info( "reading additional configuration from pipeline_conf.py" )
    execfile("pipeline_conf.py")

# <rawcell>

# @follows
# @transforms
# Third task (? if required to process BAM files [normalise, mask regions] and check quality after this):
# Prepare input files
#     process BAM files if required (mask regions, de-duplicate, normalise, etc.) (unnecessary? Check options per peakcaller better):
#         Functions: makeMask(infile,outfile) #fileter reads if required
#                    prepareBAMforPeakCalling(infiles,outfile) #remove unmapped reads, deduplicate if necessary, etc.
#                    To normalise according to file with the lowest number of reads: countReadsInBAM(infile,outfile), minReads(infiles,outfile) #gets the lowest number of reads, normalizeBAM(infiles, outfile) # generates files with the same number of reads. 
#     
#     check quality of BAM data if necessary:
#         Function: loadDuplicationStats(infiles,outfile) # Maybe unnecessary as could check during pipeline_mapping.py
#                   buildBAMStats(infile, outfile)
#                   loadBAMStats(infiles, outfile)
#                   checkDataQuality(infile, outfile) # uses peakranger in pipeline_peakcalling.py; estimates noise rate
#                   loadDataQuality(infiles, outfile)
#                   checkLibraryComplexity(infile, outfile) # uses peakranger again; do these run without calling peaks in peakranger?

# <codecell>

############################################################
############################################################
############################################################
@files(input == None, "regions.mask")
def makeMask(infile,outfile):
    '''Make a mask for filtering reads if required.
    '''
    if PARAMS["calling_filter_exons"] or PARAMS["calling_filter_regions"]:
        regions_to_filter = []

        if PARAMS["calling_filter_exons"]:
            regions_to_filter += PipelinePeakcalling.getExonLocations(PARAMS["calling_filter_exons"]) 

        if PARAMS["calling_filter_regions"]:
            regions_to_filter += Bed.iterator( IOTools.openFile( PARAMS["calling_filter_regions"]) )
            
        fh = IOTools.openFile(outfile,"w")

        for bed in itertools.chain( regions_to_filter ):
            fh.write( "%s\n" % "\t".join(map(str, (bed.contig, bed.start, bed.end))) )

        fh.close()
    else:
        P.touch(outfile)

############################################################
############################################################
############################################################
@transform( "*.genome.bam", suffix(".genome.bam"), add_inputs( makeMask), ".prep.bam")
def prepareBAMForPeakCalling(infiles, outfile):
    '''Prepare BAM files for peak calling.

        - unmapped reads are removed.

        - if the option "calling_deduplicate" is Picard.MarkDuplicates is run 
            to remove duplicate reads

        - reads may be filtered by exon or location 

           - to remove reads by exon, the option "calling_filter_exons" should specify a file containing 
             a list of ensembl gene identifiers (one per line)
           - to remove reads by location, the option "calling_filter_regions" should specify a bed file''
        
        The resulting bam file has a .prep.bam extension. Merging infiles is currently untested and the 
        methods only consider single end reads.
    '''
    bam_file, mask_file = infiles

    if PARAMS["calling_filter_exons"] or PARAMS["calling_filter_regions"]:
        mask = mask_file
    else:
        mask = None

    PipelinePeakcalling.buildBAMforPeakCalling(bam_file,outfile,PARAMS["calling_deduplicate"], mask )

############################################################
############################################################
############################################################
@merge( prepareBAMForPeakCalling, "preparation_stats.load" )
def loadDuplicationStats( infiles, outfile ):
    '''load output from Picard Deduplication step.'''
    
    PipelineMappingQC.loadPicardMetrics( infiles, outfile, 
                                         pipeline_suffix = ".prep.bam",
                                         suffix = ".picard_metrics" )
    

############################################################
############################################################
############################################################
if PARAMS["calling_normalize"]==True:
    '''Normalise the number of reads in a set of prepared bam files.

    The minimum number of reads in a prepared bam file is calculated and this
    number of reads is used as a threshold to randomly sample from each bam file 
    in order to create a set of bam files with near identical numbers of reads.
  
    This may result in considerable data loss. 

    Per experimental contrast normalisation could be preferable.
   
    Potentially usefull if using a peak caller that does not correct for tag
    count between experimental and input samples.
    '''
    # First count the number of reads in each bam
    @transform(prepareBAMForPeakCalling,suffix("prep.bam"),"prep.count")
    def countReadsInBAM(infile,outfile):
        to_cluster = True
        statement= '''samtools idxstats %s | awk '{s+=$3} END {print s}' > %s ''' % ( infile,outfile )
        P.run()

    # Get the lowest number of reads
    @merge(countReadsInBAM,"minreads")
    def minReads(infiles,outfile):
        counts = []
        countfiles = glob.glob("*.prep.count")
        for cf in countfiles:
            fh = IOTools.openFile(cf,"r")
            counts.append(int(fh.read()))
            fh.close()
        fh = IOTools.openFile(outfile,"w")
        fh.write(str(min(counts)))
        fh.close

    # Make normalised files
    @follows(minReads)
    @transform(prepareBAMForPeakCalling,
               regex(r"(.*).prep.bam"),
               inputs( (r"\1.prep.bam",r"\1.prep.count") ), 
               r"\1.call.bam")
    def normalizeBAM( infiles, outfile ):
        '''build a normalized BAM file such that all
    files have approximately the same number of 
    reads.
    '''   
        fh = IOTools.openFile("minreads")
        minreads = int(fh.read())
        fh.close
        PipelinePeakcalling.buildSimpleNormalizedBAM( infiles, 
                                                      outfile,
                                                      minreads)
else:
    @transform(prepareBAMForPeakCalling,
               suffix(".prep.bam"),
               ".call.bam")
    def normalizeBAM( infile, outfile):
        P.clone( infile, outfile )
        P.clone( infile + ".bai", outfile + ".bai" )


# <codecell>

######################################################################
######################################################################
##                                                                  ##
##                 Statistics and QC Functions                      ## 
##                                                                  ##
######################################################################
######################################################################
@follows( mkdir("readstats.dir") )
@transform( normalizeBAM,
            regex("(.*).bam"),
            r"readstats.dir/\1.readstats" )
def buildBAMStats( infile, outfile ):
    '''count number of reads mapped, duplicates, etc.
    '''

    to_cluster = True

    statement = '''python
    %(scriptsdir)s/bam2stats.py
         --force
         --output-filename-pattern=%(outfile)s.%%s
    < %(infile)s
    > %(outfile)s
    '''

    P.run()

####################################################################
@merge( buildBAMStats, "bam_stats.load" )
def loadBAMStats( infiles, outfile ):
    '''import bam statisticis.'''
    
    PipelineMappingQC.loadBAMStats( infiles, outfile )


######################################################################
@follows( normalizeBAM )
@files( [ ("%s.call.bam" % (x.asFile()), 
           "%s.data_quality" % x.asFile() ) for x in TRACKS ] )
def checkDataQuality( infile, outfile ):
    '''uses peakranger to check data quality.

    nr is the signal/noise ratio.
    '''
    
    track = P.snip( infile, ".call.bam" )
    controlfile = "%s.call.bam" % getControl(Sample(track)).asFile()

    if not os.path.exists( controlfile ):
        L.warn( "controlfile '%s' for track '%s' not found " % (controlfile, track ) )
        P.touch( outfile )
        return

    to_cluster = True
    statement = '''peakranger nr --format bam %(infile)s %(controlfile)s 
                   | awk -v FS=":" '/Estimated noise rate/ { printf("estimated_noise_rate\\n%%f\\n", $2) }' > %(outfile)s'''
    P.run()

####################################################################
@merge( checkDataQuality, "data_quality.load" )
def loadDataQuality( infiles, outfile ):
    '''load data quality information.'''

    P.concatenateAndLoad( infiles, outfile, regex_filename = "(.*).data_quality" )

####################################################################
@follows( normalizeBAM )
@files( [ ("%s.call.bam" % (x.asFile()), 
           "%s.library_complexity" % x.asFile() ) for x in TRACKS ] )
def checkLibraryComplexity( infile, outfile ):
    '''uses peakranger to check library complexity.'''
    
    to_cluster = True
    statement = '''peakranger lc --format bam %(infile)s > %(outfile)s'''
    P.run()



# <rawcell>

# @follows
# @transforms
# Fourth task:
# # Use lax thresholds if IDR analysis will follow peak calling.
# Call peaks: call peaks, load results, summarise results, load summary, load specific peak caller functions (ie motif searching for GEM, quality calls with SPP, etc.).
# 
#     Run GEM
#     Run MACS2
#     Run GeneTrack (tag count estimation, minimum number of replicate reproducibility, peak pair distance estimation, peak calling, remove singletons, peak pairing)

# <codecell>

######################################################################
######################################################################
##                                                                  ##
##                    **** Call Peaks ****                          ## 
##                                                                  ##
######################################################################
######################################################################

######################################################################
######################################################################
##                                                                  ##
##                          GEM v1                                  ## 
##                                                                  ##
######################################################################
######################################################################

##gem options:
#Required:

#--d Read_Distribution_ChIP-exo.txt #path to the read distribution file 
#--exptcond1 /ifs/projects/antoniob/chIP/exoWTCHG/GEM-runs/250-252_cond1_rep1.bam #path to the aligned reads for experiment condition 1 replicate 1
#--exptcond1 /ifs/projects/antoniob/chIP/exoWTCHG/GEM-runs/250-252_cond1_rep2.bam #path to the aligned reads for experiment condition 1 replicate 2  
#(if replicates are to be combined at this point)

#Some options:
    
#--f SAM #Read file format (default is BED)
#--k_min 2 --k_max 40 --seed AAACTCTGTCTCAA 
#--g /ifs/projects/antoniob/chIP/exoWTCHG/GEM-runs/hg19.info 
#--genome /ifs/mirror/genomes/2bit/hg19.nib/ 
#--out 250-252-multi_4 >& stout_250-252_multi_4.out &

#Some optional flags:
#--outBED #Output BED files for genome browser
#--nrf #Do not filter duplicate reads (default is to apply filtering considering its neighboring base positions)
#--k_neg_dinu_shuffle #Use di-nucleotide shuffled sequences as negative sequences for motif finding
#--fa #GEM will use a fixed user-specified alpha value for all the regions

@follows( mkdir("gem.dir"), normalizeBAM )
@files( [ ("%s.call.bam" % (x.asFile()), 
           "gem.dir/%s.gem" % x.asFile() ) for x in TRACKS ] )
def callPeaksWithGEM( infile, outfile ):
    '''run GEM for peak detection.
    output bed files are compressed and indexed.
    '''
    track = P.snip( infile, ".call.bam" )

    controlfile = "%s.call.bam" % getControl(Sample(track)).asFile()

    if not os.path.exists( controlfile ):
        L.warn( "controlfile '%s' for track '%s' not found " % (controlfile, track ) )
        controlfile = None

    PipelinePeakcalling.runGEM( infile, outfile, controlfile)

############################################################
@transform( callPeaksWithGEM,
            regex(r"(.*).gem"),
            r"\1_gem.load" )
def loadGEM( infile, outfile ):
    '''load gem results.'''
    bamfile, controlfile = getBamFiles( infile, ".gem" )
    PipelinePeakcalling.loadGEM( infile, outfile, bamfile, controlfile )

############################################################
@merge( callPeaksWithGEM, "gem.summary" )
def summarizeGEM( infiles, outfile ):
    '''summarize GEM results.''' 
    PipelinePeakcalling.summarizeGEM( infiles, outfile )

############################################################
@merge( callPeaksWithGEM, "gem_fdr.summary" )
def summarizeGEMFDR( infiles, outfile ):
    '''summarize GEM results.''' 
    PipelinePeakcalling.summarizeGEMFDR( infiles, outfile )

############################################################
@transform( summarizeGEM,
            suffix(".summary"),
            "_summary.load" )
def loadGEMSummary( infile, outfile ):
    '''load gem summary.'''
    P.load( infile, outfile, "--index=track" )

############################################################
@transform( summarizeGEMFDR,
            suffix("_fdr.summary"),
            "_fdr.load" )
def loadGEMSummaryFDR( infile, outfile ):
    '''load gem summary.'''
    P.load( infile, outfile, "--index=track", transpose="fdr" )



# <codecell>

######################################################################
######################################################################
##                                                                  ##
##                          MACS version 2                          ## 
##                                                                  ##
######################################################################
######################################################################

@follows( mkdir("macs2.dir"), normalizeBAM )
@files( [ ("%s.call.bam" % (x.asFile()), 
           "macs2.dir/%s.macs2" % x.asFile() ) for x in TRACKS ] )
def callPeaksWithMACS2( infile, outfile ):
    '''run MACS 2 for peak detection.
    output bed files are compressed and indexed.
    '''
    track = P.snip( infile, ".call.bam" )

    controlfile = "%s.call.bam" % getControl(Sample(track)).asFile()

    if not os.path.exists( controlfile ):
        L.warn( "controlfile '%s' for track '%s' not found " % (controlfile, track ) )
        controlfile = None

    PipelinePeakcalling.runMACS2( infile, outfile, controlfile)

############################################################
@transform( callPeaksWithMACS2,
            regex(r"(.*).macs2"),
            r"\1_macs2.load" )
def loadMACS2( infile, outfile ):
    '''load macs results.'''
    bamfile, controlfile = getBamFiles( infile, ".macs2" )
    PipelinePeakcalling.loadMACS2( infile, outfile, bamfile, controlfile )

############################################################
@merge( callPeaksWithMACS2, "macs2.summary" )
def summarizeMACS2( infiles, outfile ):
    '''summarize MACS results.''' 
    PipelinePeakcalling.summarizeMACS2( infiles, outfile )

############################################################
@merge( callPeaksWithMACS2, "macs2_fdr.summary" )
def summarizeMACS2FDR( infiles, outfile ):
    '''summarize MACS results.''' 
    PipelinePeakcalling.summarizeMACS2FDR( infiles, outfile )

############################################################
@transform( summarizeMACS2,
            suffix(".summary"),
            "_summary.load" )
def loadMACS2Summary( infile, outfile ):
    '''load macs2 summary.'''
    P.load( infile, outfile, "--index=track" )

############################################################
@transform( summarizeMACS2FDR,
            suffix("_fdr.summary"),
            "_fdr.load" )
def loadMACS2SummaryFDR( infile, outfile ):
    '''load macs2 summary.'''
    P.load( infile, outfile, "--index=track", transpose="fdr" )


# <codecell>

######################################################################
######################################################################
##                                                                  ##
##                              SPP                                 ## 
##                                                                  ##
######################################################################
######################################################################

@follows( mkdir("spp.dir"), normalizeBAM )
@files( [ ("%s.call.bam" % (x.asFile()), 
           "spp.dir/%s.spp" % x.asFile() ) for x in TRACKS ] )
def callPeaksWithSPP( infile, outfile ):
    '''run SPP for peak detection.'''
    track = P.snip( infile, ".call.bam" )
    controlfile = "%s.call.bam" % getControl(Sample(track)).asFile()

    if not os.path.exists( controlfile ):
        L.warn( "no controlfile '%s' for track '%s' not found " % (controlfile, track ) )
        controlfile = None

    PipelinePeakcalling.runSPP( infile, outfile, controlfile)

############################################################
@transform( callPeaksWithSPP,
            regex(r"(.*).spp"),
            r"\1_spp.load" )
def loadSPP( infile, outfile ):
    '''load spp results.'''
    bamfile, controlfile = getBamFiles( infile, ".spp" )
    PipelinePeakcalling.loadSPP( infile, outfile, bamfile, controlfile )

############################################################
@merge( callPeaksWithSPP, "spp.summary" )
def summarizeSPP( infiles, outfile ):
    '''summarize SPP results.'''
    PipelinePeakcalling.summarizeSPP( infiles, outfile )

############################################################
@transform( summarizeSPP, suffix(".summary"), "_summary.load" )
def loadSPPSummary( infile, outfile ):
    '''load SPP summary.'''
    P.load( infile, outfile, "--index=track" )

############################################################
############################################################
@follows( mkdir( os.path.join( PARAMS["exportdir"], "quality" ) ),
          mkdir( "spp.dir" ),
          normalizeBAM )
@files( [ ("%s.call.bam" % (x.asFile()), 
           "spp.dir/%s.qual" % x.asFile() ) for x in TRACKS ] )
def estimateSPPQualityMetrics( infile, outfile ):

    track = P.snip(infile, ".bam" )

    controlfile = "%s.call.bam" % getControl(Sample(track)).asFile()

    if not os.path.exists( controlfile ):
        L.warn( "no controlfile '%s' for track '%s' not found " % (controlfile, track ) )
        raise ValueError( "idr analysis requires a control")
    
    executable = P.which( "run_spp.R" )
    if executable == None:
        raise ValueError( "could not find run_spp.R" )

    statement = '''
    Rscript %(executable)s -c=%(infile)s -i=%(controlfile)s -rf -savp -out=%(outfile)s
    >& %(outfile)s.log'''
    
    P.run()

    if os.path.exists( track + ".pdf" ):
        dest = os.path.join( PARAMS["exportdir"], "quality", track + ".pdf" )
        if os.path.exists( dest ): os.unlink( dest )
        shutil.move( track + ".pdf", dest )

############################################################
############################################################
############################################################
@merge( estimateSPPQualityMetrics, "spp_quality.load" )
def loadSPPQualityMetrics( infiles, outfile ):
    '''load spp quality metrics.'''
    P.concatenateAndLoad( infiles, outfile,
                          regex_filename = "spp.dir/(.*).qual",
                          header = "track,bamfile,mapped_reads,estFragLen,corr_estFragLen,phantomPeak,corr_phantomPeak,argmin_corr,min_corr,nsc,rsc,quality")
    


# <rawcell>

# @follows
# @transforms
# Fifth task:
# Run quality control:
#     Total peaks called; read and peak distribution across the genome (Landt et al. 2012 F3A)
# 
#     Reproducibility: IDR R package # Plot Landt el al 2012 F6A, B and C; Run IDR on pseudoreplicates if no true replicates available or not appropriate; ENCODE cut-off = number of bound regions identified in IDR comparison between replicates must be at least 50% of the number of regions identified in an IDR comparison between pseudoreplicates of random partition of reads from all replicates. Apply flowchart Landt et al 2012 F7D to eliminate low quality samples.
# 
#     Background estimation: FRiP #fraction of reads in peaks (Landt et al 2012 F4C, chIP-seq cut-off >1%). 
#     Sequencing depth: call Saturation # Fraction of total peaks called vs. number of mapped reads; peak median enrichment vs. number of mapped reads (Landt et al. 2012 F3B, F3C)
# 
#     Library complexity estimation: load from mapping statistics (Picard tools) (=fraction of DNA fragments that are non-redundant [Landt et al. 2012, F4A]); show IGV shot of mapped reads; calculate the non redundant fraciont (NRF; =ratio between the number of positions in the genome that uniquely mappable reads map to and the total number of uniquely mappable reads. Aim for an NRF of 0.8 for 10 million reads [(Landt et al. 2012, F7C, E and F])
# 
#     PCR artefacts estimation: tool?
# 
#     Motif enrichment: number of peaks with recognisable motif, plot Landt et al 2012 F5F; check Johnson, Mortazavi et al Science 2007.
#     

# <codecell>

######################################################################
######################################################################
##                                                                  ##
##        IDR Analysis (using relaxed SPP peak calling)             ## 
##                                                                  ##
######################################################################
######################################################################
@follows( mkdir("idr.dir"), normalizeBAM )
@files( [ ("%s.call.bam" % (x.asFile()), 
           "idr.dir/%s.spp" % x.asFile() ) for x in TRACKS ] )
def callPeaksWithSPPForIDR( infile, outfile ):
    '''run SICER for peak detection.'''
    track = P.snip( infile, ".call.bam" )
    controlfile = "%s.call.bam" % getControl(Sample(track)).asFile()

    if not os.path.exists( controlfile ):
        L.warn( "no controlfile '%s' for track '%s' not found " % (controlfile, track ) )
        raise ValueError( "idr analysis requires a control")

    executable = P.which( "run_spp.R" )
    if executable == None:
        raise ValueError( "could not find run_spp.R" )

    statement = '''
    Rscript %(executable)s -c=%(infile)s -i=%(controlfile)s -npeak=%(idr_npeaks)s 
            -odir=idr.dir -savr -savp -rf -out=%(outfile)s
    >& %(outfile)s.log'''
    
    P.run()

    track = P.snip(infile, ".bam" )

    if os.path.exists( track + ".pdf" ):
        shutil.move( infile + ".pdf", os.path.join( PARAMS["exportdir"], "idr" ))

############################################################
@collate( callPeaksWithSPPForIDR, 
          regex( r"idr.dir/(.+)-[^-]+.spp" ),
          r"idr.dir/\1.idr")
def applyIDR( infiles, outfile ):
    '''apply IDR analysis.'''

    to_cluster = True
    chromosome_table = os.path.join(PARAMS["annotations_dir"], PARAMS_ANNOTATIONS["interface_contigs"])

    for infile1, infile2 in itertools.combinations( infiles, 2 ):
        E.info( "applyIDR: processing %s and %s" % (infile1,infile2))

        basename1 = os.path.basename( infile1 )
        basename2 = os.path.basename( infile2 )

        track1 = P.snip( basename1, ".spp" )
        control1 = getControl(Sample(track1)).asFile()
        track2 = P.snip( basename2, ".spp" )
        control2 = getControl(Sample(track2)).asFile()
        
        statement = '''
          python %(scriptsdir)s/WrapperIDR.py 
                 --action=run
                 --output-prefix=%(track1)s_vs_%(track2)s.idr
                 --chromosome-table=%(chromosome_table)s
                 idr.dir/%(track1)s.call_VS_%(control1)s.call.regionPeak.gz 
                 idr.dir/%(track2)s.call_VS_%(control2)s.call.regionPeak.gz 
          >> %(outfile)s'''

        P.run()

############################################################
@follows(mkdir( os.path.join( PARAMS["exportdir"], "idr" ) ) )
@transform( applyIDR, 
            suffix(".idr"),
            ".plot")
def plotIDR( infile, outfile ):
    '''plot IDR results.'''

    to_cluster = True

    track = P.snip( infile, ".idr")
    files = glob.glob( track + "*.idr-em.sav" )
    files = " ".join([ P.snip(x, "-em.sav" ) for x in files ])
    output_prefix = os.path.join( PARAMS["exportdir"], "idr", os.path.basename(track ) )
    statement = '''
    python %(scriptsdir)s/WrapperIDR.py
               --action=plot
               --output-prefix=%(output_prefix)s
               %(files)s
    > %(outfile)s'''

    P.run()

############################################################



# <rawcell>

# @follows
# @transforms
# Sixth task:
# Collate peak calling results
#     Generate bed files
#     Peaks: number, height?, location, p-values, FDR values
#     Test/Report reproducibility between replicates

# <codecell>

######################################################################
######################################################################
##                                                                  ##
##            ___ Collate peak calling results ___                  ## 
##                                                                  ##
######################################################################
######################################################################

CALLINGTARGETS, SUMMARYTARGETS = [], []
mapToCallingTargets = { 'macs2' : loadMACS2,
                        'spp': loadSPP,
			'gem': loadGEM,
			'GeneTrack': loadGeneTrack,
                        }

mapToSummaryTargets = { 'macs2': [loadMACS2Summary, loadMACS2SummaryFDR],
                        'spp' : [loadSPPSummary],
                        'gem': [loadGEMSummary, loadGEMSummaryFDR],
			'GeneTrack': [loadGeneTrackSummary, loadGeneTrackSummaryFDR],
                        }

for x in P.asList( PARAMS["peakcallers"]):
    CALLINGTARGETS.append( mapToCallingTargets[x] )
    if x in mapToSummaryTargets:
        SUMMARYTARGETS.extend( mapToSummaryTargets[x] )

###################################################################
###################################################################
###################################################################
@follows( *(CALLINGTARGETS + SUMMARYTARGETS) )
def calling(): pass

############################################################
############################################################
############################################################
@follows( mkdir( os.path.join( PARAMS["exportdir"], "bedfiles" ) ) )
@transform( CALLINGTARGETS, regex("(.*)/(.*).load"), 
            (os.path.join( PARAMS["exportdir"], "bedfiles", r"\2.peaks.bed.gz" ),
             os.path.join( PARAMS["exportdir"], "bedfiles", r"\2.regions.bed.gz" ),
             os.path.join( PARAMS["exportdir"], "bedfiles", r"\2.summits.bed.gz" ) ) )
def exportIntervalsAsBed( infile, outfiles ):
    '''export all intervals as bed files.'''

    outfile_peaks, outfile_regions, outfile_summits = outfiles
    track = P.snip( os.path.basename(infile), ".load" ) 

    #PipelinePeakcalling.exportIntervalsAsBed( infile, outfile_regions, "%s_regions" % P.quote(track) )

    dbh = connect()
    tablename = "%s_peaks" % P.quote(track) 
    if tablename in Database.getTables( dbh ):
        PipelinePeakcalling.exportIntervalsAsBed( infile, outfile_peaks, tablename )
    else:
        E.warn( "no table %s - empty bed file output" % tablename )
        P.touch( outfile_peaks )

    dbh = connect()
    tablename = "%s_regions" % P.quote(track) 
    if tablename in Database.getTables( dbh ):
        PipelinePeakcalling.exportIntervalsAsBed( infile, outfile_regions, tablename )
    else:
        E.warn( "no table %s - empty bed file output" % tablename )
        P.touch( outfile_regions )

    dbh = connect()
    tablename = "%s_summits" % P.quote(track) 
    if tablename in Database.getTables( dbh ):
        PipelinePeakcalling.exportIntervalsAsBed( infile, outfile_summits, tablename )
    else:
        E.warn( "no table %s - empty bed file output" % tablename )
        P.touch( outfile_summits )

###################################################################
###################################################################
###################################################################
# Targets for the annotation of intervals.
###################################################################
@split( exportIntervalsAsBed, os.path.join( PARAMS["exportdir"], "bedfiles", "*.bed.gz" ) )
def flattenBedFiles( infile, outfile ):
    '''dummy target - merge all files in exportIntervalsAsBed'''


def getPeakShift( track, method ):
    '''return peak shift for track and method.'''
    dbh = connect()
    result = Database.executewait( dbh, "SELECT shift FROM %(method)s_summary where track = '%(track)s'" % locals())
    return result.fetchone()[0]

###################################################################
###################################################################
###################################################################
@follows( mkdir( "peakshapes.dir" ) )
@transform( flattenBedFiles,
            regex(".*/(.*).bed.gz"),
            r"peakshapes.dir/\1.peakshape.tsv.gz" )
def buildPeakShapeTable( infile, outfile ):
    '''build a table with peak shape parameters.'''
    
    to_cluster = True

    # compute suffix (includes method name)
    track, method, section = re.match( "(.*)_(.*)\.(.*).bed.gz", os.path.basename(infile) ).groups()
    
    suffix = "_%s.%s.bed.gz" % (method, section)
    bamfile, controlfile = getBamFiles( infile, suffix )
   
    shift = getPeakShift( track, method )

    if shift:
        E.info( "applying read shift %i for track %s" % (shift, track ) )

    options = []
    if controlfile:
        options.append( "--control-file=%s" % controlfile )
    options = " ".join( options )

    statement = '''python %(scriptsdir)s/bam2peakshape.py
                      --window-size=%(peakshape_window_size)i
                      --bin-size=%(peakshape_bin_size)i
                      --output-filename-pattern="%(outfile)s.%%s"
                      --force
                      --shift=%(shift)i
                      --sort=peak-height
                      --sort=peak-width
                      %(options)s
                      --log=%(outfile)s.log
                      %(bamfile)s %(infile)s
                   | gzip
                   > %(outfile)s
                '''
    P.run()

###################################################################
@transform( buildPeakShapeTable, suffix(".tsv.gz"), ".load" )
def loadPeakShapeTable( infile, outfile ):
    '''load peak shape information.'''
    P.load( infile, outfile, "--ignore-column=bins --ignore-column=counts --allow-empty" )

############################################################
############################################################
############################################################
## targets to do with the analysis of replicates
############################################################
# dummy task - flatten the nested list of result files created
# by exportIntervalsAsBed
@split( exportIntervalsAsBed, os.path.join( PARAMS["exportdir"], "bedfiles", "*.bed.gz" ) )
def allIntervalsAsBed( infile, outfile ): pass

############################################################
############################################################
############################################################
@follows( mkdir( "reproducibility.dir"))
@collate( allIntervalsAsBed,
          regex( os.path.join( PARAMS["exportdir"], "bedfiles", r"(.+)_(.+)\.(.+).bed.gz" )),
          r"reproducibility.dir/\1.\3.reproducibility")
def makeReproducibilityOfMethods( infiles, outfile ):
    '''compute overlap between intervals.

    Note the exon percentages are approximations assuming that there are
    not more than one intervals in one set overlapping one in the other set.
    '''
    PipelinePeakcalling.makeReproducibility( infiles, outfile )

############################################################
############################################################
############################################################
@follows( mkdir( "reproducibility.dir"))
@collate( allIntervalsAsBed, 
          regex( os.path.join( PARAMS["exportdir"], "bedfiles", r"(.+)-[^-]+_(.+)\.(.+).bed.gz" )),
          r"reproducibility.dir/\1-\2.\3.reproducibility")
def makeReproducibilityOfReplicates( infiles, outfile ):
    '''compute overlap between intervals.

    Note the exon percentages are approximations assuming that there are
    not more than one intervals in one set overlapping one in the other set.
    '''
    PipelinePeakcalling.makeReproducibility( infiles, outfile )

############################################################
############################################################
############################################################
@transform( (makeReproducibilityOfMethods, makeReproducibilityOfReplicates), suffix(".reproducibility"), "_reproducibility.load" )
def loadReproducibility( infile, outfile ):
    '''load Reproducibility results
    '''
    P.load( infile, outfile, options="--allow-empty" )

############################################################
############################################################
############################################################
@follows( loadReproducibility )
def reproducibility(): pass
    
###################################################################
@follows( loadBAMStats, loadDuplicationStats, loadSPPQualityMetrics )
def qc(): pass

###################################################################
@follows( calling, exportIntervalsAsBed, qc )
def full(): pass


# <rawcell>

# @follows
# Seventh task:
#     Compare ChIP-seq to ChIP-Exo:
#         Number of shared peaks, differences in location, peak height (fold enrichment)
#         Background noise comparison, read depth required, resolution, motifs found
#         

# <codecell>


# <rawcell>

# @follows
# Eighth task:
# Create report and publish (Landt et al 2012 Reporting guidelines (Box 4)):
#     
# Metadata
#     Investigator, organism, or cell line, experimental protocol (or reference to a known protocol).
#     Indication as to whether an experiment is a technical or biological replicate.
#     Catalog and lot number for any antibody used. If not a commercial antibody, indicate the precise source of the antibody.
#     Information used to characterize the antibody, including summary of results (images of immunoblots, immunofluorescence, list of proteins identified by mass spec, etc.).
#     Peak calling algorithm and parameters used, including threshold and reference genome used to map peaks.
#     A summary of the number of reads and number of targets for each replicate and for the merged data set.
#     Criteria that were used to validate the quality of the resultant data (i.e., overlap results or IDR).
#     Experimental validation results (e.g., qPCR).
#     Link to the control track that was used (if applicable).
#     An explanation if the experiment fails to meet any of the standards.
# 
# High-throughput sequencing data
#     Image files from sequencing experiments do not need to be stored.
#     Raw data (FASTQ files) should be submitted to both GEO and SRA.
#     Each replicate should be submitted independently.
#     Target region and peak calling results to GEO.
# 
# Point source peaks
# 
#     Peak position, defined as a single base pair.
#     Start and end positions, defined as specific base pairs.
#     Signal value (e.g., fold enrichment) using an algorithm chosen by the submitter.
#     Significance/accuracy measures:
#         ○ P-value determined using a method chosen by the submitter.
#         ○ Q-value (false discovery rate correction) determined using a method chosen by the submitter.
#     Other methods and data as applicable

# <codecell>

###################################################################
###################################################################
###################################################################
@follows( mkdir( "report" ) )
def build_report():
    '''build report from scratch.'''

    E.info( "starting documentation build process from scratch" )
    P.run_report( clean = True )

###################################################################
###################################################################
###################################################################
@follows( mkdir( "report" ) )
def update_report():
    '''update report.'''

    E.info( "updating documentation" )
    P.run_report( clean = False )

###################################################################
###################################################################
###################################################################
@follows( mkdir( "%s/bamfiles" % PARAMS["web_dir"]), 
          mkdir("%s/genesets" % PARAMS["web_dir"]),
          mkdir("%s/classification" % PARAMS["web_dir"]),
          mkdir("%s/differential_expression" % PARAMS["web_dir"]),
          update_report,
          )
def publish():
    '''publish files.'''
    # publish web pages
    P.publish_report()

    # publish additional data
    web_dir = PARAMS["web_dir"]
    project_id = P.getProjectId()

    # directory, files
    exportfiles = {
        "bamfiles" : glob.glob( "*.accepted.bam" ) + glob.glob( "*.accepted.bam.bai" ),
        "genesets": [ "lincrna.gtf.gz", "abinitio.gtf.gz" ],
        "classification": glob.glob("*.class.tsv.gz") ,
        "differential_expression" : glob.glob( "*.cuffdiff.dir" ),
        }
    
    bams = []

    for targetdir, filenames in exportfiles.iteritems():
        for src in filenames:
            dest = "%s/%s/%s" % (web_dir, targetdir, src)
            if dest.endswith( ".bam"): bams.append( dest )
            dest = os.path.abspath( dest )
            if not os.path.exists( dest ):
                os.symlink( os.path.abspath(src), dest )
    
    # output ucsc links
    for bam in bams: 
        filename = os.path.basename( bam )
        track = P.snip( filename, ".bam" )
        print """track type=bam name="%(track)s" bigDataUrl=http://www.cgat.org/downloads/%(project_id)s/bamfiles/%(filename)s""" % locals()
        

# <rawcell>

# Final step:
# exit

# <codecell>

if __name__== "__main__":
    sys.exit( P.main(sys.argv) )

# <rawcell>

# Next steps:
#     
# Motif finding --> Annotation and integration --> Study specific questions (allele specific binding)

# <codecell>



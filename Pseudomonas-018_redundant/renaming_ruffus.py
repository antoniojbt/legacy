###################################################
##Recursive renaming 018 Pseudomonas and CF files##
###################################################

from ruffus import *

import sys, os, csv, glob

import CGAT.IOTools as IOTools

import CGAT.Pipeline as P

import CGAT.Experiment as E


@transform(glob.glob('*/*.gz'), regex(r"(.*)/(.*).gz"), r"\1\2.gz")

def rename(infile, outfile): 
    to_cluster = False

    statement = ''' cd %(infile)s;
                    rename 's/130906_0085_000000000-A5A34_1_SA-PE-0.._/%(outfile)s-1-1\.fastq\./' *;
                    rename sanfastq.gz gz *;
                    			
		''' % locals()
    P.run()

if __name__=="__main__":

    sys.exit(P.main(sys.argv))

###################################################

import fnmatch
import os
import re


def listFiles(dir):
    
    rootdir = dir
    
    rootdir = '/ifs/projects/sftp/backup/proj018/2013096_Bundy_Jake/raw_reads_renamed/'
    for root, subFolders, files in os.walk(rootdir):
        for subFolder in subFolders:
            print subFolder
            
#        for file in files:
            yield os.path.join(root,file)
    return

if __name__=="__main__":

    sys.exit(sys.argv)


#print listFiles('/ifs/projects/sftp/backup/proj018/2013096_Bundy_Jake/raw_reads_renamed')

for f in listFiles(r'/ifs/projects/sftp/backup/proj018/2013096_Bundy_Jake/raw_reads_renamed'):
    read_pair = re.search('PE-0(..)_(.)', filename)
    os.rename(filename, 'sample_'+ dirnames + '-condition_1-replicate_1.fastq.'+ read_pair.group(1) + '.gz')


    #if f[-5].isalpha():
    #    os.rename(f,f[:-5]+f[-5].lower() + ".JPG")
    #    print "Renamed " +  "---to---" + f[:-5]+f[-5].lower() + ".JPG"  


###################################################

import re
import os

# Top directory to start the recursive renaming
DIR = '/ifs/projects/sftp/backup/proj018/2013096_Bundy_Jake/raw_reads_renamed/'

# Regex to search for
rexpr   = '([\S]*)[.]middle[.]txt([\S]*)'
pattern = re.compile(rexpr)

# Short function to rename all files
def rename_files(dir):
    for fname in os.listdir(dir):
        match = pattern.search(fname)
        if match:
            oldname = os.path.join( dir, fname )
            newname = os.path.join( dir, match.group(1) + '.eph' + match.group(2) )
            print oldname + ' => ' + newname
            #os.rename( oldname, newname )

for root, dirs, files in os.walk(DIR):
    rename_files(root)

###################################################

import os

def fileRenaming(dir):

    raw_input = 'Directory path?'

    dir = raw_input
    subFolders = os.walk(rootdir)
    for subFolder in subFolders:
        print subFolder
   
call to main 

if __name__ == "__main__":
    sys.exit( main( sys.argv) )

import sys
import os
import re

top = '/ifs/projects/sftp/backup/proj018/2013096_Bundy_Jake/raw_reads_renamed'

def renaming(top):
    temp = os.walk(top, topdown=False)
    for dirpath, dirnames, filenames in temp:
        for f in filenames:
            f = os.path.join(dirnames, filenames)
            read_pair_number = re.search('PE-0(..)_(.)', f)
            read_pair = read_pair_number.group(1)
            new_name = str('sample_'+ dirnames + '-condition_1-replicate_1.fastq.'+ read_pair.group(1) + '.gz')
            os.rename(f, new_name)
    return

print renaming(top)


#Old file:
#raw_reads_renamed/'dir'/130906_0085_000000000-A5A34_1_SA-PE-00?_?.sanfastq.gz 

#New file:
#raw_reads_renamed/'dir'/sample_'dir'-condition_1-replicate_1.fastq.#read_pair.gz

#if __name__=="__main__":
#    sys.exit()
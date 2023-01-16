"""
Created on Tue Oct 23 17:41:59 2018

@author: mqbpkjs2
reated on Tue Oct 23 17:41:59 2018

@author: mqbpkjs2
"""

import csv
import sys
import os

script, sampleTable,cores = sys.argv

def getFastqFiles(sampleTable,cores):
    
    with open(sampleTable) as tsvfile:
        reader = csv.DictReader(tsvfile, delimiter='\t')
        for row in reader:
            fileNames=str.split(row["File"],"|")
    
            os.system('vdb-config -s "/repository/user/main/public/root=./"')
            

            for fileName in fileNames:   
                command = "prefetch -O ./ " + fileName
                os.system(command)                
                fileName = fileName + ".sra"
                command = "parallel-fastq-dump --split-3 --gzip --sra-id {0} --threads {1} --tmpdir ./".format(fileName, cores)
                os.system(command)
                os.system("rm " +fileName)

                   
                    
if __name__ == "__main__":
    getFastqFiles(sampleTable,cores)



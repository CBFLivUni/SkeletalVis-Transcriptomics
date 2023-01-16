#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import sys
import os

sys.argv.pop(0)
single = sys.argv.pop(0)
samples = sys.argv


if single.lower() == "true":
    for sample in samples:
        command = "trimmomatic SE -threads 4 -phred33 {0} trimmed_{0} LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:25".format(sample)
        os.system(command)
else:
    for i in range(0,len(samples)-1,2):
        sample1 = samples[i]
        sample2 = samples[i+1]
        command = "trimmomatic PE -threads 4 -phred33 {0} {1} trimmed_{0} unpaired_{0} trimmed_{1} unpaired_{1} LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:25".format(sample1, sample2)
        os.system(command)
        
        

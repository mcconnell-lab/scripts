#!/usr/bin/env python
# encoding: utf-8
"""
dupeStats.py
Created by Michael Lindberg on 2010-12-16.
Modified by Jim Havrilla & Alex Koeppel on 2014-01-11
Minimal annotations by Will Chronister 2019-09-06
"""

#Calculates duplicate read statistics from flagstat files in main pipeline:
#sts = call("python " + workingpath + "dupeStats_MOD.py -i " + path + bamhead + "_sorted.flagstat -d" +  path + bamhead + "_sorted_rmdup.flagstat >> " + output +"/sample_stats/" + "All_FlagStats.txt", shell=True)

import sys
import os
from optparse import OptionParser
import math
import re


#Jim's regex version of the function
def func(normFile, dupeFile):
	norm = open(normFile, 'r')
	for line in norm:
		m=re.match('(\d*)\s*\S.*?(total|properly paired|read1)',line)
		n=re.match('(\d*)\s*.*?(mapped)\s*\((\d*\.\d*)%.*\)',line)
		if m!=None:
#			print m.group(0)
			if m.group(2)=="total":
				totalReads=float(m.group(1)) 
			if m.group(2)=="read1":
				pairs=float(m.group(1))
			if m.group(2)=="properly paired":
				conc=m.group(1)
 		if n!=None:
			mappedReads=float(n.group(1)) 
			percentParse = n.group(3)
	dupe = open(dupeFile, 'r')
	for line in dupe:
		n=re.match('(\d*)\s*.*?(mapped)\s*\((\d*\.\d*)%.*\)',line)
		if n!=None:
			deduped=float(n.group(1)) 
	dedupeRate = (1-(deduped/mappedReads))*100    
	print normFile + "\t" + ("%d\t%d\t%d\t%s" % (pairs,totalReads,mappedReads,percentParse)) + "%\t" + conc + "\t" + ("%f" % (dedupeRate)) + "%\t" + ("%d" % (deduped))
    
def main():
    usage = \
    """%prog -i /<prededupe flagstat>/ -d /<post dedupe flagstat>/

    """
    parser = OptionParser(usage)
    
    parser.add_option("-i", dest="normFile", help="Flagstat of predudupe", metavar="FILE")
    parser.add_option("-d", dest="dupeFile", help="Flagstat of dedupe", metavar="FILE")

    # Grab the command line options
    (opts, args) = parser.parse_args()
    
    # Get to work.
    if opts.normFile and opts.dupeFile:
        func(opts.normFile, opts.dupeFile)
    else:
        parser.print_help()
        sys.exit(2)

if __name__ == '__main__':
    main()


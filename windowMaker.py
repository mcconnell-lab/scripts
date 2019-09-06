#!/usr/local/bin/python2.6
# encoding: utf-8
"""
windowMaker.py
Created by Michael Lindberg on 2011-03-10.
Minimal annotations by Will Chronister, 2019-09-06.
"""
#This script is invoked from main CNV pipeline with following command:
#python /path/to/windowMaker.py -c /path/to/human_g1k_v37.chrlist -n 500000 -f /path/to/b37_full.maskedMap.40mer.fa > all.window.40mer.b37.map.500kb.bed
#Creates genomic bins/windows of specified length of mappable sequence
#Requires already masked FASTA genome (masking based on 40mer mappability track)


import sys
import os
from optparse import OptionParser
import pysam

#read in the FASTA Genome
def readFA(fastaFile, chrFile, windowSize):
	chrList = []
	for line in open(chrFile, 'r'):
		chrome = line.strip('\n')
		if len(chrome) > 0:
			chrList.append(chrome)
	fafile = pysam.Fastafile(fastaFile)
	for el in chrList:
		if fafile.fetch(region=el):
			currChr = fafile.fetch(region=el)
			calcWindow(currChr, windowSize, el)

def calcWindow(chrom, windowMax, chrName):
	currPos = 0
	GCCount = 0
	currWindow = 0
	NCount = 0
	startPos = 0
	for base in chrom:
		if (currWindow == 0 and NCount == 0):
			startPos = currPos
		currPos += 1
		if (currWindow <= int(windowMax) - 1):
			if base == 'N':
				NCount += 1
			else:
				if base == 'G' or base == 'C' or base == 'g' or base == 'c':
					GCCount += 1
				currWindow += 1

		if currWindow == int(windowMax):
			writeBed(chrName, startPos, currPos, GCCount, NCount, windowMax)
			NCount = 0
			GCCount = 0
			currWindow = 0
	if currWindow < int(windowMax):
		writeBed(chrName, startPos, currPos, GCCount, NCount, windowMax)
	

def writeBed(chrName, startPos, currPos, GCCount, NCount, windowMax):
	print chrName + "\t" + str(startPos) + "\t" + str(currPos) + "\t" +  str(float(GCCount)/float(windowMax)) + "\t" + str(float(NCount)/float(NCount+int(windowMax)))


		


def main():
	usage = """%prog -f <fasta genome file> -n <window size>\n
	Calculates window sizes masking repeats and Ns on reference genome while reporting GC content and N content"""

	parser = OptionParser(usage)
	parser.add_option("-f", "--fasta", dest="fastaFile",
	                                help="fasta input file", metavar="FILE")
	parser.add_option("-n", "--window", dest="windowSize",
	                                help="window size", default=1000, metavar="INT")	
	parser.add_option("-c", "--chroms", dest="chrFile",
			                help="chrome list to sample windows from, single column from identifiers in genome fasta", metavar="FILE")

			  
	#Get to work
	(options, args) = parser.parse_args()
	
	if options.fastaFile:
		readFA(options.fastaFile, options.chrFile, options.windowSize)
	
	else:
	        parser.print_help()
	 #       no gib input, output the help menu
	
if __name__ == "__main__":
	main()

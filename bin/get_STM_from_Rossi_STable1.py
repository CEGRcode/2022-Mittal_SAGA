import sys, getopt
#import pysam
import re
import numpy as np
import argparse

usage = """
Usage:
This script processes STM occupancy reference out of Rossi STable1.

Example: python
"""

def getParams():
	'''Parse parameters from the command line'''
	parser = argparse.ArgumentParser(description='This script processes STM occupancy reference out of Rossi STable1.')
	parser.add_argument('-i','--input', metavar='table', required=True, help='The Rossi 2021 Supplemenatry Table 1 file as a tab-delimited flat text file')
	parser.add_argument('-p','--peaks', metavar='peaks', required=True, help='Sua7 peaks as backup ref')
	parser.add_argument('-o','--output', metavar='outbed', required=True, help='bed output')
	parser.add_argument('-s','--shift', metavar='numbp', dest='shift', type=int, default=0, help='upstream bp shift when imputing default from TSS')

	args = parser.parse_args()
	return(args)

def load_peaks(peakfile):
	'''Load peaks keyed on id'''
	peak_info = {}
	reader = open(peakfile,'r')
	# Update sequence with each variant
	for line in reader:
		if(line.find("#")==0):
			continue
		tokens = line.strip().split("\t")
		if(tokens[3] in peak_info.keys()):
			print("Error!!! redundant id in peak file:")
			print(line)
			sys.exit(1)
		peak_info.update({tokens[3]:(tokens[0],tokens[2])})
	reader.close()
	return peak_info

if __name__ == "__main__":
	'''Write CDT file in new order.'''
	args = getParams()

	sua7_peaks = load_peaks(args.peaks)

	writer = open(args.output,'w')
	reader = open(args.input,'r')
	# Update sequence with each variant
	for line in reader:
		tokens = line.strip().split("\t")
		if(line.find("Chrom")==0):
			continue
		if(tokens[3] not in ["01_RP","02_STM","03_TFO","04_UNB"]):
			continue
		chr = "chrZ"
		coord = "-9999"
		id = tokens[6]
		score = "NaN"
		strand = tokens[1]

		if(tokens[34]!="." and tokens[34]!=""):
			refpt_info = tokens[34].split("_")
			chr = refpt_info[2]
			coord = refpt_info[3]
			score = "%s_%s" % (refpt_info[0], refpt_info[1])
		elif(id in sua7_peaks.keys()):
			chr,coord = sua7_peaks[id]
			score = "Sua7_CX_processed"
		else:
			chr = tokens[0]
			coord = int(tokens[14])
			if(strand=="-"):
				coord = str(coord+args.shift)
			else:
				coord = str(coord-args.shift)
			score = "imputed"
		writer.write("\t".join([chr,coord,coord,id,score,strand]) + "\n")
	reader.close()
	writer.close()

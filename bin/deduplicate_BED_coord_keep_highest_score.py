import sys, getopt
import re
import argparse

usage = """
Usage:
This script takes a peak BED and BAM file and tallies up an occupancy score of the tags in a strand-specific manner (to take advantage of XO profile).
-Tags mapping to sense strand are tallied in a +/-100bp window around the point -36bp (upstream) of the peak BED midpoint.
-Tags mappint to anti strand are tallied +/-100bp window around the point +36bp (downstream) of the peak BED midpoint.
-Sense Tags + Anti Tags = occupancy score replaced in peak BED file for output

Example: python deduplicate_BED_coord_keep_highest_score.py -i PEAK.bed -b NGS.bam -o PEAK_wscore.bed
"""

ALL_NT = ["A","T","C","G"]

def getParams():
	'''Parse parameters from the command line'''
	parser = argparse.ArgumentParser(description='This script takes a BED with a distance score and keeps the entry prioritizing upstream distances, then by closest distance.')
	parser.add_argument('-i','--input', metavar='bedfile', required=True, help='a BED file of coords to uniq (keep smallest absolute score)')
	parser.add_argument('-o','--output', metavar='outfile', required=True, help='a BED file of deduplicated entries')

	#parser.add_argument('-t','--title', metavar='figure_title', dest='title', required=True, help='')

	args = parser.parse_args()
	return(args)

def loadBED(bed_file):
	'''Load coordinate intervals bed_formatted_tokens'''
	all_peaks = []
	reader = open(bed_file,'r')
	# Update sequence with each variant
	for line in reader:
		tokens = line.strip().split("\t")
		#chr, midpoint, strand
		midpoint = int(tokens[1]) + ((int(tokens[2])-int(tokens[1])/2))
		dist_score = int(tokens[4])
		all_peaks.append(("%s_%i" % (tokens[0],midpoint),dist_score,line))
	reader.close()
	return all_peaks

if __name__ == "__main__":
	'''Write BED file in new order.'''
	args = getParams()

	# Load all BED file entries into memory
	COORDS = loadBED(args.input)
	idx_to_filter_out = []

	# Get filter indexes
	for i in range(len(COORDS)):
		if(i in idx_to_filter_out):
			continue
		if(COORDS[i][1] < -100 or COORDS[i][1] > 60):
			idx_to_filter_out.append(i)
			continue
		duplicates = [ j for j in range(len(COORDS)) if(COORDS[j][0] == COORDS[i][0]) ]
		if(len(duplicates) == 1):
			continue
		# Add all duplicates
		distances = [abs(COORDS[j][1]) for j in duplicates]
		del duplicates[distances.index(min(distances))]
		idx_to_filter_out.extend(duplicates)

	# Write to keep oritinal BED order
	writer = open(args.output,'w')
	for i in range(len(COORDS)):
		if(i in idx_to_filter_out):
			continue
		writer.write(COORDS[i][2])
	writer.close()

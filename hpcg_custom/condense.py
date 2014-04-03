#!/usr/bin/python

# This program condenses several lines of EDTs (1 line per DB) 
# into a single line of EDT (all DBs together)

import sys

size_array = {}
ic_array = {}
fp_array = {}
reads_array = {}
writes_array = {}

def ProcessFile(infilename, outfilename):
	with open(infilename, "r") as r, open(outfilename, "w") as w:
		for line in r:
			edt, count, total, size, ic, fp, reads, writes  = line.split(' ')[:8]

# If edt is not already there, add it
			if edt not in size_array:
				size_array[edt] = int(size, 16)
				ic_array[edt] = int(ic, 16)
				fp_array[edt] = int(fp, 16)
				reads_array[edt] = int(reads, 16)
				writes_array[edt] = int(writes, 16)
			else:
# If edt is already there, 
#    add to its existing counts
				size_array[edt] += int(size, 16)
				reads_array[edt] += int(reads, 16)
				writes_array[edt] += int(writes, 16)
# print once we've read all DBs
			if (count == total):
				w.write("%s %d %d %d %d %d\n" % (edt, size_array[edt], 
					ic_array[edt], fp_array[edt],
					reads_array[edt], writes_array[edt]))
				del size_array[edt]
				del ic_array[edt]
				del fp_array[edt]
				del reads_array[edt]
				del writes_array[edt]

ProcessFile(sys.argv[1], sys.argv[2])

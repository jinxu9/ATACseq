#!/usr/local/bin/python
import sys
import re
import pprint
from optparse import OptionParser
import os
import itertools

def main(argv):
    parser = OptionParser()
    usage = "usage: %prog [options] [inputs]"
    opts = OptionParser(usage=usage)
    opts.add_option("-o", default="qcTable.txt", help="OutputFile for qc")
    opts.add_option("-i", help="input configuration file")
    options, arguments = opts.parse_args()

    outTableFl = open(options.o, 'w')
    outTableFl.write("Sample" + "\t" + "TotalRawReads" + "\t" + "OverallAlignmentRate" + "\t" + "FinalMappedReads" + "\t" + "FinalMapped%" + "\t" + 
    "chrM%" + "\t" + "BlackListReads%" + "\t" + "MAPQFiltered%" + "\t" + "Duplicate%" + "\n")
    totalReads = 0.0
    chrMCount = 0
    configFl = open(options.i, 'r') 
    sampleLines = configFl.readlines()[4:]
    for line in sampleLines:
        line = line.rstrip()
	getQCForSample(line, outTableFl)

def getQCForSample(line, outFl):
    items = line.split("\t")
    logFile = items[2] + '/' + items[3] + '.sh.err'
    outFl.write(items[3] + "\t")
    #raw reads
    logFl = open(logFile, 'r')
    lines = logFl.readlines()
    oneLine = lines[19].rstrip()
    words = oneLine.split()
    outFl.write(words[0] + "\t")
    totalReads = int(words[0])

    #Overall Alignment Rate
    oneLine = lines[33].rstrip()
    words = oneLine.split()
    outFl.write(words[0] + "\t")
    logFl.close()

    logFile = items[2] + '/' + items[3] + '.sh.log'
    logFl = open(logFile, 'r')
    lines = logFl.readlines()

    #Final mapped reads
    oneLine = lines[25].rstrip()
    words = oneLine.split()
    mappedReads = int(words[0]) / 2
    outFl.write(str(mappedReads) + "\t")
    mappedPercent = mappedReads / float(totalReads) * 100
    outFl.write("%.2f" % mappedPercent + "\t")

    #chrM %
    oneLine = lines[23].rstrip()
    words = oneLine.split()
    chrMCount = int(words[0])
    chrMPercent = chrMCount / float(totalReads) / 2 * 100
    outFl.write("%.2f" % chrMPercent + "\t")

    #black list
    oneLine = lines[35].rstrip()
    words = oneLine.split()
    BLCount = int(words[0])
    BLPercent = BLCount / float(totalReads) / 2 * 100
    outFl.write("%.2f" % BLPercent + "\t")


    #QC filter
    oneLine = lines[27].rstrip()
    words = oneLine.split()
    beforeQC = int(words[0])

    oneLine = lines[21].rstrip()
    words = oneLine.split()
    afterQC = int(words[0])
    filterQC = beforeQC - afterQC - chrMCount
    qcPercent = filterQC / float(totalReads) /2 * 100
    outFl.write("%.2f" % qcPercent + "\t")

    #Duplicate %
    logFile = items[2] + '/QC/' + items[3] + '.Picard.log'
    logFl = open(logFile, 'r')
    lines = logFl.readlines()
    for line in lines:
    	oneLine = line.rstrip()
	words = oneLine.split()
	#if words[4] == 'Marking':
	keyword = 'Marking'
	if keyword in words:
	    duplicateCnt = int(words[5])
	    dupPercent = duplicateCnt / float(totalReads) / 2  * 100
	    outFl.write("%.2f" % dupPercent + "\n")


if __name__ == "__main__":
    main(sys.argv)

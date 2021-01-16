#!/usr/bin/python

#########################################
# Author: Yunjiang Qiu <yuq003@eng.ucsd.edu>
# File: bin2DI.py
# Create Date: 2015-03-22 16:50:15
#########################################

import sys
import argparse
import numpy as np
import gzip

def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("-i", "--input", dest="infile", required=True, help="input bin pair file")
    parser.add_argument("-g", "--gsize", dest="genome", required=True, help="genome size")
    parser.add_argument("-b", "--bsize", dest="bin_size", required=True, help="bin size")
    parser.add_argument("-c", "--chrom", dest="chrom", required=True, help="chrom")
    parser.add_argument("-w", "--wsize", dest="win_size", required=True, help="win size")
    parser.add_argument("-o", "--out", dest="out", required=True, help="output prefix")
    args = parser.parse_args()

    bin_size = args.bin_size.replace('Kb','000')
    bin_size = bin_size.replace('Mb','000000')
    win_size = args.win_size.replace('Kb','000')
    win_size = win_size.replace('Mb','000000')
    chrom = args.chrom.replace("chr", "").replace("X", "23")
    
    try:
        bin_size = int(bin_size)
        win_size = int(win_size)
    except ValueError:
        sys.exit("Unknown bin size %s or win size %s, please double check." % (args.bin_size, args.win_size))

    bins = win_size / bin_size

    with open(args.genome, 'r') as f:
        for line in f:
            key, value = line.rstrip().split('\t')
            if key == args.chrom:
                print("{0}".format(key))
                length = int(value)/bin_size + 1
                genome_matrix = np.zeros((length, length))

    with gzip.open(args.infile, 'r') as f:
        for line in f:
            line = line.rstrip().split('\t')
            val = float(line[2])
            if (np.isnan(val)):
                continue
            try:
                genome_matrix[int(line[0])/bin_size][int(line[1])/bin_size] = val
                genome_matrix[int(line[1])/bin_size][int(line[0])/bin_size] = val
            except (KeyError, ValueError) as e:
                pass

    with open(args.out, 'w') as f:
        for i in range(genome_matrix.shape[0]):
            if np.sum(genome_matrix[i]) > 0:
                starter = i
                break
            else:
                starter = i + 1
                #f.write('%s\t%d\t%d\t0\n' % (key, i * bin_size, (i+1) * bin_size))
            
        for i in range(starter, genome_matrix.shape[0]):

            A = 0
            B = 0
            
            for z in range(i-bins, i):
                if z >= starter and z < genome_matrix.shape[0]:
                    A += genome_matrix[i,z]

            for z in range(i+1, i+bins+1):
                if z >= starter and z < genome_matrix.shape[0]:
                    B += genome_matrix[i,z]

            E = (A + B) / 2
            if E == 0 or A == B:
                DI = 0
            else:
                DI = (B - A) /  abs(B - A) * ((A - E)**2/E + (B - E)**2/E)
            
            f.write('%s\t%d\t%d\t%.8g\n' % (chrom, i * bin_size, (i+1) * bin_size, DI))
            #f.write('%s\t%d\t%d\t%.8g\t%.8g\t%.8g\t%.8g\n' % (key,i * bin_size, (i+1) * bin_size, DI, A, B, E))

if __name__ == "__main__":
    sys.exit(main())


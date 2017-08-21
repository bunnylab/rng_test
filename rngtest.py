#!/usr/bin/env python
import numpy as np 
from scipy import stats
import argparse
import os

def monobits(sequence):
	"""return the proportion of 1/0, chisq statistic and p value for sequence"""
	freq = np.array([0,0])
	freq[0] = np.sum(sequence)
	freq[1] = sequence.size - freq[0]
	chisq, p = stats.chisquare(freq)
	print(freq)
		
	return freq[0]/(freq[0]+freq[1]), chisq, p
	
	
def mbits(sequence, blocks):
	"""returns the proportion of p<.05 for chisq test of all blocks of size m"""
	pvals = []
	for block in blocks:
		if(sequence.size//block > 0):
			pvals.append(np.empty(sequence.size//block))
			
	for sindex, size in enumerate(blocks):
		buffer = np.empty(size, dtype=np.uint8)
		freq = np.array([0,0])
		
		for index, bit in enumerate(sequence):
			buffer[index%size] = bit
			if ( (index+1)%size == 0 and index > 0):
				freq[0] = np.sum(buffer)
				freq[1] = buffer.size - freq[0]
				chisq, p = stats.chisquare(freq)
				pvals[sindex][((index+1)//size)-1] = p
							
	
	prop = []
	for index, block in enumerate(blocks):
		prop.append( (block, (pvals[index] <= 0.05).sum()/pvals[index].size) )
		
	return prop
		
		
	
parser = argparse.ArgumentParser(description='RNGTest: python reimplementation of nist rng stat tests')
parser.add_argument('keyfile', nargs=1, type=str, help='binary file to run tests on')
#parser.add_argument('-m', '--mono', nargs="?", type=int, default=80, help='target port')
#parser.add_argument('-s', '--sockets', nargs="?", type=int, default=150, help='# of sockets, default 150')
args = parser.parse_args()


def main(args):
	keyfile = open(args.keyfile[0],"rb")
	stat = os.stat(args.keyfile[0])
	bits = np.empty(stat.st_size*8, dtype=np.uint8)
	
	counter = 0
	byte = keyfile.read(1)
	while byte:
		for b in byte:
			for i in range(8):
				bits[counter] = (b >> i) & 1
				counter+=1
		byte = keyfile.read(1)
						
	# Run our tests
	# Monobits
	pr, chi, p = monobits(bits)
	print("Monobits:")
	print("proportion of 1s={:.3} chi={} p={:.3}\n".format(pr, chi, p))

	# MBits
	prop = mbits(bits, [8, 16, 32, 64])
	print("Mbits:")
	for pr in prop:
		print("{} bit bin - % pval<.05: {}".format(pr[0], pr[1]))

if __name__ == "__main__":
	main(args)
	
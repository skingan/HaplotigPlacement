#!/usr/bin/env python

# Sarah B. Kingan
# 20 September 2017
# 
# Pacific Biosciences
# Applications Lab
#
# Convert nucmer.coords files for all haplotigs aligned
# to primary contigs into ncbi placement file
#
###########################################################################################################

# import libraries
import numpy as np
import subprocess
import argparse

###########################################################################################################

desc ='''Make haplotig placement file from individual coords file of all haplotigs to primary'''
parser = argparse.ArgumentParser(description=desc)
parser.add_argument("infile", help="coords file")
args = parser.parse_args()



# print header
header=['hID','pID', 'ori', 'pstart', 'pstop', 'p_alnL', 'pL', 'hstart', 'hstop', 'h_alnL', 'hL', 'h_start_tail', 'h_stop_tail']
print '\t'.join(str(h) for h in header)

def file_len(fname):
	#http://stackoverflow.com/questions/845058/how-to-get-line-count-cheaply-in-python
	p = subprocess.Popen(['wc', '-l', fname], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	result, err = p.communicate()
	if p.returncode != 0:
		raise IOError(err)
	return int(result.strip().split()[0])


coords_file = args.infile

# htigs that align to primary
if file_len(coords_file) >= 5:

	d=np.loadtxt(open(coords_file, "rb"), skiprows=4, dtype="str", ndmin=2)

	haplotigs=list(set(list(d[:,10])))

	for h in haplotigs:
		output=['hID','pID', '+', 'pstart', 'pstop', 'p_alnL', 'pL', 'hstart', 'hstop', 'h_alnL', 'hL', '0', '0']
		subset = []
		for i in range(0,d.shape[0]):
			if d[i,10] == h:
				subset.append(i)
		output[0]=d[subset[0],10]	#htig name
		output[1]=d[subset[0],9]	#pcontig name
		output[6]=d[subset[0],7]	#pcontig length
		output[10]=d[subset[0],8]	#hcontig length
		h=np.array(d[subset[0],[2,3]],dtype=int) # coax into right data type to use amin
		p=np.array(d[subset[0],[0,1]],dtype=int)
		output[7]=np.amin(h) # hstart
		output[8]=np.amax(h) # hstop
		output[3]=np.amin(p) # pstart
		output[4]=np.amax(p) # pstop
		output[5]=int(output[4])-int(output[3])+1	# p aln L
		output[9]=int(output[8])-int(output[7])+1	# h aln L
# tail regions
		if int(output[7]) != 1: 		# if aln doesn't start at begining of htig
			output[11]=int(output[7])-1
		if int(output[8]) != int(output[10]): 	# if aln doesn't stop at end of htig
			output[12]=int(output[10])-int(output[8])
# alignment orientation
		if np.argmin(p) > np.argmax(p):
			output[2] = '-'
# htigs that don't align to primary
# print to stdout
		print '\t'.join(str(o) for o in output)


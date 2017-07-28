#!/usr/bin/env python

# Sarah B. Kingan
# 28 January 2017
# Modified by Ivan Sovic
# 24 July 2017
# Pacific Biosciences
# Applications Lab
#
# Convert directory of nucmer.coords files for haplotigs aligned
# to primary contiga into ncbi placement file
#
# Modified version for collaboration with dovetail
# To improve their HiRise scaffolding program
#
#
###########################################################################################################

# import libraries
import glob
import numpy as np
import subprocess

###########################################################################################################

# process each coords file
# coords files created with shell script:
#
#REF_DIR=/path/to/dir/containing/reference/file
#REF=$REF_DIR/referenceFilePrimaryContigsAndHaplotigs.fasta
#HAP_LIST=`cat $REF_DIR/listOfHaplotigIDs.txt`
#
#
#module load samtools
#module load mummer
#
#for h in $HAP_LIST
#do
#	samtools faidx $REF $h | sed 's/|quiver|arrow|arrow//g' > hap.fa  	# query file of haplotig seq
#	p=`echo $h | sed 's/_[0-9]\{3\}//'`					# associated primary contig id
#	samtools faidx $REF $p | sed 's/|quiver|arrow|arrow//g' > prim.fa	# reference file of primary contig seq
#	prefix=`echo $h | sed 's/|quiver|arrow|arrow//'`			# prefix for nucmer run
#	nucmer -maxmatch prim.fa hap.fa -prefix $prefix				# nucmer run
#	delta-filter -g $prefix.delta > $prefix.global.delta			# 1-1 global aln without rearragments
#	show-coords -qTl $prefix.global.delta > $prefix.coords			# coords format, tab delimited, sorted by query (htig), with ref and query lengths
#done;
#
# bulk run: batch_run_nucmer.sh 
# creates chunked lists of haplotigs,
# sends mummer4.0.0 jobs to scheduler with run_nucmer4.sh
#



# print header
# header=['#alt_asm_name','prim_asm_name','alt_scaf_name','parent_type','parent_name','ori','alt_scaf_start','alt_scaf_stop','alt_scaf_len','parent_start','parent_stop','parent_len','alt_start_tail','alt_stop_tail']
header=['#secondary_contig_ID','primary_contig_ID','ori','secondary_contig_start','secondary_contig_stop','secondary_contig_len','primary_contig_start','primary_contig_stop','primary_contig_len','secondary_contig_start_tail','secondary_contig_stop_tail']
print '\t'.join(str(h) for h in header)


def file_len(fname):
	#http://stackoverflow.com/questions/845058/how-to-get-line-count-cheaply-in-python
	p = subprocess.Popen(['wc', '-l', fname], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	result, err = p.communicate()
	if p.returncode != 0:
		raise IOError(err)
	return int(result.strip().split()[0])


fileList = glob.glob('*.coords')
for file in fileList:
# htigs that align to primary
	if file_len(file) >= 5:
		output=['000000F_000', '000000F', '+', 'hstart', 'hstop', 'hlen', 'pstart', 'pstop', 'plen', '0', '0']
		d=np.loadtxt(open(file, "rb"), skiprows=4, dtype="str", ndmin=2)
		output[0]=d[0,10]	#htig name
		output[1]=d[0,9]	#pcontig name
		h=np.array(d[:,[2,3]],dtype=int) # coax into right data type to use amin
		p=np.array(d[:,[0,1]],dtype=int)
		output[3]=np.amin(h) # hstart
		output[4]=np.amax(h) # hstop
		output[5]=int(d[0,8])     # hlen
		output[6]=np.amin(p) # pstart
		output[7]=np.amax(p)# pstop
		output[8]=int(d[0,7])    # plen
# tail regions
		if int(output[3]) != 1: 		# if aln doesn't start at begining of htig
			output[9]=int(output[3])-1
		if int(output[4]) != int(d[0,8]): 	# if aln doesn't stop at end of htig
			output[10]=int(d[0,8])-int(output[4])
# alignment orientation
		if np.argmin(p) > np.argmax(p):
			output[2] = '-'
# htigs that don't align to primary
	else:
		output=[file.split('.')[0],'na','na','na','na','na','na','na','na','na','na']
# print to stdout
	print '\t'.join(str(o) for o in output)


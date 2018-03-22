source /mnt/software/Modules/current/init/bash
module load samtools
REF=$1
DIR=$PWD/HaplotigPlacement
grep '_' $REF.fai | cut -f1 > htigs.txt
for h in `cat htigs.txt`
do
	hbase=$(echo $h | sed 's/|arrow//')			# base name of htig
	samtools faidx $REF $h | sed 's/|arrow//' > $hbase.fa     # query file of haplotig seq
	p=`echo $h | sed 's/_[0-9]\{3\}//'`                             # associated primary contig id
	pbase=$(echo $p | sed 's/|arrow//')			# base name of primary
	samtools faidx $REF $p | sed 's/|arrow//g' > $pbase.fa    # reference file of primary contig seq
	qsub $DIR/run_nucmer4.sh $pbase $hbase				# run nucmer and filtering
done;

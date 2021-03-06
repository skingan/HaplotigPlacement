#!/bin/bash
#$ -S /bin/bash
#$ -N ncmr
#$ -cwd
#$ -q def66
#$ -pe smp 2

module load mummer/4.0.0

REF=$1
QRY=$2
PREFIX=$QRY

nucmer --maxmatch -l 100 -c 500 -t 2 -p $PREFIX $REF.fa $QRY.fa
delta-filter -g $PREFIX.delta > $PREFIX.g.delta
show-coords -qTl -L 5000 $PREFIX.g.delta > $PREFIX.coords

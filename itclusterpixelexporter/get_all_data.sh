#!/bin/bash

cd ../

set -e

cmsenv

user=$USERNAME
initial=${user:0:1}
cernbox_dir=/eos/home-$initial/$user/

for file in `ls ./itdigiexporter/data/datasets/*.txt`
do
	f=${file%.txt}
	f1=${f##*/}

	output_file=output_$f1.root
	
	echo $f1 run starting

	cmsRun itdigiexporter/python/ITdigiExporter.py dataset=$f1 xrdredirector=file: output=$output_file

	mv ./$output_file $cernbox_dir

done


#!/bin/bash

module load r/4.0.2
module load boost


for m in {1..100}
do
	fn_out="seed""$m.out"
        jb_out="seed""$m.o"
        fn_err="seed""$m.e"
		subCommand="Rscript euler_scripts/cluster_euler_seed""$m.R"
        echo sbatch -n 1 --time=72:00:00 --mem-per-cpu=64000 -o $jb_out -e $fn_err --wrap='"'${subCommand}'"'
#        sbatch -n 1 --time=72:00:00 --mem-per-cpu=64000 -o $jb_out -e $fn_err --wrap='"'${subCommand}'"'
		eval "$subCommand"

done



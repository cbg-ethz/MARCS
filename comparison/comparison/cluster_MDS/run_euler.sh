#!/bin/bash

module load gcc/8.2.0 r/4.2.2

echo "module load gcc/8.2.0 r/4.2.2" > launcher.sh

for m in {1..100}
do
	fn_out="seed""$m.out"
        jb_out="seed""$m.o"
        fn_err="seed""$m.e"
		subCommand="Rscript --vanilla euler_scripts/cluster_euler_seed""$m.R"
        echo sbatch -n 1 --time=12:00:00 --mem-per-cpu=64000 -o $jb_out -e $fn_err --wrap='"'${subCommand}'"'  --nodefile=../nodelist.txt >> launcher.sh
        #sbatch -n 1 --time=96:00:00 --mem-per-cpu=64000 -o $jb_out -e $fn_err --wrap='"'${subCommand}'"' 
	#eval "$subCommand"

done



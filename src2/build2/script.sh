#!/bin/bash

base_dir=${PWD}
executable_path="/home/lorenzo/Op2 with Martino/Operational_Research_2/src2/build2/"

nnodes=(500)
models=("Cplex")
methods=("Benders")

for model in ${models[@]}
do
	for method in ${methods[@]}
	do
		for nnode in ${nnodes[@]}
		do
			for ((seed=0; seed<=9; seed++))
			do
				echo "${model}_${method}_${seed},"
				cd "${executable_path}" && ./tsp -n ${nnode} -tl 3600 -seed ${seed} -model ${model} -method ${method} -patch 0 -post 0 #>> ${base_dir}/${model}_${method}_${nnode}_${seed}.txt
			done
		done
	done
done
printf "END\n"





#!/bin/bash 
for i in `seq 1 94`
do
  cp -v tcf_DUMMY.pbs tcf_${i}.pbs
  sed -i "s#DUMMY#${i}#g" tcf_${i}.pbs
  cp -v tcfint_DUMMY.f ./tcfint_${i}.f
  sbatch tcf_${i}.pbs
done

exit

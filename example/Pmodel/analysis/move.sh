#!/bin/bash/

mkdir extras

mv *.log extras

mv *.pbs extras

mv *.exe extras

mv *.f extras

for i in {24..45}
	do
		j=$(($i-1))
		echo ${i} "->" ${j}
		mv m1CAsimint_${i} m1CAsimint_${j}
	
done

mv m1CAsimint_* ../tcfdata

exit

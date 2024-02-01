#!/bin/bash

FILE=treemix_input.gz
ENDi=10
echo "**** RUNNING TREEMIX ****"
	
for i in $(seq 1 $ENDi)
do
treemix -i $FILE -noRoot -m $i -o $FILE_${i} > treemix_${i}.log &
done


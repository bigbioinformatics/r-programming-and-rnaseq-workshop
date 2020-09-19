#!/bin/bash

file="srrs_to_analyze.txt"

while read line
do
  echo $line
  prefetch $line
  fasterq-dump --split-files -e 12 -p -O $line $line/$line.sra
  
done < $file
wait


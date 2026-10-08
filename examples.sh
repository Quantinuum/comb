#!/bin/bash
false_positives=""
	true_negatives=""
	for f in ./examples/fault-tolerant/*.run.gz
	do
  gunzip -k $f -c > /tmp/test-run
	LD_LIBRARY_PATH=./m4ri/lib ./comb /tmp/test-run 1
	if [ $? -eq 0 ] 
	then
	echo "Success"
	else
	true_negatives="$true_negatives $f"
	fi
	done
	for f in ./examples/non-fault-tolerant/*.run.gz
	do
  gunzip -k $f -c > /tmp/test-run
	LD_LIBRARY_PATH=./m4ri/lib ./comb /tmp/test-run 1
	if [ $? -eq 1 ] 
	then
	echo "Success"
	else
	true_negatives="$false_positives $f"
	fi
	done
	rm -f *.run
  echo True negatives: $true_negatives
  echo False positives: $false_positives

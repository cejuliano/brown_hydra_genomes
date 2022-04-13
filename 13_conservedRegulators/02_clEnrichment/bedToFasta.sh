#!/bin/bash

cd mgBeds

for arg in *bed
do
	newName="${arg/bed/fa}"
	
	bedtools getfasta -nameOnly -s -fi ../../../clytiaG.fa -bed $arg > "$newName"
	
	gsed -i 's/(.*)//g' "$newName"
done


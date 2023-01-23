#!/bin/bash

cd "$1"

while read -ra array; do
  ar1+=("${array[0]}")
  ar2+=("${array[1]}")
done < "$1".genome

for i in "${ar1[@]}"; do
	for j in "${ar1[@]}"; do
		if [[ "$i" != "$j" ]]; then

			java -Xms512m -Xmx2048m -jar ../juicer_tools.jar \
				dump observed KR inter.*hic "$i" "$j" \
				BP 100000 > "$i"_"$j".txt
		fi
	done
done
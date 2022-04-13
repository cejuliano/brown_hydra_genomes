#! /bin/bash -l

shopt -s nullglob

# The only argument needed for this script is the unique sample ID
prefix="$1"

conda activate deepEnv

#index bam files
echo "indexing bam files"
for arg in "$prefix"[^_]*_final.bam;
do
	echo "$arg"
	samtools index "$arg"
done

echo "done"

#namesort bam files
echo "shifting reads"

for arg in "$prefix"[^_]*_final.bam;
do
	rep="${arg/_final.bam/}"
	echo "$rep"

	alignmentSieve -b "$arg" -o "$rep"_final_shift.bam -p 16 --ATACshift

	echo "resorting shifted reads"

	samtools sort -T "$rep".sort -o "$rep"_final_shift.sort.bam "$rep"_final_shift.bam

	mv "$rep"_final_shift.sort.bam "$rep"_final_shift.bam

	samtools index "$rep"_final_shift.bam
done

echo "done"

#make self pseudoreplicates
echo "making self-pseudoreplicates"

for arg in "$prefix"[^_]*_final_shift.bam;
do
	rep="${arg/_final_shift.bam/}"
	echo "$rep"

	echo "generating first psuedoreplicate"
	samtools view -s 1234.5 -b -@ 16 -o "$rep"_PR1_final_shift.bam "$arg"

	echo "generating second psuedoreplicate"
	samtools view "$rep"_PR1_final_shift.bam | cut -f 1 > "$rep"_PR1_qname.txt

	java -Xmx32g -jar resources/picard.jar FilterSamReads \
		I="$arg" \
		O="$rep"_PR2_final_shift.bam \
		READ_LIST_FILE="$rep"_PR1_qname.txt \
		VALIDATION_STRINGENCY=SILENT \
		FILTER=excludeReadList \
		QUIET=true

	rm "$rep"_PR1_qname.txt

done

#pool replicates

echo "pooling replicates"

samtools merge "$prefix"_MG_final_shift.bam "$prefix"[1-9]_final_shift.bam
samtools index "$prefix"_MG_final_shift.bam


#split pooled rep into psuedoreps (same number as total reps)

echo "generating pseudoreplicates from pooled counts"

numReps=("$prefix"[^_]*_final.bam)
numReps=${#numReps[@]}

echo "splitting into $numReps files"

count=1

cp "$prefix"_MG_final_shift.bam "$prefix"_MG_final_shift.sub.bam

while [ $numReps -gt 1 ]
do

	subSampleValue="$(Rscript resources/generateSubsampleValue.R $numReps | cut -f 2)"
	echo "$numReps"
	echo "$count"
	echo "$subSampleValue"
	samtools view -s "$subSampleValue" -b -@ 16 \
		-o "$prefix"_MG_PR"$count"_final_shift.bam "$prefix"_MG_final_shift.sub.bam

	samtools view "$prefix"_MG_PR"$count"_final_shift.bam | cut -f 1 > "$prefix"_PR_qname.txt


	java -Xmx32g -jar resources/picard.jar FilterSamReads \
		I= "$prefix"_MG_final_shift.sub.bam \
		O= "$prefix"_MG_final_shift.sub.tmp.bam \
		READ_LIST_FILE="$prefix"_PR_qname.txt \
		VALIDATION_STRINGENCY=SILENT \
		FILTER=excludeReadList \
		QUIET=true

	rm "$prefix"_MG_final_shift.sub.bam

	mv "$prefix"_MG_final_shift.sub.tmp.bam "$prefix"_MG_final_shift.sub.bam

	count=$(( $count + 1 ))
	numReps=$(( $numReps - 1 ))

done

mv "$prefix"_MG_final_shift.sub.bam "$prefix"_MG_PR"$count"_final_shift.bam

rm "$prefix"_PR_qname.txt

echo "Calling peaks"

for arg in "$prefix"*_final_shift.bam;
do
	rep="${arg/_final_shift.bam/}"
	echo "$rep"
	macs2 callpeak \
		-t "$arg" -f BAMPE -n "$rep" -g 9e8 -p 0.1 \
		--nomodel --keep-dup all

	sort -k 8gr,8gr "$rep"_peaks.narrowPeak \
		| awk 'BEGIN{OFS="\t"}{$4="Peak_"NR ; print $0}' > "$rep".narrowPeak

	rm -f "$rep"_peaks.narrowPeak "$rep"_peaks.xls

done

rm *summits*

echo "done"

#perform idr for all reps
echo "Performing IDR on biological replicates"

for i in "$prefix"[1-9].narrowPeak
do
	for j in "$prefix"[1-9].narrowPeak
	do
		if [[ "$i" < "$j" ]]; then
			echo "$i"
			echo "$j"

			inputFile1="$i"
			prefix1="${inputFile1/.narrowPeak/}"

			inputFile2="$j"
			prefix2="${inputFile2/.narrowPeak/}"

			idr --samples "$i" "$j" \
				--peak-list "$prefix"_MG.narrowPeak --input-file-type narrowPeak \
				--output-file "$prefix1"_"$prefix2".idr --rank p.value --soft-idr-threshold 0.1

			IDR_THRESH_TRANSFORMED=$(awk -v p=0.1 'BEGIN{print -log(p)/log(10)}')

			awk 'BEGIN{OFS="\t"} $12>='"${IDR_THRESH_TRANSFORMED}"' \
				{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}' "$prefix1"_"$prefix2".idr | \
					sort | \
					uniq | \
					sort -k7n,7n > "$prefix1"_"$prefix2".IDR.narrowPeak
		fi
	done
done

echo "done"

#perform idr for all self-pseudoreps
echo "Perfoming IDR on self-pseudoreplicates"

for arg in "$prefix"[1-9]_final_shift.bam
do
	rep="${arg/_final_shift.bam/}"
	echo "$rep"

	inputFile1="$rep"_PR1.narrowPeak
	prefix1="${inputFile1/.narrowPeak/}"

	inputFile2="$rep"_PR2.narrowPeak
	prefix2="${inputFile2/.narrowPeak/}"

	idr --samples "$inputFile1" "$inputFile2" \
		--peak-list "$rep".narrowPeak --input-file-type narrowPeak \
		--output-file "$prefix1"_"$prefix2".idr --rank p.value --soft-idr-threshold 0.1

	IDR_THRESH_TRANSFORMED=$(awk -v p=0.1 'BEGIN{print -log(p)/log(10)}')

	awk 'BEGIN{OFS="\t"} $12>='"${IDR_THRESH_TRANSFORMED}"' \
		{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}' "$prefix1"_"$prefix2".idr | \
			sort | \
			uniq | \
			sort -k7n,7n > "$prefix1"_"$prefix2".IDR.narrowPeak
done

echo "done"

#perform idr for all pooled pseudoreps
echo "Perfoming IDR on pooled pseudoreplicates"

for i in "$prefix"_MG_PR*.narrowPeak
do
	for j in "$prefix"_MG_PR*.narrowPeak
	do
		if [[ "$i" < "$j" ]]; then
			echo "$i"
			echo "$j"

			inputFile1="$i"
			prefix1="${inputFile1/.narrowPeak/}"

			inputFile2="$j"
			prefix2="${inputFile2/.narrowPeak/}"

			idr --samples "$i" "$j" \
				--peak-list "$prefix"_MG.narrowPeak --input-file-type narrowPeak \
				--output-file "$prefix1"_"$prefix2".idr --rank p.value --soft-idr-threshold 0.1

			IDR_THRESH_TRANSFORMED=$(awk -v p=0.1 'BEGIN{print -log(p)/log(10)}')

			awk 'BEGIN{OFS="\t"} $12>='"${IDR_THRESH_TRANSFORMED}"' \
				{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}' "$prefix1"_"$prefix2".idr | \
					sort | \
					uniq | \
					sort -k7n,7n > "$prefix1"_"$prefix2".IDR.narrowPeak
		fi
	done
done

echo "done"

conda activate deepEnv

#generate bigwig files
echo "Generating bigwig tracks"

for arg in "$prefix"[1-9]_final_shift.bam
do
	echo "$arg"
	bamCoverage -b "$arg" -o "${arg/.bam/.bw}" \
		-of "bigwig" -bs 10 -p 16 --normalizeUsing "CPM"
done

echo "$prefix"_MG_final_shift.bam
bamCoverage -b "$prefix"_MG_final_shift.bam -o "$prefix"_MG_final_shift.bw \
	-of "bigwig" -bs 10 -p 16 --normalizeUsing "CPM"

echo "done"

echo "cleaning up"

rm regions_"$prefix"*ATAC.bed Values_"$prefix"*.txt "$prefix"*matrix.gz "$prefix"*.idr

echo "done"



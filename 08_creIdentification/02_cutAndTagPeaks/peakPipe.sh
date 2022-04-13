#! /bin/bash -l

shopt -s nullglob

# The only argument needed for this script is the unique sample ID
prefix="$1"

#make self pseudoreplicates
echo "making self-pseudoreplicates"

for arg in "$prefix"_[0-9]*.final.bam;
do
	rep="${arg/.final.bam/}"
	echo "$rep"

	echo "generating first psuedoreplicate"
	samtools view -s 1234.5 -b -@ 16 -o "$rep"_PR1.final.bam "$arg"

	echo "generating second psuedoreplicate"
	samtools view "$rep"_PR1.final.bam | cut -f 1 > "$rep"_PR1_qname.txt

	java -Xmx60g -jar resources/picard.jar FilterSamReads \
		I="$arg" \
		O="$rep"_PR2.final.bam \
		READ_LIST_FILE="$rep"_PR1_qname.txt \
		VALIDATION_STRINGENCY=SILENT \
		FILTER=excludeReadList \
		QUIET=true

	rm "$rep"_PR1_qname.txt

done

#pool replicates

echo "pooling replicates"

samtools merge -f "$prefix"_MG.final.bam "$prefix"_[1-9].final.bam
samtools index "$prefix"_MG.final.bam


#split pooled rep into psuedoreps (same number as total reps)

echo "generating pseudoreplicates from pooled counts"

numReps=("$prefix"_[0-9].final.bam)
numReps=${#numReps[@]}

echo "splitting into $numReps files"

count=1

cp "$prefix"_MG.final.bam "$prefix"_MG_final.sub.bam

while [ $numReps -gt 1 ]
do

	subSampleValue="$(Rscript resources/generateSubsampleValue.R $numReps | cut -f 2)"
	echo "$numReps"
	echo "$count"
	echo "$subSampleValue"
	samtools view -s "$subSampleValue" -b -@ 16 \
		-o "$prefix"_MG_PR"$count".final.bam "$prefix"_MG_final.sub.bam

	samtools view "$prefix"_MG_PR"$count".final.bam | cut -f 1 > "$prefix"_PR_qname.txt


	java -Xmx60g -jar resources/picard.jar FilterSamReads \
		I= "$prefix"_MG_final.sub.bam \
		O= "$prefix"_MG_final.sub.tmp.bam \
		READ_LIST_FILE="$prefix"_PR_qname.txt \
		VALIDATION_STRINGENCY=SILENT \
		FILTER=excludeReadList \
		QUIET=true

	rm "$prefix"_MG_final.sub.bam

	mv "$prefix"_MG_final.sub.tmp.bam "$prefix"_MG_final.sub.bam

	count=$(( $count + 1 ))
	numReps=$(( $numReps - 1 ))

done

mv "$prefix"_MG_final.sub.bam "$prefix"_MG_PR"$count".final.bam

rm "$prefix"_PR_qname.txt

echo "generating bigwigs"

conda activate deepEnv

for arg in "$prefix"_[0-9].final.bam
do
	rep="${arg/.final.bam/}"
	echo "$rep"

	bamCoverage -b "$arg" \
		-o "$rep".bw \
		-of "bigwig" \
		-bs 10 \
		-p 24 \
		--normalizeUsing CPM \
		--exactScaling
done

bamCoverage -b "$prefix"_MG.final.bam \
	-o "$prefix"_MG.bw \
	-of "bigwig" \
	-bs 10 \
	-p 24 \
	--normalizeUsing CPM \
	--exactScaling \

conda deactivate

echo "Calling peaks"

for arg in 1 2 3 "MG" 
do
		rep="$prefix"_"$arg"
		echo "$rep"

		samtools sort -T $rep.sort -n -o $rep.ns.bam "$rep".final.bam

		bedtools bamtobed -bedpe -i $rep.ns.bam > $rep.bed

		cut -f 1,2,6 $rep.bed | \
			sort -k1,1 -k2,2n -k3,3n > $rep.fragments.bed

		bedtools genomecov -bg -i $rep.fragments.bed -g resources/aep.genome > $rep.bg

		rm $rep.bed $rep.fragments.bed $rep.ns.bam

		bash resources/SEACR/SEACR_1.3.sh $rep.bg IGG_"$arg".bg norm relaxed $rep

		bash resources/SEACR/SEACR_1.3.sh $rep.bg IGG_"$arg".bg norm stringent $rep

		Rscript resources/refBed.R "$rep".relaxed.bed

		Rscript resources/refBed.R "$rep".stringent.bed

		rm "$rep".relaxed.bed "$rep".stringent.bed

		mv "$rep".relaxed.rfmt.bed "$rep".relaxed.bed

		mv "$rep".stringent.rfmt.bed "$rep".stringent.bed

		for prR in "$rep"_PR*.final.bam
		do
			prep="${prR/.final.bam/}"
			echo "$prep"

			samtools sort -T $prep.sort -n -o $prep.ns.bam "$prep".final.bam

			bedtools bamtobed -bedpe -i $prep.ns.bam > $prep.bed

			cut -f 1,2,6 $prep.bed | \
				sort -k1,1 -k2,2n -k3,3n > $prep.fragments.bed

			bedtools genomecov -bg -i $prep.fragments.bed -g resources/aep.genome > $prep.bg

			rm $prep.bed $prep.fragments.bed $prep.ns.bam

			bash resources/SEACR/SEACR_1.3.sh $prep.bg IGG_"$arg".bg norm relaxed $prep

			bash resources/SEACR/SEACR_1.3.sh $prep.bg IGG_"$arg".bg norm stringent $prep

			Rscript resources/refBed.R "$prep".relaxed.bed

			Rscript resources/refBed.R "$prep".stringent.bed

			rm "$prep".relaxed.bed "$prep".stringent.bed

			mv "$prep".relaxed.rfmt.bed "$prep".relaxed.bed

			mv "$prep".stringent.rfmt.bed "$prep".stringent.bed
		done

	done

echo "done"

conda activate deepEnv

#perform idr for all reps

echo "Performing IDR on biological replicates"

for i in "$prefix"_[0-9].relaxed.bed
do
	for j in "$prefix"_[0-9].relaxed.bed
	do
		if [[ "$i" < "$j" ]]; then
			echo "$i"
			echo "$j"

			inputFile1="$i"
			prefix1="${inputFile1/.relaxed.bed/}"

			inputFile2="$j"
			prefix2="${inputFile2/.relaxed.bed/}"

			idr --samples "$i" "$j" \
				--peak-list "$prefix"_MG.relaxed.bed --input-file-type bed \
				--output-file "$prefix1"_"$prefix2".idr \
				-i 0.1 --rank 5

			bedtools intersect -e -f 0.25 -F 0.25 -c \
				-a "$prefix"_MG.relaxed.bed \
				-b "$prefix1"_"$prefix2".idr > "$prefix1"_"$prefix2".int.bed

		fi
	done
done

echo "done"

#perform idr for all reps
echo "Performing IDR for all self-pseudoreps"

for arg in "$prefix"_[0-9].relaxed.bed
do
	rep="${arg/.relaxed.bed/}"
	echo "$rep"

	inputFile1="$rep"_PR1.relaxed.bed
	prefix1="${inputFile1/.relaxed.bed/}"

	inputFile2="$rep"_PR2.relaxed.bed
	prefix2="${inputFile2/.relaxed.bed/}"

	idr --samples "$inputFile1" "$inputFile2" \
		--peak-list "$rep".relaxed.bed --input-file-type bed \
		--output-file "$prefix1"_"$prefix2".idr \
		-i 0.1 --rank 5

done

echo "done"

#perform idr for all pooled pseudoreps
echo "Perfoming IDR on pooled pseudoreplicates"

for i in "$prefix"_MG_PR[0-9].relaxed.bed
do
	for j in "$prefix"_MG_PR[0-9].relaxed.bed
	do
		if [[ "$i" < "$j" ]]; then
			echo "$i"
			echo "$j"

			inputFile1="$i"
			prefix1="${inputFile1/.relaxed.bed/}"

			inputFile2="$j"
			prefix2="${inputFile2/.relaxed.bed/}"

			idr --samples "$inputFile1" "$inputFile2" \
				--peak-list "$prefix"_MG.relaxed.bed --input-file-type bed \
				--output-file "$prefix1"_"$prefix2".idr \
				-i 0.1 --rank 5

			bedtools intersect -e -f 0.25 -F 0.25 -c \
				-a "$prefix"_MG.relaxed.bed \
				-b "$prefix1"_"$prefix2".idr > "$prefix1"_"$prefix2".int.bed

		fi
	done
done

echo "done"

conda deactivate

#generating consensus peaklists

Rscript resources/getCon.R "$prefix" "$prefix"*_[0-9].int.bed

echo "# of peaks in biologically reproducible peak set:"

wc -l consensus"$prefix".bed

Rscript resources/getCon.R "$prefix"_PR "$prefix"*_PR[0-9].int.bed

echo "# of peaks in ideal reproducible peak set:"

wc -l consensus"$prefix"_PR.bed

echo "# of peaks across self-pseudoreplicates"

wc -l "$prefix"_[0-9]_PR[0-9]_"$prefix"_[0-9]_PR[0-9].idr

echo "done"

echo "calculating FRiP"

conda activate deepEnv

plotEnrichment -b "$prefix"_[1-9].final.bam \
	--BED "$prefix"_MG.relaxed.bed consensus"$prefix".bed \
	--regionLabels "Permissive Peaks" "High Confidence Peaks" \
	-o "$prefix"_FRiP.pdf -p 24

conda deactivate

echo "done"

echo "plotting size distrubution"

for arg in "$prefix"_[0-9].final.bam
do
	rep="${arg/.final.bam/}"
	echo "$rep"

	java -Xmx60g -jar resources/picard.jar CollectInsertSizeMetrics \
		I="$arg" \
		O="$rep"_insert_size_metrics.txt \
		H="$rep"_insert_size_histogram.pdf

done

echo "done"


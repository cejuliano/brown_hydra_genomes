#!/bin/bash -l

source ~/.bash_profile

#The only input for the pipeline script is the name of a file 
#within a subdirectory called query.fa.split/
fileN="$1"

fileN="${fileN/query.fa.split\//}"

#here we define a run name, derived from the input file name
runN="${fileN/.fa/}"

#we'll make an output folder named after the run name to put our output
mkdir "$runN"

cd "$runN"

echo "$runN"

echo "$fileN"

#many of the fasta sequences in the input sequences have long headers with characters
#that can cause parsing errors. These long headers have useful information though,
#so we just pull out the full headers to set aside, then truncate the headers in our actual
#input files
grep ">" ../query.fa.split/"$fileN" > headers.txt

sed 's/|.*//g' ../query.fa.split/"$fileN" > inSeqs.r.fa

#using cd-hit to get rid of any possibly redundant sequences in the input
cd-hit -i inSeqs.r.fa -o inSeqs.fa -c 0.95 -G 0 -aL 0.8

#now we split the multi-fasta into individual files
#which we'll then iterate through one by one
seqkit split --quiet -i -f -O subseqs inSeqs.fa

#initialize empty output files
echo -n > fullRes.gff3

echo -n > rawExo.txt

echo -n > reformatted.gff

#iterate through each sequence from the input file
for arg in subseqs/*fa
do
	#extract the sequence ID from it's file name
	name=$(echo "$arg" | sed "s/subseqs\/inSeqs.id_//g;s/.fa//g")
	echo $name
	
  #Exonerate's cdna2genome algorithm can produce high quality alignments, but
  #it's also prohibitively slow when it has to deal with a genome-sized search
  #space. We use basically the same solution as the MAKER pipeline, and use an
  #initial BLAST search to produce a rough alignment that we can use to
  #identify the general coordinates for the target gene. Then we can generate a
  #high quality Exonerate alignment based on a much smaller search space.
	blastn -query $arg -db ../AEPgenome -outfmt '17 SR' -max_target_seqs 1 > out.sam
	
	#if BLAST produced no alignments, move on to the next query sequence
	if [ ! -s out.sam ]
	then
		continue
	fi
	
	#take the SAM alignments from blast and output the genomic sequence that spans
	#the entirety of the alignment
	samtools view -b out.sam | \
		bedtools bamtobed -tag AS -i - | \
		sort -k1,1 -k2,2n | \
		#merge nearby alignments into a single chunk
		#have to be within 20Kb (max intron size)
		bedtools merge -c 5 -o sum -d 20000 -s -i - | 
		sort -k 4 | \
		tail -n 1 | \
		#add 20 Kb on either side of chunk for good measure
		bedtools slop -i - -g ../aep.genome -b 20000 |
		bedtools getfasta -fi ../aep.genome.fullsoft.fa -bed - \
		> searchSpace.fa
	
  #pull all possible ORF coordinates from transcript sequence
	getorf $arg orfList.txt -find 1
	
	#pull out coordinates for the longest ORF
	seqkit sort --quiet -l -r orfList.txt | \
		head -n 1 | \
		sed 's/.*\[//g;s/\]//g;s/- //g' > coords.txt
	
  #reformat coordinates to work with Exonerate (need to include stop codon in total length)
	/usr/bin/Rscript ../addStopCoord.R coords.txt
		
	rm orfList.txt
	
  #final ORF coords reformatting
	echo $name "+" | \
		cat - coords.txt | \
		tr '\n' ' ' > annot.txt
		
	rm coords.txt
	
	#run the actual exonerate algorithm
	exonerate -q $arg \
		-t searchSpace.fa \
		-E TRUE \
		-m cdna2genome \
		--percent 25 \
		--showalignment FALSE \
		--showvulgar FALSE \
		--showtargetgff TRUE \
		-n 1 \
		-S FALSE \
		--annotation annot.txt \
		--softmasktarget TRUE \
		--seedrepeat 4 \
		--geneseed 250 | \
		sed 's/utr3b/utr3/g;/^#[^#]/d;/^#$/d;/similarity/d' > exo.txt
	#exonerate outputs a lot of additional text beyond just the GFF that we need to get rid of
	csplit -f hit exo.txt '/gff-version/' '{*}'
		
	cat exo.txt >> rawExo.txt
	
	rm hit00
	
	#if exonerate produced no alignments then move on to the next query sequence
	if [ ! -f hit01 ]
	then
		continue
	fi
	
	#this script fixes the all the formatting problems with the exonerate GFF output
	/usr/bin/Rscript ../reformatGff.R hit01
	
	sed -i '/^##/d' hit01.gff
	
	cat hit01.gff >> reformatted.gff

	echo "##gff-version 2" | cat - hit01.gff > geneRes.gff
	
	#delete UTR rows, we'll add them back in later using AGAT 
	#(does a better job with formatting/accuracy)
	sed '/utr5/d;/utr3/d' geneRes.gff > geneRes.utr.gff
	
	conda activate agatEnv

	agat_convert_sp_gxf2gxf.pl -g geneRes.utr.gff \
		-c "gene_id" \
		-gvi 2 \
		-gvo 3 \
		-o geneRes.gff3
	
	conda deactivate
	
  #final tweaks to add proper parent ID to gene model
	/usr/bin/Rscript ../fixParents.R geneRes.gff3
		
	rm hit0*
	
	#add result to full list of gene models
	cat geneRes.pfix.gff3 >> fullRes.gff3
			
done

#################################################################################################
#                   Mapping total RNA-seq data & lncRNA detection                               #                    
#################################################################################################



######    TOOLS   #####
# TopHat v2.1.1
# bedtools v2.27.0
# samtools v1.9
# R v3.4.4



#####   Mapping total RNA-seq data with TopHat    #####

# N. castellii genome downloaded from https://fungi.ensembl.org

# bowtie 2 genome index 

genome_index=(...) # path to genome index

indexes=(...) # fastq indexes

# mapping paired-end with TopHat, filtering multi-hits, sorting and indexing with samtools 

for index in ${indexes[*]};
do 
 
	R1="$index.R1.fastq.gz" # read 1
	R2="$index.R2.fastq.gz" # read 2
 
	# tophat mapping
	tophat2 -g 2 --library-type fr-firststrand -I 2000 --read-edit-dist 3 -r 155 --mate-std-dev 80 --no-coverage-search -p 6 --keep-fasta-order -N 3 -o "$index" "$genome_index" "$R1" "$R2"

	# filter multi-hits and sort
	samtools view -h "$index/accepted_hits.bam" | grep -E '^@|NH:i:1' | samtools view -bh - | samtools sort -o "$index/$index.unique.sorted.bam" 
  
	# index with samtools
	samtools index "$index/$index.unique.sorted.bam"
done



######    Annotation of lncRNA    ##### 

## list of bams
bams=
for index in ${indexes[*]}; 
do 
  bams=$(echo "$b" "$inded/$index.unique.sorted.bam"); done
done

## merge all bams (and remove ERCC, only for data from Szachnowski et al.), sort and index 
## only WT and rrp6 strains were used for Alcid & Tsukiyama data

samtools merge - $(echo ${bams[*]}) | samtools view -h - | awk 'OFS="\t" $3 ~ /Chr/ && $7 ~ /=/ || $1 ~ /^@/' | grep -Ev '^@SQ.*ERCC.*' | samtools view -bh - | samtools sort -o "$mergedBams"  
samtools index $mergedBams

## signal segmentation

# zinar_log.R : 
#   input : 
#     - bam file
#     - L : window size (min, max, step)
#     - s : threshold (min, max, step)
#   output :
#     for each window size and threshold, a gff file is produced with coordinates of segments for which 
#     the mean of signal (log2) over the window is above threshold

Rscript zinar_log.R -c 6 -L 5 -L 200 -L 5 -s 1.44 -s 432 -s 1.44 -v 1 -t 200000 -l "inverse" -o "$zinarDir" "$mergedBams"

## change scientific to decimal notation in gff files 
ls "$zinarDir" | grep "gff" | while read file; do
   Rscript sci_to_dec.R "$zinarDir/$file"
   rm "$zinarDir/$file.old";
done 

##  script R sci_to_dec.R
##  -----------------------------
##	file = commandArgs(trailingOnly=TRUE)[1]
##	gff <- read.csv(file, sep = "\t", header = F)
##	write.table(gff, file = paste0(file,".old"), quote = FALSE, sep = "\t", row.names = F, col.names = F)
##	gff[,4] <- format(gff[,4], scientific = FALSE)
##	gff[,5] <- format(gff[,5], scientific = FALSE)
##	write.table(gff, file = file, quote = FALSE, sep = "\t", row.names = F, col.names = F)
##  -----------------------------

## prepare data for lncRNA detection : 

## get counts for normalization :

## for data from data Szachnowski et al., count done on ERCC spike-in
cat "ERCC.gff" | awk 'OFS="\t" {print $1,$1,$4,$5,$6}' > annot_for_norm.saf

## for data from Alcid & Tsukiyama, count done on coding sequences
# cat "Naumovozyma_castellii_cbs_4309.ASM23734v1.38.gff3" | grep -Ev "^#" | awk -F'[\t;]' 'OFS = "\t" {if($3=="exon") print $10,$1,$4,$5,$7}' | grep "CCC" | sed 's/Name=//g' > annot_for_norm.saf

featureCounts -a annot_for_norm.saf -F "SAF" -B -C -p -s 2 -T 8 -O -o "counts_for_norm.tab" $(echo ${bams[*]}) 

## for all strains, get fragment number for FPKM calculation
for bam in ${bams[*]}; 
do
  index=$(basename $bam)
  samtools view "$bam" | wc -l | awk '{print $1/2}' > "nRead_${index}.txt"
done

## get N.castellii annotation with mRNA, sn(o)RNA, tRNA and rRNA for filtering transcripts aready annotated
grep -v "#" "Naumovozyma_castellii_cbs_4309.ASM23734v1.38.gff3" | awk '{if (($3 == "mRNA") || ($3 == "rRNA") || ($3 == "snoRNA") || ($3 == "snRNA") || ($3  == "tRNA")) print $0}' > "known_transcripts.gff"

## function to count reads and execute differential expression analysis 

detect_lncRNA () {

	## get gff file
	segmentation=$1

	## segmentation name
	name=$(basename $segmentation | sed -r 's/.gff//g')

	outSeg="DE_new_transcripts/$name"
	mkdir "$outSeg"

	## select new transcripts : size > 200 nt and no overlap with mRNA, sn(o)RNA, rRNA or tRNA
	# in : segmentation
	#      transcripts known
	# out : new_[segmentation name].gff
	
	known_transcripts=known_transcripts.gff
	awk 'OFS="\t" {if (($5-$4+1) >= 200) print $0}' "$segmentation" | bedtools intersect -nonamecheck -wao -s -a stdin -b "$known_transcripts" | awk 'OFS="\t" {if ($10 == ".") print $9, $1, $4, $5, $7}' | sed -r 's/ID=//g' > "$outSeg/new_$name.saf"
	
	## counting reads on new transcripts using featureCounts
	# in : annotation new transcripts selected and bam files
	# out : count table
  
  	## for data from Szachnowski et al., count were done on WT, xrn1-delta and dicer-delta strains (XUT and DUT detection) 
  	featureCounts -a "$outSeg/new_$name.saf" -F "SAF" -B -C -p -s 2 -T 4 -O -o "$outSeg/counts.tab" $WT_1 $WT_2 $xrn1_1 $xrn1_2 $dicer_1 $dicer_2
  
  	## for data from Alcid & Tsukiyama, count were done on WT and rrp6-delta strains (CUT detection)
  	# featureCounts -a "$outSeg/new_$name.saf" -F "SAF" -B -C -p -s 2 -T 4 -O -o "$outSeg/counts.tab" $WT_1 $WT_2 $rrp6_1 $rrp6_2
  
  	## Differential expression analysis, lncRNA detection
  
  	# in : count table
  	#      count for normalisation
  	#      nb reads tot (for each replicate and each strain) 
  	#      out directory
  
  	# out : table results from differential analysis (DESeq), 
  	#       gff file of lncRNA detected
	
  	## for data from Szachnowski et al., XUT, DUT and SUT detection
  
  	$nReadWT_1=... 
  	$nReadWT_2=... 
  	$nReadxrn1_1=... 
  	$nReadxrn1_2=... 
  	$nReaddicer_1=... 
  	$nReaddicer_2=...
  
  	Rscript DESeq_XUT_DUT.R "$outSeg/counts.tab" "counts_for_norm.tab" $nReadWT_1 $nReadWT_2 $nReadxrn1_1 $nReadxrn1_2 $nReaddicer_1 $nReaddicer_2 "$outSeg/"
  
  	## for data from Szacnowski et al., CUT and SUT detection
  
  	# $nReadWT_1=...
  	# $nReadWT_2=... 
  	# $nReadrrp6_1=... 
  	# $nReadrrp6_2=...
  
  	# Rscript DESeq_CUT.R "$outSeg/counts.tab" "counts_for_norm.tab" $nReadWT_1 $nReadWT_2 $nReadrrp6_1 $nReadrrp6_2 "$outSeg/"
}

export -f detect_lncRNA

## to store results
mkdir DE_new_transcripts

## parallelize execution of the function

find "$zinarDir" -name "*.gff" | xargs -P 4 -n 6 -I file bash -c "detect_lncRNA file" 

## summarize results for segmentations comparison : 
# count : 
#	- number of mRNA retrieved (75% reciprocal overlap with segment)
#	- number of new lncRNA detected

## for XUT and DUT detection
echo -e "segmentation\tmRNA\tnew_XUT\tnew_DUT\tnew_SUT" > "summary_new_lncRNA.tab"

## for CUT detection
#echo -e "segmentation\tmRNA\tnew_CUT\tnew_SUT" > "summary_new_lncRNA.tab"

ls DE_new_transcripts | while read segmentation; 
do 

	# nb of segments
	nbSeg=$(wc -l "$zinarDir/$segmentation.gff" | awk '{print $1}')
 
	# nb of mRNA retrieved 
	nbmRNA=$(awk 'OFS="\t" {if ($3 == "mRNA") print $0}' "known_transcripts.gff" | bedtools intersect -nonamecheck -s -f 0.75 -F 0.75 -wao -a stdin -b "$zinarDir/$segmentation.gff" | wc -l)
 
  	# nb of new SUT
	nbNewSUT=$(wc -l "./DE_new_transcripts/$segmentation/new_SUT.gff" | cut -d' ' -f1)
 
  	## for XUT and DUT detection
	nbNewXUT=$(wc -l "./DE_new_transcripts/$segmentation/new_XUT.gff" | cut -d' ' -f1)
	nbNewDUT=$(wc -l "./DE_new_transcripts/$segmentation/new_DUT.gff" | cut -d' ' -f1)

  	# summary
	echo -e "$segmentation\t$nbmRNA\t$nbNewXUT\t$nbNewDUT\t$nbNewSUT" >> "./summary_new_lncRNA.tab"

  	## for CUT detection
  	# nbNewCUT=$(wc -l "./DE_new_transcripts/$segmentation/new_CUT.gff" | cut -d' ' -f1)

  	# summary
	# echo -e "$segmentation\t$nbmRNA\t$nbNewXUT\t$nbNewDUT\t$nbNewSUT" >> "./summary_new_lncRNA.tab"

done

### Segmentation chosen for XUT/DUT detection is transcrits_319.gff
### Segmentation chosen for CUT detection is transcrits_159.gff

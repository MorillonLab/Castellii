#############################################################################
# Processing small RNA-seq in WT, xrn1, dicer, xrn1 dicer and Dcr1-GFP data #
#############################################################################

#####   TOOLS   ######
# bowtie2-align-s version 2.3.5
# bedtools v2.27.0
# samtools v1.9

## N. castellii genome and annotation downloaded from https://fungi.ensembl.org

## S. pombe genome and annotation downloaded from https://www.pombase.org
# Lock A, Rutherford K, Harris MA, Hayles J, Oliver SG, Bähler J; Wood V.
# PomBase 2018: user-driven reimplementation of the fission yeast database provides rapid and intuitive access to diverse, interconnected information.
# Nucleic Acids Res. 2018 (Database issue): gky961 Epub 2018 Oct 13

## genome indexes

# bowtie2 index for N. castellii & S. pombe genome
genomeIndexNcasSpom=...

# bowtie2 index for S. pombe centromeric regions
genomeIndexSpomCen=...

# bowtie2 index for N. castellii genome
genomeIndexNcas=...

## make file with annotation of S. pombe centromere
grep regional_centro schizosaccharomyces_pombe.chr.gff3 > annot_centromere_pombe.gff

## make file with annotation of N. castellii tRNA & rRNA 
awk 'OFS="\t"{if (($3 == "tRNA")||($3=="rRNA")) print $0}' Naumovozyma_castellii_cbs_4309.ASM23734v1.38.gff3 | grep -Ev "^#" > "annot_tRNA_rRNA_castellii.gff"

## chromosome names
chr_Ncas=("Chr_1" "Chr_2" "Chr_3" "Chr_4" "Chr_5" "Chr_6" "Chr_7" "Chr_8" "Chr_9" "Chr_10")
chr_Spom=("I" "II" "III")

##	fichier description_data.csv
##	index1,strain1
##	index2,strain2
##	...

description_data=description_data.csv

## for each sample
cat "$description_data" | sed '1d' | while read line; do 
  
	index=$(echo "$line" | cut -d',' -f1)
  
	fastq="$index.R1.fastq.gz"
  
  	## create mapping directory
  
  	mappingDir="$index"
  	mkdir "$mappingDir"
  
  	## for mapping stats
	mapping_log="$mappingDir/mapping_log_${index}.txt"

	## mapping with bowtie2 on N. castellii and S. pombe reference
	echo "# mapping on N. castellii and S. pombe" > "$mapping_log"
	
	(zcat "$fastq" | bowtie2 -p 8 --un "$mappingDir/$index.unmapped.Ncas.Spom.fastq" -x "$genomeIndexNcasSpom" -U - 2>> "$mapping_log") | samtools view -hF 0x4 - | samtools sort -@ 8 -o "$mappingDir/$index.mapped.all.sorted.bam" -
  	
	samtools index "$mappingDir/$index.mapped.all.sorted.bam"
  
  	# count reads mapping on N. castellii & S. pombe
  	echo -n "reads mapping on N. castellii : " >> "$mapping_log"
  	samtools view "$mappingDir/$index.mapped.all.sorted.bam" ${chr_Ncas[*]} | wc -l >> "$mapping_log"
  	
  	echo -n "reads mapping on S. pombe : " >> "$mapping_log"
  	samtools view "$mappingDir/$index.mapped.all.sorted.bam" ${chr_Spom[*]} | wc -l >> "$mapping_log"
  	
	## get uniquely mapped on N. castellii
	samtools view -h "$mappingDir/$index.mapped.all.sorted.bam" ${chr_Ncas[*]} | grep -E "^@|NM:" | grep -Ev "XS:i" | samtools view -bhS - | samtools sort -@ 8 -o "$mappingDir/$index.mapped.Ncas.unique.sorted.bam" -
	
	## get 1st base and size distribution (no tRNA and rRNA)
	## reads mapping on tRNA and rRNA excluded with bedtools
	bedtools intersect -nonamecheck -v -a "$mappingDir/$index.mapped.Ncas.unique.sorted.bam" -b "annot_tRNA_rRNA_castellii.gff" | samtools view - | awk 'OFS="\t"{if((length($10) >= 18) && (length($10) <= 30)){print $0}}' | \
  awk 'function comp(b){ 
  if (b == "A"){return "T"}; if (b == "T"){return "A"}; if(b == "C"){return "G"}; if(b == "G"){return "C"}; if(b == "N"){return "N"}
  } OFS= "\t" {if ($2 == 0) {print length($10), substr($10,1,1)} else if ($2 == 16) {print length($10), comp(substr($10,length($10),1))}}' | sort | uniq -c | awk 'OFS="\t" {print $1,$2,$3}' > "${index}_first_base_small_no_tRNA_rRNA.tab"
	
	## get small RNA 22-23bp for analysis
	samtools view -h "$mappingDir/$index.mapped.Ncas.unique.sorted.bam"  | awk 'OFS="\t"{if ($1~/@/) {print $0} else {if ((length($10) == 22) || (length($10) == 23)){print $0}}}' | samtools view -bhS - > "$mappingDir/$index.mapped.Ncas.unique.sorted.22.23.bam"

  	## index bam
	samtools index "$mappingDir/$index.mapped.Ncas.unique.sorted.22.23.bam"

	## count reads 22-23nt mapping uniquely on N. castellii 
	echo -n "reads 22-23 nt long mapping on N. castellii only : " >> "$mapping_log"
	samtools view "$mappingDir/$index.mapped.Ncas.unique.sorted.22.23.bam" | wc -l >> "$mapping_log"

	## get reads from S. pombe centromere :
	# select reads mapping on S. pombe
 	# filter out reads mapping on N. castellii
  	# map S. pombe specific reads on centromeres	
  
	samtools view -h $mappingDir/$index.mapped.all.sorted.bam ${chr_Spom[*]} | samtools view -bh | bedtools bamtofastq -i stdin -fq /dev/stdout | bowtie2 -p 8 -x "$genomeIndexNcas" -U - | samtools view -h -f 0x4 | samtools view -bh | bedtools bamtofastq -i stdin -fq /dev/stdout > "$mappingDir/$index.mapped.only.pombe.fa" 
	
	echo "# mapping on S. pombe reads on centromere" >> "$mapping_log"
	
	(cat "$mappingDir/$index.mapped.only.pombe.fa" | bowtie2 -p 8 -x "$genomeIndexSpomCen" -U - 2>> "$mapping_log") | samtools view -h -F 0x4 | samtools sort -o "$mappingDir/$index.mapped.Spom.centromere.sorted.bam"
	
	## index bam
	samtools index "$mappingDir/$index.mapped.Spom.centromere.sorted.bam"
	
	## select and count reads 22-23nt long mapping on S. pombe centromere
	samtools view -h "$mappingDir/$index.mapped.Spom.centromere.sorted.bam" | awk 'OFS="\t"{if ($1~/@/) {print $0} else {if ((length($10) == 22) || (length($10) == 23)){print $0}}}' | samtools view -bh - > "$mappingDir/$index.mapped.Spom.centromere.22.23.sorted.bam"
	
	## index bam
	samtools index "$mappingDir/$index.mapped.Spom.centromere.22.23.sorted.bam"
	
  	## count reads 22-23 nt long mapping on S. pombe centromere 
	echo -n "reads 22-23 nt long mapping on S. pombe centromere : " >> "$mapping_log"
	samtools view "$mappingDir/$index.mapped.Spom.centromere.22.23.sorted.bam" | wc -l >> "$mapping_log"
	
done

## make a file to summarize mapping stats

mapping_stats=mapping_stats_smallRNA.tab

echo -e "index\tnReads.mapped\tnReads.mapped.Ncas\tnReads.unique.Ncas\tnReads.unique.Ncas.22.23t\tnReads.mapped.Spom\tnReads.mapped.Spom.specific\tnReads.Spom.cen\tnReads.Spom.cen.22.23" > "$mapping_stats"

## chromosome names
chr_Ncas=("Chr_1" "Chr_2" "Chr_3" "Chr_4" "Chr_5" "Chr_6" "Chr_7" "Chr_8" "Chr_9" "Chr_10")
chr_Spom=("I" "II" "III" "MT" "AB325691" "mating_type_region")

cat "$description_data" | sed '1d' | while read line; do 
  
	index=$(echo "$line" | cut -d',' -f1)
	
	echo "$index"
	
	nReadsMapped=$(samtools view "$mappingDir/$index.mapped.all.sorted.bam" | wc -l) 
	nReadsMappedNcas=$(samtools view "$mappingDir/$index.mapped.all.sorted.bam" ${chr_Ncas[*]} | wc -l)
	nReadsUniqNcas=$(samtools view "$mappingDir/$index.mapped.Ncas.unique.sorted.bam" | wc -l) 
	nReadsUniqNcas22=$(samtools view "$mappingDir/$index.mapped.Ncas.unique.sorted.22.23.bam" | wc -l) 
	
	nReadsMappedSpom=$(samtools view "$mappingDir/$index.mapped.all.sorted.bam" ${chr_Spom[*]} | wc -l)
	nReadsMappedSpomSpecific=$(zcat "$mappingDir/$index.mapped.only.pombe.fa"  | wc -l)
	nReadsMappedSpomCentro=$(samtools view "$mappingDir/$index.mapped.Spom.centromere.sorted.bam" | wc -l)
	nReadsMappedSpomCentro22=$(samtools view "$mappingDir/$index.mapped.Spom.centromere.22.23.sorted.bam" | wc -l) 
	
	echo -e "$index\t$nReadsMapped\t$nReadsMappedNcas\t$nReadsUniqNcas\t$nReadsUniqNcas22\t$nReadsMappedSpom\t$nReadsMappedSpomSpecific\t$nReadsMappedSpomCentro\t$nReadsMappedSpomCentro22" >> "$mapping_stats"
	
done

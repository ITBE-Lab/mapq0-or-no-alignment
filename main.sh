#!/bin/bash
#SBATCH -p slim18 -c 18


#config
threads=18
dwgsim="~/workspace/DWGSIM/dwgsim"
bwa="bwa"
min_mapq=1
num_reads=100000

#data config
lister_427_fna="~/data/genomes/lister427/HGAP3_Tb427v10_phased_diploid_core.fasta"
lister_427_a_fna="~/data/genomes/lister427/HGAP3_Tb427v10_phased_diploid_coreA.fasta"
true927_fna="~/data/genomes/treu927/GCF_000002445.2_ASM244v1_genomic.fna"


# static
tmp="tmp"
read_prefix="$tmp/reads"
read_a_prefix="$tmp/reads_a"
read_true_prefix="$tmp/reads_true"
true927_idx="$tmp/true927_idx"
lister_427_idx="$tmp/lister_427_idx"
lister_427_a_idx="$tmp/lister_427_a_idx"
alignments_true="$tmp/alignments_true927.sam"
alignments_true927_a="$tmp/alignments_true927_a.sam"
alignments_lister="$tmp/alignments_lister427.sam"
alignments_lister_a="$tmp/alignments_lister427_a.sam"
alignments_true_self="$tmp/alignments_true927.sam"

# main
mkdir $tmp

# gen reads
eval $dwgsim -r 0 -y 0 -N $num_reads $lister_427_fna $read_prefix
eval $dwgsim -r 0 -y 0 -N $num_reads $lister_427_a_fna $read_a_prefix
eval $dwgsim -r 0 -y 0 -N $num_reads $true927_fna $read_true_prefix

# gen indices
eval $bwa index -p $true927_idx $true927_fna
eval $bwa index -p $lister_427_idx $lister_427_fna
eval $bwa index -p $lister_427_a_idx $lister_427_a_fna

# gen alignments
eval $bwa mem -t $threads $true927_idx "$read_prefix.bwa.read1.fastq.gz" "$read_prefix.bwa.read2.fastq.gz" > $alignments_true
eval $bwa mem -t $threads $true927_idx "$read_a_prefix.bwa.read1.fastq.gz" "$read_a_prefix.bwa.read2.fastq.gz" > $alignments_true927_a
eval $bwa mem -t $threads $lister_427_idx "$read_prefix.bwa.read1.fastq.gz" "$read_prefix.bwa.read2.fastq.gz" > $alignments_lister
eval $bwa mem -t $threads $lister_427_a_idx "$read_a_prefix.bwa.read1.fastq.gz" "$read_a_prefix.bwa.read2.fastq.gz" > $alignments_lister_a
eval $bwa mem -t $threads $true927_idx "$read_true_prefix.bwa.read1.fastq.gz" "$read_true_prefix.bwa.read2.fastq.gz" > $alignments_true_self

# compute percentages
num_with_good_map_q_true=$(awk '$5 > '$min_mapq' {print}' $alignments_true | wc -l)
num_with_good_map_q_true_a=$(awk '$5 > '$min_mapq' {print}' $alignments_true927_a | wc -l)
num_with_good_map_q_lister=$(awk '$5 > '$min_mapq' {print}' $alignments_lister | wc -l)
num_with_good_map_q_lister_a=$(awk '$5 > '$min_mapq' {print}' $alignments_lister_a | wc -l)
num_with_good_map_q_true_self=$(awk '$5 > '$min_mapq' {print}' $alignments_true_self | wc -l)


num_missing_true=$(samtools view -f 4 $alignments_true | wc -l)
num_missing_true_a=$(samtools view -f 4 $alignments_true927_a | wc -l)
num_missing_lister=$(samtools view -f 4 $alignments_lister | wc -l)
num_missing_lister_a=$(samtools view -f 4 $alignments_lister_a | wc -l)
num_missing_true_self=$(samtools view -f 4 $alignments_true_self | wc -l)


percent_true=$((100-100*$num_with_good_map_q_true/(2 * $num_reads)))
percent_true_a=$((100-100*$num_with_good_map_q_true_a/(2 * $num_reads)))
percent_lister=$((100-100*$num_with_good_map_q_lister/(2 * $num_reads)))
percent_lister_a=$((100-100*$num_with_good_map_q_lister_a/(2 * $num_reads)))
percent_true_self=$((100-100*$num_with_good_map_q_true_self/(2 * $num_reads)))

echo "mapq0 alignment or no alignments from LISTER 427 diploid to TRUE 927: $percent_true%"
echo "mapq0 alignment or no alignments from LISTER 427 diploid to itself: $percent_lister% <- assembly is diploid so a large portions of mapq0 is expected"
echo "mapq0 alignment or no alignments from LISTER 427 core a to TRUE 927: $percent_true_a%"
echo "mapq0 alignment or no alignments from LISTER 427 core a to itself: $percent_lister_a%"
echo "mapq0 alignment or no alignments from TRUE 927 to itself: $percent_true_self%"


percent_missing_true=$((100*$num_missing_true/(2 * $num_reads)))
percent_missing_true_a=$((100*$num_missing_true_a/(2 * $num_reads)))
percent_missing_lister=$((100*$num_missing_lister/(2 * $num_reads)))
percent_missing_lister_a=$((100*$num_missing_lister_a/(2 * $num_reads)))
percent_missing_true_self=$((100*$num_missing_true_self/(2 * $num_reads)))

echo "no allignment from LISTER 427 diploid to TRUE 927: $percent_missing_true%"
echo "no allignment from LISTER 427 diploid to itself: $percent_missing_lister%"
echo "no allignment from LISTER 427 core a to TRUE 927: $percent_missing_true_a%"
echo "no allignment from LISTER 427 core a to itself: $percent_missing_lister_a%"
echo "no allignment from TRUE 927 to itself: $percent_missing_true_self%"

# cleanup
rm -r $tmp
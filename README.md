# Percentage Lost Reads

The script in this repository evaluates the amount of reads that
 - have no alignment
 - have only alignments with mapping quality zero
 
when aligning simulated reads from one assembly to another.

In order to run this, you need to install samtools, [DWGSIM](https://github.com/nh13/DWGSIM), and [BWA](https://github.com/lh3/bwa).
Further, please set the variables in the '\#config' and '\#data config' sections at the top of of 'main.sh'.

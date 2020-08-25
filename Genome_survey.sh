###1. Data preparation
#Here uses the data of Staphylococcus aureusa from GAGE ​​for practice.
mkdir Staphylococcus_aureu && cd Staphylococcus_aureus
mkdir genome
curl ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/013/425/GCF_000013425.1_ASM1342v1/GCF_000013425.1_ASM1342v1_genomic.fna.gz> genome/Saureus.fna.gz
mkdir -p raw-data/{lib1,lib2}
curl http://gage.cbcb.umd.edu/data/Staphylococcus_aureus/Data.original/frag_1.fastq.gz> raw-data/lib1/frag_1.fastq.gz
curl http://gage.cbcb.umd.edu/data/Staphylococcus_aureus/Data.original/frag_2.fastq.gz> raw-data/lib2/frag_2.fastq.gz
curl http://gage.cbcb.umd.edu/data/Staphylococcus_aureus/Data.original/shortjump_1.fastq.gz> raw-data/lib2/shortjump_1.fastq.gz
curl http://gage.cbcb.umd.edu/data/Staphylococcus_aureus/Data.original/shortjump_2.fastq.gz> raw-data/lib2/shortjump_2.fastq.gz



###2. Quality Control
#2.1 trim_galore paired-end sequencing to remove the adapter
dir='/home/raw-data'
fq1='raw-data/lib1/frag_1.fastq.gz'
fq2='raw-data/lib1/frag_2.fastq.gz'
trim_galore -q 25 --phred33 --length 36 -e 0.1 --stringency 3 --paired -o $dir $fq1 $fq2
#dir is the result output folder. --paired For paired-end sequencing results, one read is discarded due to quality problems, and the corresponding one cannot be used. --length Length threshold, less than this length is discarded. --strency sets the number of bases that can be tolerated by the overlap between the front and rear adapters, the default is 1. -e quality control number, default is 0.1, error rate greater than 10% is discarded.

#2.2.1 fastp low-quality reads filtering and removing connectors
fastp -i data/frag_1.fastq.gz -I data/frag_2.fastq.gz -o frag_1.fq.gz -O frag_2.fq.gz
#-i specifies the input of read1, -I specifies the input of read2, -o specifies the output of read1, -O specifies the output of read2
#Generated JSON report: fastp.json and HTML report: fastp.html report is in the fastp software installation folder
#2.2.2 BFC wrong base correction
bfc -s 3m -t 16 data/frag_1.fq.gz | gzip -1> data/corrected_1.fq.gz
bfc -s 3m -t 16 data/frag_2.fq.gz | gzip -1> data/corrected_2.fq.gz
#parameter:
-s: The approximate size of the specified genome k/m/g/ means KB/MB/GB respectively
-t: indicates the number of threads



###3.K-mer analysis
#3.1 Based on software GCE
#3.1.1 kmer analysis
###k-mers generally choose 17, for highly repetitive genomes or genomes that are too large, you can choose 19 or even 31. But it is not the bigger the better, because if there is an error site in a read, the larger the k-mers will lead to an increase in the number of k-mers containing the error site.
###Use kmer_freq_hash to count k-mer frequencies
# Staphylococcus_aureus project root directory
mkdir genome_survey && cd genome_survey
#Provide read location information, write the absolute path of the fastq file (plain text fastq or fastq compressed by gz) into the text file reads.list
ls raw-data/lib1/frag_*.fastq.gz> genome_survey/reads.list
kmer_freq_hash -h
#Common parameter option description:
-k, set the length of the selected k-mer, the value range is 9-27, and the default value is 17;
-l, the name of a text document that records the absolute path of the fastq file, each line in the document is the path of a fastq file;
-c, value -1 (default), do not filter k-mer; or 0-0.99, filter low-quality k-mer based on base quality, the higher the value, the stricter;
-q, the phred format of fastq files, the default is 64, or 33 or 63 can be selected;
-r is used to filter short reads according to the set length. Reads whose length is lower than this value will not be counted. By default, all reads will be used without filtering;
-a, ignore the base of this length before reads, default none;
-d, ignore the bases of this length after reads, default none;
-i, set the size of the hash table used to store the k-mer sequence, the default is 1048576;
-p, set the prefix of the output file;
-o, whether to output k-mer sequence, 1 (default) is yes, 0 is no;
-t, the number of threads required for the program to run, single threaded by default;
-L is used to filter long reads according to the set length. Reads with a length higher than this value will not be counted. The default is 100bp.
### k-mer_freq_hash statistics, the result files "sa.freq.stat" and "kmer_freq.log" are obtained after running. "Sa.freq.stat" is the k-mer frequency statistics table, separated by tabs, the first column is the k-mer depth, that is, the frequency of each k-mer; the second column is the total number of k-mer fragments with this frequency . "Kmer_freq.log" is the log file of the program running. It counts the total number of k-mer fragments, the number of k-mer types, the average frequency of k-mer, the total number of bases, the average length of reads, and the rough estimate of genome size. Pay attention to Kmer_individual_num The data is used for the input parameters of gce.
~/opt/biosoft/gce-1.0.0/kmerfreq/kmer_freq_hash/kmer_freq_hash -k 17 -L 150 -l genome_survey/reads.list -t 10 -o 0 -p genome_survey/sa &> genome_survey/kmer_freq.log
#3.1.2 Draw a k-mer curve based on the statistical results of k-mer frequency to initially check the genome characteristics.
#R Script example
getwd()
setwd("/genome_survey")
kmer <- read.table('sa.freq.stat')
kmer <- subset(kmer, V1 >=5)
Frequency <- kmer$V1
Number <- kmer$V2
png('kmer_plot.png')
plot(Frequency, Number, type ='l', col ='blue')
dev.off()
#3.1.3 Genomic feature evaluation
gce -h
#Common parameter option description:
-f, k-mer frequency statistics table;
-c, the frequency corresponding to the "main peak" of the k-mer frequency statistics result, that is, the abscissa corresponding to the main peak of the k-mer curve;
-g, the total number of k-mer fragments. By default, the input k-mer frequency statistics table file is used to calculate the total number of k-mer;
-b, whether the data is (1) or not (0, default) has bias; when K> 19, you need to set -b 1.
-H, use heterozygous mode (1), or haploid mode (0, default);
-m, the choice of estimation model, discrete (0, default) or continuous (1);
-M, set the maximum k-mer frequency supported during calculation, the default is 256;
-D, the expected value accuracy, the default is 1; if the -m parameter selects 1, it is recommended to set the value to 8
#(1) Examples of non-hybrid mode
~/opt/biosoft/gce-1.0.0/gce -f genome_survey/sa.freq.stat -c 16 -g 108366227 -m 1 -D 8 -b 0> genome_survey/sa.table 2> genome_survey/sa .log
# -c corresponds to the depth of the main peak
# -g uses the corresponding value of Kmer_individual_num
# -m Choose estimation model, real data choose 1, which means continuous
#Result note: test.gce.stat" is the output of intermediate calculation results, and generally does not need to be concerned. "test.gce.log" is the log file of the program running, and records the evaluation results of the species genome characteristics.
raw_peak, the frequency of the "main peak" of the k-mer frequency statistics;
now_node, the number of k-mer types;
low_kmer, the number of k-mers with low coverage;
now_kmer, the total number of k-mers after filtering the number of k-mers with low coverage;
cvg, the estimated average depth of sequencing;
genome_size, the estimated genome size, genome_size = now_kmer / cvg;
a[1], the ratio of the number of k-mer species that only appears once to the total number of k-mer species;
b[1], the ratio of the number of k-mer fragments that appear only once to the total number of k-mer fragments; this value can be used to characterize the proportion of sequences with copy number 1 in the genome, then 1-b[1] can be regarded as Repetitive sequence ratio.

#(2) Examples of heterozygous patterns
#If the genome of the sequenced species is highly heterozygous (such as judged by the k-mer curve, for example, there is a heterozygous peak in the curve), we can consider running the heterozygous mode, and the results will be better than the non-heterozygous mode. . Of course, if the genome of the sequenced species is only a simple genome (no heterozygous), it cannot be run in heterozygous mode, and the result of forcibly running the heterozygous mode will be quite unreliable.
#The peak value corresponds to the position of the k-mer frequency depth. The total number of k-mer fragments is no longer specified. Instead, the input k-mer frequency statistics table file is used to calculate the total number of k-mers by default, and the hybrid mode is used.
~/opt/biosoft/gce-1.0.0/gce -f genome_survey/sa.freq.stat -c 16 -H 1 -m 1 -D 8 -b 0> genome_survey/sa.table 2> genome_survey/sa .log
#Result note: Here are the results of heterozygosity evaluation, among which:
a[1/2], the ratio of the number of heterozygous k-mer species to the total number of k-mer species; at this time, kmer-species heterozygous ratio = a[1/2] / (2- a[1/2]) , Heterozygosity rate = kmer-species heterozygous ratio / k-mer length (kmer_size);
b[1/2], the ratio of the number of heterozygous k-mer fragments to the total number of k-mer fragments; at this time, 1-b[1]-b[1/2] can be regarded as the proportion of repeated sequences.



#3.2 Based on Jellyfish+GenomeScope
#3.2.1 Run jellyfish software
jellyfish count -m 21 -s 20G -t 20 -o 21mer_out -C <(zcat test_1.fq.gz) <(zcat test_2.fq.gz)
jellyfish histo -o 21mer_out.histo 21mer_out
#Parameter comment:
# -m k-mers of K; -s Hash size, determined according to file size; -t thread; -o output prefix; -C statistics positive and negative chain
#3.2.2 Import the data into the R language for drawing
    getwd()
setwd("/genome_survey")
pdf("21_mer.out.pdf")
    dataframe19 <- read.table("21mer_out.histo")
   plot(dataframe19[1:200,], type="l")
    dev.off()
#3.2.3 Upload the obtained reads.histo to http://qb.cshl.edu/genomescope/
#parameter settings:
The predicted size of the genome is closely related to the Max kmer coverage in the first page. I set 1000 and 10000, the difference in genome is 30M. The authors explanation is that GenomeScope will filter out kmers that appear more than 1000 times by default to avoid the influence of the organelle genome. If you think the genome is small, then adjust the value to a larger value.
#Output kmer distribution map results:
Genome size (len), degree of heterozygosity (het), repeat sequence ratio (1-uniq) and other information.
#3.2.4 Use the genomescope software to map the obtained reads.histo
#Download genomescope software
git clone https://github.com/schatzlab/genomescope
#Unzip
jar xvf genomescope-master.zip
#In the installation directory of the software, the genomescopre.R file is the core running script. The usage is as follows: the first parameter kmer.hist is the kmer frequency distribution data generated by the jellyfish software, the second parameter 31 represents the length of kmer, and the third The first parameter 150 represents the sequence read length, and the fourth parameter test represents the name of the output directory.
Rscript genomescope.R kmer.hist 31 150 test
#Output result: het means heterozygosity, which is 2.36%, and len means genome size; the output file usually pays attention to the two files summary.txt and plot.png.



#3.3 Based on KmerGenie software
#Diploid mode test
#If the genome of the sequenced species is highly heterozygous (such as judging by the k-mer curve, for example, there is a heterozygous peak in the curve), we can consider running the diploid mode. At this time, the running result will be better than the haploid mode. good. Of course, if the genome of the sequenced species is only a simple genome (no heterozygosity), it cannot be run in diploid mode. Forcibly running diploid mode will not get useful results.
mkdir result_haploid
#Write the absolute path of the sequencing file into the fastq_list.txt text
ls -R /home/minminli/genomesurvey/Baci/Bacillus_subtilis.raw_R*.fastq.gz> fastq_list.txt
#--diploid, use diploid mode; choose appropriate parameters as appropriate, here start with k-mer=15, stop at 121, and 10 as the interval; the number of program running threads is 4. The result is output in result_diploid. "Diploid.log1.txt" and "diploid.log2.txt" are the correct/error output logs when the program is running. After running, you can view the result in "result_diploid". Similarly, you can focus on the result report "report.html".
kmergenie fastq_list.txt -o ./result_diploid/diploid -l 15 -k 105 -s 10 -t 4 --diploid >diploid.log1.txt 2>diploid.log2.txt
#Result interpretation: k-mer curve, the red curve is the observed k-mer curve; the blue curve is the heterozygous k-mer curve; the green curve is the homozygous k-mer curve. Single k-mer curve shape in haploid operation mode, conventional single curve display mode.
#Note one point: KmerGenie software defaults to log10 conversion of the ordinate of the k-mer frequency curve, so the k-mer curve drawn by the software by default and the conventional k-mer curve (the ordinate is not log-transformed) are visually There are intuitive differences.
#If you want KmerGenie to display the original ordinate, that is, without logarithmic conversion of the number of k-mers, you can find it in KmerGenie's drawing script (search in the path of the KmerGenie program, "kmergenie-1.7051/scripts/plot_histogram.r") to modify. Locate to line 110 in "plot_histogram.r" and delete "log='y'". #In addition, locate the 110th line, consider blocking the area with Abundance <5. Replace the command covNormalized with covNormalized[-c(1:5)] and test again. The k-mer curve here is consistent with the k-mer curve we usually see.
##Assessment of optional k-mer values ​​for second-generation genome assembly
#KmerGenie has evaluated a "best k-mer", which provides a reference for the k-mer value selected when the sequencing data is used for the formal assembly of the genome.
#According to the k-mer frequency curve obtained by KmerGenie, some other information can be obtained. For example, when the genome is assembled, the value of k-mer is affected by the sequencing depth. If the sequencing depth is higher, we can choose a higher k-mer to try assembly to obtain a longer and more complete contigs sequence; but If a higher k-mer is used for assembly in low depth sequencing mode, a higher error rate will be introduced.



###4. Genome Assembly
#4.1 SOAPdonovo2 software
#4.1.1 Prepare the configuration file config_file
#Reference configuration file: https://github.com/aquaskyline/SOAPdenovo2/blob/master/example.config
#template:
max_rd_len=100
[LIB]
avg_ins=200
reverse_seq=0
asm_flags=3
rd_len_cutoff=100
rank=1
pair_num_cutoff=3
map_len=32
q1=/home/minminli/genomesurvey/Staphylococcus_aureu/raw-data/lib1/corrected_1.fq
q2=/home/minminli/genomesurvey/Staphylococcus_aureu/raw-data/lib1/corrected_2.fq
#Template basic parameters:
#maximal read length
#This value is generally set to be slightly shorter than the actual read length, cut off the last part of the sequencing, the specific length depends on the sequencing quality
#max_rd_len=100
#文库信息Begin with this
#[LIB]
#average insert size
#Library average insertion length, generally take the library size given in the insert fragment distribution diagram
#avg_ins=200
#if sequence needs to be reversed
#Whether the sequence needs to be reverse complemented, for pair-end data, reverse complementarity is not required, set to 0; for mate-pair data, reverse complementarity is required, set to 1
#reverse_seq=0
#in which part(s) the reads are used
#1 means only contig is assembled, 2 means only scaffold is assembled, 3 means contig and scaffold are assembled at the same time, 4 means only gap is added. The short insert (<2K) is set to 3, which is used to construct contig and scaffold at the same time, the long insert (>=2k) is set to 2, not used to construct contig, only used to construct scaffold, 454single long reads are only used for complement hole.
#asm_flags=3
#use only first 100 bps of each read
#Sequence length threshold, the function is the same as max_rd_len, sequences larger than this length will be cut to this length
#rd_len_cutoff=100
#in which order the reads are used while scaffolding
#Set the priority order of different library data, the value range is an integer, and multiple libraries with the same rank value will be used at the same time when assembling the scaffold. Generally, the short insert is set to 1; 2k is set to 2; 5k is set to 3; 10k is set to 4; when a file has a large amount of data, it can also be divided into multiple files. Similarly, when a file When the amount of data is not enough, you can combine data from multiple files to build a scaffold. The amount of data mentioned here is considered from the two aspects of sequencing coverage and physical coverage of the file.
#rank=1
#cutoff of pair number for a reliable connection (at least 3 for short insert size)
#contig or the minimum overlap number before scaffold. For pair-end data, the default value is 3; for mate-paird data, the default value is 5
#pair_num_cutoff=3
#minimum aligned length to contigs for a reliable read location (at least 32 for short insert size)
#The minimum threshold of the comparison length. For pair-end data, the default value is 32; for mate-pair data, the default value is 35
#map_len=32
#a pair of fastq file, read 1 file should always be followed by read 2 file
#q1=/home/minminli/genomesurvey/Staphylococcus_aureu/raw-data/lib1/corrected_1.fq
#q2=/home/minminli/genomesurvey/Staphylococcus_aureu/raw-data/lib1/corrected_2.fq

#4.1.2 Assembly
SOAPdenovo-63mer all -K 31 -p 24 -s config_file -o soapdenovo_asm
#Basic parameters:
#-K Specify Kmer size; -p thread; -s configuration file; -o output file prefix
#Output result: *.contig, contig sequences without mate pair information; *.scafSeq, the final assembly sequence result of SOAPdenovo software, which can be used for follow-up research; *.scafStatistics, contigs and scaffolds final statistical information.
# Extract the sequence larger than 300bp.
bioawk -c fastx -v name=1'{if(length($seq)>300) print ">"name "\n" $seq;name+=1}' assembly/graph.scafSeq >soapdenovo_asm.scaffolds. fa



#4.1.3 depth_GC distribution graph
#Post the original sequence back to soapdenovo_asm.scaffolds.fa, and use samtools+bcftools to perform snp calling. Count the percentage of bases that have been altered in the total.
mkdir -p index
#Create index, index is the folder name; soapdenovo is the index prefix.
bwa index soapdenovo_asm.scaffolds.fa -p index/soapdenovo
samtools faidx ./soapdenovo_asm.scaffolds.fa
#bwa Comparison> bam
bwa mem -v 2 -t 10 index/contig read_1.fq read_2.fq | samtools view -S -b-> soapdenovo.bam
#bam sort samtools sort [options] <in.bam> <out.prefix>. In the result file, "soapdenovo.sort.bam" is the sorted comparison result file, which records the ID of each sequencing read, its alignment position in the genome, the degree of alignment, sequencing quality and other information. The so-called "sorting" is to sort the results based on the alignment of the sequencing reads and the bases in the genome, and the sequence comparison of these bases in each contigs/scaffolds sequence of the genome.
samtools sort -@ 4 soapdenovo.bam soapdenovo.sort
rm soapdenovo.bam
#Single base site sequencing depth statistics, "soapdenovo.depth" is the sequencing coverage depth of each base site in the genome obtained by Samtools statistics. The content is as follows. The first column, the sequence ID of each contigs/scaffolds in the assembled genome; the second column, each base position in the sequence, the value represents the position of the base in the sequence, and only the positions with read coverage are displayed Point; the third column, the depth of coverage of reads for each base, or the frequency at which the base at that site is sequenced. Through the content of the base site sequencing depth statistics file, we can know which sites have PCR preferential amplification and other information.
samtools depth soapdenovo.sort.bam> soapdenovo.depth
#View the comparison result bam file
samtools view soapdenovo.sort.bam | less
#View base sequencing depth statistics file
less soapdenovo.depth
#Genomic base content and sliding window statistics of sequencing reads coverage
#Sequencing depth, genomic base content sliding window statistics (here set 2000bp as a sliding window length)
#depth_base_stat.py script (custom script used in this article and an example of an important result file (extraction code hc3l): https://pan.baidu.com/s/1WYccQluDZYR079HjkVvNPA)
python3 depth_base_stat.py -g Bacillus_subtilis.scaffolds.fasta -d Bacillus_subtilis.depth -s depth_base.stat.txt -l 2000
#Output file "depth_base.stat.txt" content: the first column (seq_ID), the sequence ID of each contigs/scaffolds in the assembled genome; the second column (seq_start) and the third column (seq_end), each sliding window is in contigs/ The start/end position in the scaffolds sequence; the fourth column (depth), the average sequencing depth of each sliding window sequence; the fifth column and after (GC, A, T, G, C), within each sliding window sequence GC content percentage, A/T/G/C four base content percentages. Compared with the single-base site sequencing depth statistics file "Bacillus_subtilis.depth" obtained above, although the content of this file is reduced, its readability is enhanced, and it is easier for us to find heterogeneous amplified sequencing regions in the genome and high GC Area etc.
less -S depth_base.stat.txt
#depth_GC Scatter plot
#depth_GC_plot.r script (custom script used in this article and important result file example (extraction code hc3l): https://pan.baidu.com/s/1WYccQluDZYR079HjkVvNPA)
#Rscript depth_GC_plot.r -h
Rscript depth_GC_plot.r -i depth_base.stat.txt -o depth_GC

#4.1.4 Count the variation sites.
#Compare the sequencing reads to the assembly results, and use mutation detection software (such as GATK, etc.) to detect the "SNP site" in the genome, that is, "self call SNP". This step is modeled on the "call" in the resequencing analysis. The "SNP" idea can be used to evaluate the quality of sequencing reads in turn. Since this is not a re-sequencing analysis, the genome sequence is originally assembled from sequencing reads, so when the sequencing reads are compared back to their own assembly results, if SNP sites are detected, two situations are considered. For haploid species, SNP sites should not be detected theoretically, and the detected variant sites must be base mismatches caused by sequencing errors. For diploid or polyploid species, in addition to base sequencing errors, it may also be caused by heterozygous sites between two or more homologous chromosomes. If for non-haploid species, the detected variant type is "homozygous SNP", it can be inferred that it is caused by a sequencing error; if it is a "heterozygous SNP", it may be caused by two or more identical SNPs. The base composition between the source chromosomes is heterozygous.
samtools mpileup -d 8000 -guSDf ./soapdenovo_asm.scaffolds.fa ./soapdenovo.sort.bam | bcftools view -cvNg-> var.raw.vcf
#On the one hand, because the SOAPdenovo assembly process will make mistakes, on the other hand, samtools also has a high false positive in mutation detection, so you must first filter a batch of false positives according to depth and quality.
bcftools view -i 'DP> 30 && MQ> 30' -H variants.vcf.gz | wc -l
#Result number/genome size = mutation rate.



#4.2 ABySS software
abyss-pe k=31 name=stap in='/home/minminli/genomesurvey/Staphylococcus_aureu/raw-data/lib1/corrected_1.fq /home/minminli/genomesurvey/Staphylococcus_aureu/raw-data/lib1/corrected_2.fq '
#Basic parameters :
#-K specifies the Kmer size; name sets the output file prefix; in specifies the input file
#Multiple paired-end library assembly:
abyss-pe k=64 name=ecoli lib='lib1 lib2' lib1='lib1_1.fq lib1_2.fq' lib2='lib2_1.fq lib2_2.fq'
#Output result: stap-contigs.fa and stap-scaffolds.fa are the genome sequence after assembly; stamp-stats.csv is the statistical information of the assembly result

#4.3 Velvet Software


###5.Quast genome assembly quality assessment


###6.Bandage assembly result visualization


###7. Genome annotation


###8. Evolution Analysis



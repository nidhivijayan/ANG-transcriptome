# Two symbiotic organs and a squid
Project proposal

To identify the differences in de-novo and genome guidede transcriptome assembly


## Sample Information

|  | Nyholm Lab (unpub) | [Pankey et al. 2014](https://www.pnas.org/content/111/44/E4736.short) | [Moriano-Gutierrez et al. 2019](https://www.pnas.org/content/116/16/7990#sec-13) |
| ------ | ------ | ------ | ------ | 
| Organs | ANG (n=3) | E. scolopes ANG (n=3), and LO (n=3) | LO (n=3) |
| Time of day | 11am | Morning | Night |
| Extraction method | Trizol extraction | RiboPure kit (Ambion) | RNEasy columns (Qiagen) |
| Library prep | TruSeq Stranded mRNA Library Prep Kit +polyA selection | TruSeq Stranded mRNA Library Prep Kit | TruSeq stranded mRNA + polyA selection |
| Platform | NextSeq 500 (2x150) | Hiseq 2000 (2x150) | Hiseq 2000 (2x125) |

### I had the following objective:
1. To compare the qualities of the genome guided and de novo transcriptome assemblies

## Processing RNAseq data
All the following commands were written as batch scripts and executed on the Xanadu cluster. 

The samples sequenced by the Nyholm lab were done so in four lanes. So each forward(R1) and reverse(R2) reads were concatenated separately to one file using:
```ruby
cat E8-ANG-2_S5_L001_R1_001.fastq.gz E8-ANG-2_S5_L002_R1_001.fastq.gz E8-ANG-2_S5_L003_R1_001.fastq.gz E8-ANG-2_S5_L004_R1_001.fastq.gz > E8_total_R1_raw.fastq.gz
```

### 1. I ran trimmomatic on the fastq files:
```ruby
module load trimmomatic/0.36
java -jar $Trimmomatic PE \
 ../Concatenated_reads/E8_total_R1_raw.fastq.gz \
 ../Concatenated_reads/E8_total_R2_raw.fastq.gz \
  E8ANG_forward_paired.fq.gz E8ANG_forward_unpaired.fq.gz \
 E8ANG_reverse_paired.fq.gz E8ANG_reverse_unpaired.fq.gz \
  ILLUMINACLIP:TruSeq3-PE:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:45
```

### 2. To check the quality of the reads before and after trimming, fastqc was run on the fastq files
```ruby
module load fastqc/0.11.5
mkdir /fastqc/before_trim/
mkdir /fastqc/after_trim/
fastqc --outdir  /fastqc/before_trim/ E8_euk_R1_raw.fastq.gz
fastqc --outdir /fastqc/before_trim/ E8_euk_R2_raw.fastq.gz

fastqc --outdir  /fastqc/after_trim/ E8ANG_forward_paired_trimmed.fq.gz
fastqc --outdir /fastqc/after_trim/ E8ANG_reverse_paired_trimmed.fq.gz
```

To compare all the quality reads I ran multiqc. I moved all the zipped outputs from fastqc to a folder called zipped
```ruby
multiqc --outdir multiqc_before_after fastqc/zipped/
```
Overall the qualities between before and after trimming were not significantly different, and we did not lose many sequences to trimming.

![multiqc_plot_SeqsNumber](https://user-images.githubusercontent.com/80131639/116816286-5868ae00-ab2f-11eb-8c49-dee223ede6c4.png)

Summary of the quality metrics all the reads by MultiQC: Most of the metrics show good quality reads. We do see abnormal results for per base sequence content, sequence duplication and overrepresented sequences. However, these are normal for transcriptomic sequences.
![multiQC](https://user-images.githubusercontent.com/80131639/116816765-4d168200-ab31-11eb-9108-e4d3aa4bf93f.png)


### 3. Genome guided Transcriptome assembly
I used the genome sequenced by [Belcaid et al. 2019](https://www.pnas.org/content/116/8/3030). We first indexed the genome using hisat
```ruby
module load hisat2/2.1.0
hisat2-build eup_scolopes_assembly_v1.0.fa /hisat_genome/Es_genome
```
I then aligned the trimmed reads to the indexed genome
```ruby
module load hisat2/2.1.0
module load samtools/1.9
genome=/hisat_genome/Es_genome
hisat2 -p 8 -x ${genome} -1 ../ANG_Nyholm_trimmed_reads/E8ANG_forward_paired_trimmed.fq.gz -2 ../ANG_Nyholm_trimmed_reads/E8ANG_reverse_paired_trimmed.fq.gz -S ANG_E8.sam --rg PL:ILLUMINA --rna-strandness FR
```
I converted sam to bam files
```ruby
samtools view -S -h -u -@ 8 E8_euk_ANG.sam | samtools sort -@ 8  > ../bam/ANG_E8_euk.bam
```
Since the bam files are larged, I indexed the bam files like that of the genome. This makes it easier to search through when aligning
```ruby
samtools index ANG_E8.bam ANG_E8.bai
```
An alignment summary of the sam files were determined with samtools
```ruby
samtools flagstat /E8_euk_ANG.sam > stat_E8_ANG.txt
```
This is the summary of the alignment statistics of all the files
| SAM FILES ID | Mapped (%) | Total |
| ------ | ------ | ------ |
|J1| 77.27|58428307|
|J8|72.11|52105737|
|E16|30.85|55894617|
|E8|47.78|49674356|
|Es1|69.69|46429461|
|Es2|67.53|48746023|
|LO_Es1|71.04|49674566|
|LO_Es2|68.55|40302143|
|LO_Es3|66.68|59298719|
|WT_LO1|77.68|34823107|
|WT_LO2|68.93|44084300|
|WT_LO3|63.06|44511853|

Genome guided trinity was executed for each aligned bam file independently. I used an intron size of 6000 based on the intron size calculated from other cephalopods like the Octopus. 
```ruby
Trinity --genome_guided_bam /bam/ANG_E8_euk.bam --genome_guided_max_intron 6000 --max_memory 50G --SS_lib_type FR --output ../trinity/trinity_GG_E8_euk
```
I edited the output names of the trinity assemblies
```ruby
out=/trinity
prefix=Trinity_prefix
#sed "s/>/>E8_/g" trinity/trinity_GG_E8/Trinity-GG.fasta > /trinity/Trinity_prefix_E8.fasta  
```

I combined the output of the trinity assemblies:
```ruby
cat Trinity_prefix_* > trinity_combined_GG.fasta
```

### 4. De novo transcriptome assembly
After some trials I later learned that "forward" and "Reverse" read names are now incompatible with trinity. So we modified the files using sed to remove those words
```ruby
sed -e 's/\_forward\>//g' /Es_paired/E8ANG_forward_paired.fq > /Es_paired/E8ANG_R1.fq
```

I assembled the trimmed reads using de novo genome assembly with trinity

```ruby
module load trinity
module load samtools
input=/Es_paired/
r_p=R2
f_p=R1
Trinity --seqType fq --left ${input}/E8ANG_${f_p}.fq --right  ${input}/E8ANG_${r_p}.fq --min_contig_length 300 --max_memory 100G --output trinity_Es1ANG --full_cleanup --CPU 20
```
As it was done for the genome guided assemblies, I prefixed and concatenated the files from de novo guided assemblies.

### Identifying the coding regions

I identified the coding regions using Transdecoder
```ruby
module load TransDecoder/5.3.0
TransDecoder.LongOrfs -t /trinity_combined_denovo.fasta
```
-t is the inuput of the concatenated fasta files. I did this for genome guided (trinity_combined_GG.fasta) and de novo (shown above). This creates the folder "trinity_combine.fasta.transdecoder_dir"

Since there may be many open reading frames for the transcripts, I used pfam database to identify homology to known proteins
```ruby
hmmscan --cpu 30 --domtblout pfam.domtblout /isg/shared/databases/Trinotate/Pfam-A.hmm /trinity_combined_denovo.fasta.transdecoder_dir/longest_orfs.pep > pfam.log
```
It took about 2 days to complete hmmscan

I then use transdecoder again with the output from hmmscan, to identify which ORFs are real
```ruby
TransDecoder.Predict --cpu 16 -t trinity_combined_denovo.fasta --retain_pfam_hits pfam.domtblout
```

I then created fasta file of the representative transcripts by clustering them with vsearch
```ruby
vsearch --threads 8 --log LOGFile --cluster_fast trinity_combined_denovo.fasta.transdecoder.cds --id 0.90 --centroids centroids.fasta --uc clusters.uc
```

## 5. Compare Qualities of assemblies
For comaprison between assemblies, 5 ANGs, 2 from Nyholm and 3 from Pankey were assembled and compared. Similarly 5 LOs were assembled separately and compared. This was done due to time constraints. 

### 5.1 Bowtie2 
I examined the read composition of the assemblies using Bowtie2 to assess if the reads were properly paired. 

I first build the index of the combined fasta file"
```ruby
bowtie2-build trinity_combined_denovo.fasta trinity_combined_denovo.fasta
```

Then aligned the reads back to the indexed fasta file
```ruby
input1=/ANG_Nyholm_trimmed_reads/
input2=/Es_paired
fileList="Es1LO Es2LO Es3LO WT_LO1 WT_LO2 WT_LO2 J1ANGeuk J8ANGeuk K3ANGeuk K7ANGeuk E16_ANG2 Es1ANG Es2ANG Es3ANG"

module load bowtie2
module load samtools

for file in ${fileList}
	do
		PREFIX=`echo ${file} | cut -d "." -f 1`
		bowtie2 -p 10 -q --no-unal -k 20 -x trinity_combined_GG.fasta -1 ${input2}/${file}_forward_paired.fq -2 ${input2}/${file}_reverse_paired.fq \
			2>align_stats_${PREFIX}.txt |samtools view -@10 -Sb -o bowtie2_${PREFIX}.bam
done | tee -a log_bowtie_2.txt
```
I found all the assemblies were 100% paired with different overall alignment rate:

| ID | Overall alignment genome guided (%) | Overall alignment de novo (%) |
| ------ | ------ | ------ |
|K7ANG|89.85|92|
|K3ANG|91.78|93.6|
|E8|93.29|95.2|
|J1|90.2|92.4|
|Es1ANG|94.03|97.05|
|Es2ANG|93.54|96.3|
|Es3ANG|94.23|97.27|
|Es1LO|93|95.07|
|Es2LO|92.78|94.12|
|Es3LO|92.56|95.66|
|WT_LO1|91.5|96.09|
|WT_LO2|92.54|96.57|

### 5.2 RNAquast 

I compared the basic metrics data between genome guided and de novo assemblies of the ANG (n=5)

```ruby
module load rnaQUAST/1.5.2
rnaQUAST.py --transcripts centroids.fasta --gene_mark --threads 8 --output_dir RNAquast_denovo
```

| Metric | genome guided | de novo |
| ------ | ------ | ------ |
|Transcripts|30822|38099|
|Transcript > 500 bp|15203|17081|
|Transcripts > 1000 bp|6654|7341|
|Transcript N50|918|2637|

![rnaquast_compare](https://user-images.githubusercontent.com/80131639/116941185-36f5e800-ac3d-11eb-886a-401921396b7c.png)

There are more reads and better N50 score for the de novo assembly compared to the genome guided. This could be the result of possible bacterial contamination  However, interestingly I get more number of transcripts/isoforms of longer length in genome guided than de novo assembly. This appears to be a result of few (about 10s) longer length transcripts in genome guided. 

### 5.3 BUSCO

I compared the BUSCO scores using the Metazoa database of Busco v.4.0.2 between genome guided and de novo assemblies of ANG (n=5)
```ruby
module load busco/4.0.2
export AUGUSTUS_CONFIG_PATH=../../augustus/config

busco -i centroids.fasta -o busco_gg_metazoa \
	-m tran -l /isg/shared/databases/BUSCO/odb10/lineages/metazoa_odb10 -f
```

| Metric | Genome guided| de novo |
| ------ | ------ | ------ |
|Complete (%)| 82 | 92.3 |
|Single copy (%)| 80.2 | 90.8 |
|Duplicated (%)| 1.7 | 1.5 |
|Fragmented (%)| 11.5 | 3.5 |
|Missing(%)| 6.6 | 4.2 |

De novo guided assembly gave a better BUSCO score compare to the genome guided. 


## Work Cited:

* Alegado, R. A., Brown, L. W., Cao, S., Dermenjian, R. K., Zuzow, R., Fairclough, S. R., et al. (2012). A bacterial sulfonolipid triggers multicellular development in the closest living relatives of animals. Elife 1. doi:10.7554/eLife.00013.
* Gromek, S. M., Suria, A. M., Fullmer, M. S., Garcia, J. L., Gogarten, J. P., Nyholm, S. V, et al. (2016). Leisingera sp. JC1, a Bacterial Isolate from Hawaiian Bobtail Squid Eggs, Produces Indigoidine and Differentially Inhibits Vibrios   . Front. Microbiol.   7, 1342. Available at: https://www.frontiersin.org/article/10.3389/fmicb.2016.01342.
* Huang, J.-D., Lee, S.-Y., Chiang, T.-Y., Lu, C.-C., and Lee, M.-F. (2018). Morphology of reproductive accessory glands in female Sepia pharaonis (Cephalopoda: Sepiidae) sheds light on egg encapsulation. J. Morphol. 279, 1120–1131. doi:10.1002/jmor.20835.
* Kerwin, A. H., Gromek, S. M., Suria, A. M., Samples, R. M., Deoss, D. J., O’Donnell, K., et al. (2019). Shielding the Next Generation: Symbiotic Bacteria from a Reproductive Organ Protect Bobtail Squid Eggs from Fungal Fouling. MBio 10, e02376-19. doi:10.1128/mBio.02376-19.
* McFall-Ngai, M. J., and Ruby, E. G. (2000). Developmental biology in marine invertebrate symbioses. Curr. Opin. Microbiol. 3, 603–607. doi:https://doi.org/10.1016/S1369-5274(00)00147-8.
* Moriano-Gutierrez S, Koch EJ, Bussan H, Romano K, Belcaid M, Rey FE, Ruby EG, McFall-Ngai MJ. Critical symbiont signals drive both local and systemic changes in diel and developmental host gene expression. Proceedings of the National Academy of Sciences. 2019 Apr 16;116(16):7990-9.
* Nyholm, S. V, and McFall-Ngai, M. (2004). The winnowing: establishing the squid–vibrio symbiosis. Nat. Rev. Microbiol. 2, 632–642. doi:10.1038/nrmicro957.
* Pankey MS, Minin VN, Imholte GC, Suchard MA, Oakley TH. Predictable transcriptome evolution in the convergent and complex bioluminescent organs of squid. Proceedings of the National Academy of Sciences. 2014 Nov 4;111(44):E4736-42.
* Suria, A. M., Tan, K. C., Kerwin, A. H., Gitzel, L., Abini-Agbomson, L., Bertenshaw, J. M., et al. (2020). Hawaiian Bobtail Squid Symbionts Inhibit Marine Bacteria via Production of Specialized Metabolites, Including New Bromoalterochromides BAC-D/D′. mSphere 5, e00166-20. doi:10.1128/mSphere.00166-20.

![image](https://user-images.githubusercontent.com/80131639/110220931-ef4e1c80-7e96-11eb-97a4-ea1685980e19.png)



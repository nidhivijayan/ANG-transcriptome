# ANG-transcriptome
Project proposal

Most eukaryotes maintain a symbiotic relationship with different prokaryotes; the term used to describe an entity encompassing both host and microbes is called holobionts. The evolution of holobionts may have originated in the sea (McFall-Ngai and Ruby, 2000). For example, animal multicellularity may have emerged from marine protists. Protists like choanoflagellates need a signaling cue produced by Bacteroidetes to switch from solitary to colonial animals (Alegado et al., 2012). The Hawaiian bobtail squid (*Euprymna scolopes*) have two symbiotic organs, the light organ (LO) and the female-specific accessory nidamental gland (ANG). 

As a defense mechanism, the nocturnal bobtail squid recruits bioluminescent *Vibrio fischeri* in a specialized light organ to counterilluminate against the moonlight. The mechanism of colonization and interaction between *E. scolopes* and *V. fischeri* has been well studied and is a model system for host-microbe interactions (Nyholm and McFall-Ngai, 2004).

The ANG is a reservoir for antimicrobial bacteria (Kerwin et al., 2019). It's hypothesized that mucopolysaccharides from the nidamental gland (NG) and bacteria from the ANG coat the fertilized eggs in multiple layers of jelly coats at the connection between the NG and ANG (Fig 1B, C). The difference in post-spawning care of the cephalopods with and without an ANG supports this hypothesis; Cephalopods that lack an ANG, like octopods and some pelagic squid, guard and aerate or brood their eggs (i.e. excessive post-spawning care) whereas squid and cuttlefish that possess an ANG provide little to no post-spawn care as they lay their eggs on a rocky substrate, cover it with sand and abandon their eggs. Furthermore, antibacterial and antifungal activities were observed in eggs and bacterial isolates from the ANG and eggs of *Doryteuthis pealei* (Barbieri et al. 1997), *Uroteuthis duavaceli* (Gomathi et al. 2010), and *E. scolopes* (Gromek et al., 2016; Kerwin et al., 2019; Suria et al., 2020). Connections between the tubules containing bacteria and the nidamental gland were observed in *Sepia pharaonis* using histology (Huang et al., 2018). We do not yet know how the ANG maintains a stable microbiome.

The bobtail squid provides a unique system to study two stable symbiotic organs with different functions. Most animals, including humans, have different organs with different microbiomes; for example, gut, skin, and vaginal microbiomes. These communities are dynamic and vary with several factors, making it difficult to study their effects on the host or how the host selects for them. The light organ of *E. scolopes* selectively maintains *V. fischeri* using various immune proteins like PGRPs and chemicals like halide peroxidase and nitric oxide. We do not know the immune and metabolic response to bacteria in the ANG. 

Transcriptomics may help us understand what mechanisms the ANG uses to actively maintain and respond to bacteria. Transcriptomes can be assembled de-novo or genome-guided. There is a The is 5.1 Mb genome sequence of *E. scolopes*. The assembly has an N50 of 3,724kb with 35% gaps (Belcaid et al. 2019). This project will:
1. Compare the qualities of de novo and genome-guided ANG (n=8) and LO transcriptomes (n=5).
2. Analyze the pattern of differential expression between the ANG and LO from the de novo and genome-guided transcriptomes. If no significant difference is found, Nidhi will use the genome assembled transcriptomes for further analysis of the physiology of the ANG (to minimize issues with future authorship), while analysis of the de-novo assembled transcripts will continue for the project. If there are major differences, these will be explored. 
	* This will closely follow the tutorials for RNA-ses data analysis [Model Marine RNA-Seq and Non-ModelPlant RNA-Seq]
		* HISAT2 for the genome-guided transcriptome
		* Trinity for the de-novo assembled transcriptome
	* Comparisons between the transcriptomes will be made using mesures such as RNA-Quast and BUSCO scores to detrmine the completeness and size of the transcriptome. 
	* Differentially expressd genes will likely be called using DESeq2.
	* Differentially expressed genes will be annotated and enriched ontology terms will be identified.
3. identify and annotate orthologous genes between ANG from *E. scolopes* and *Uroteuthis edulis*. Many, but not all, squid possess an ANG, therefore, it is possible some immune genes are conserved between species with this organ. 

## Sample Information

|  | Nyholm Lab (unpub) | Pankey et al. 2014 | Moriano-Gutierrez et al. 2019 |
| ------ | ------ | ------ | ------ | 
| Organs | ANG (n=3) | E. scolopes ANG (n=3), and LO (n=3) | LO (n=3) |
| Time of day | 11am | Morning | Night |
| Extraction method | Trizol extraction | RiboPure kit (Ambion) | RNEasy columns (Qiagen) |
| Library prep | TruSeq Stranded mRNA Library Prep Kit +polyA selection | TruSeq Stranded mRNA Library Prep Kit | TruSeq stranded mRNA + polyA selection |
| Platform | NextSeq 500 (2x150) | Hiseq 2000 (2x150) | Hiseq 2000 (2x125) |

##Processing RNAseq data
All the following commands were written as batch scripts and executed on the Xanadu cluster. 

The samples sequenced by the Nyholm lab were done so in four lanes. So each forward(R1) and reverse(R2) reads were concatenated separately to one file using:
```ruby
cat E8-ANG-2_S5_L001_R1_001.fastq.gz E8-ANG-2_S5_L002_R1_001.fastq.gz E8-ANG-2_S5_L003_R1_001.fastq.gz E8-ANG-2_S5_L004_R1_001.fastq.gz > E8_total_R1_raw.fastq.gz
```

###1. We ran trimmomatic on the fastq files:
```ruby
module load trimmomatic/0.36
java -jar $Trimmomatic PE \
 ../Concatenated_reads/E8_total_R1_raw.fastq.gz \
 ../Concatenated_reads/E8_total_R2_raw.fastq.gz \
  E8ANG_forward_paired.fq.gz E8ANG_forward_unpaired.fq.gz \
 E8ANG_reverse_paired.fq.gz E8ANG_reverse_unpaired.fq.gz \
  ILLUMINACLIP:TruSeq3-PE:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:45
```

###3.To check the quality of the reads before and after trimming, fastqc was run on the fastq files
```ruby
module load fastqc/0.11.5
mkdir /fastqc/before_trim/
mkdir /fastqc/after_trim/
fastqc --outdir  /fastqc/before_trim/ E8_euk_R1_raw.fastq.gz
fastqc --outdir /fastqc/before_trim/ E8_euk_R2_raw.fastq.gz

fastqc --outdir  /fastqc/after_trim/ E8ANG_forward_paired_trimmed.fq.gz
fastqc --outdir /fastqc/after_trim/ E8ANG_reverse_paired_trimmed.fq.gz
```

To compare all the quality reads we ran multiqc. I moved all the zipped outputs from fastqc to a folder called zipped
```ruby
multiqc --outdir multiqc_before_after fastqc/zipped/
```
Overall the qualities between before and after trimming were not significantly different, and we did not lose many sequences to trimming.

![multiqc_plot_SeqsNumber](https://user-images.githubusercontent.com/80131639/116816286-5868ae00-ab2f-11eb-8c49-dee223ede6c4.png)

Summary of the quality metrics all the reads by MultiQC
![multiQC](https://user-images.githubusercontent.com/80131639/116816765-4d168200-ab31-11eb-9108-e4d3aa4bf93f.png)


###Genome guided Transcriptome assembly
We used the genome sequenced by [Belcaid et al. 2019] (https://www.pnas.org/content/116/8/3030). We first indexed the genome using hisat
```ruby
module load hisat2/2.1.0
hisat2-build eup_scolopes_assembly_v1.0.fa /hisat_genome/Es_genome
```
We then aligned the trimmed reads to the indexed genome
```ruby
module load hisat2/2.1.0
module load samtools/1.9
genome=/hisat_genome/Es_genome
hisat2 -p 8 -x ${genome} -1 ../ANG_Nyholm_trimmed_reads/E8ANG_forward_paired_trimmed.fq.gz -2 ../ANG_Nyholm_trimmed_reads/E8ANG_reverse_paired_trimmed.fq.gz -S ANG_E8.sam --rg PL:ILLUMINA --rna-strandness FR
```
We converted sam to bam files
```ruby
samtools view -S -h -u -@ 8 E8_euk_ANG.sam | samtools sort -@ 8  > ../bam/ANG_E8_euk.bam
```
Since the bam files are larged, we indexed the bam files like that of the genome. This makes it easier to search through when aligning
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
|J8|72.11|52105737
![image](https://user-images.githubusercontent.com/80131639/116817064-b4810180-ab32-11eb-81d1-af126e2dafbc.png)


Work Cited:

* Alegado, R. A., Brown, L. W., Cao, S., Dermenjian, R. K., Zuzow, R., Fairclough, S. R., et al. (2012). A bacterial sulfonolipid triggers multicellular development in the closest living relatives of animals. Elife 1. doi:10.7554/eLife.00013.
* Gromek, S. M., Suria, A. M., Fullmer, M. S., Garcia, J. L., Gogarten, J. P., Nyholm, S. V, et al. (2016). Leisingera sp. JC1, a Bacterial Isolate from Hawaiian Bobtail Squid Eggs, Produces Indigoidine and Differentially Inhibits Vibrios   . Front. Microbiol.   7, 1342. Available at: https://www.frontiersin.org/article/10.3389/fmicb.2016.01342.
* Huang, J.-D., Lee, S.-Y., Chiang, T.-Y., Lu, C.-C., and Lee, M.-F. (2018). Morphology of reproductive accessory glands in female Sepia pharaonis (Cephalopoda: Sepiidae) sheds light on egg encapsulation. J. Morphol. 279, 1120–1131. doi:10.1002/jmor.20835.
* Kerwin, A. H., Gromek, S. M., Suria, A. M., Samples, R. M., Deoss, D. J., O’Donnell, K., et al. (2019). Shielding the Next Generation: Symbiotic Bacteria from a Reproductive Organ Protect Bobtail Squid Eggs from Fungal Fouling. MBio 10, e02376-19. doi:10.1128/mBio.02376-19.
* McFall-Ngai, M. J., and Ruby, E. G. (2000). Developmental biology in marine invertebrate symbioses. Curr. Opin. Microbiol. 3, 603–607. doi:https://doi.org/10.1016/S1369-5274(00)00147-8.
* Nyholm, S. V, and McFall-Ngai, M. (2004). The winnowing: establishing the squid–vibrio symbiosis. Nat. Rev. Microbiol. 2, 632–642. doi:10.1038/nrmicro957.
* Suria, A. M., Tan, K. C., Kerwin, A. H., Gitzel, L., Abini-Agbomson, L., Bertenshaw, J. M., et al. (2020). Hawaiian Bobtail Squid Symbionts Inhibit Marine Bacteria via Production of Specialized Metabolites, Including New Bromoalterochromides BAC-D/D′. mSphere 5, e00166-20. doi:10.1128/mSphere.00166-20.

![image](https://user-images.githubusercontent.com/80131639/110220931-ef4e1c80-7e96-11eb-97a4-ea1685980e19.png)


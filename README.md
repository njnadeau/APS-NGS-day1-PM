**_Advanced Data Analysis - Introduction to NGS data analysis_**<br>
*Department of Animal and Plant Sciences, University of Sheffield*

# Aligning Illumina RNA-seq data
#### Alison Wright, Nicola Nadeau

The aim of this practical is to learn how to align Illumina RNA-seq data to a reference genome and assemble transcripts. We will be using a dataset of expression data of multiple individuals of *Heliconius melpomene*.

## Table of contents

1. Prepare reference genome

2. Align RNA-seq reads to reference

3. Visualise alignments

4. Assess mapping quality

5. Assemble transcripts

---

## PRACTICAL - Initial set up
First of all, this tutorial must be run using an interactive session in Bessemer. You will also submit jobs to Bessemer. For that, you should log in into Bessemer with `ssh`, and then request an interactive session with `srun --pty bash -l`. Your shell prompt should show `bessemer-nodeXXX` (XXX being a number between 001 and 172) and not `@bessemer-login1` nor `@bessemer-login2`.

For this particular tutorial, we are going to create and work on a directory called `align` in your /fastdata/$USER directory. Remember you need to replace $USER with your username.

        cd /fastdata/$USER
        mkdir 1.align
        cd 1.align

Remember you can check your current working directory anytime with the command `pwd`.
It should show something like:

        pwd
        /fastdata/$USER/1.align

---

## READ - Programmes

There are a number of different programs that can be used to align RNA-seq reads to a reference genome. There are two types of aligners: Splice-unaware and Splice-aware. Splice-unaware aligners are not aware of exon/intron junctions and are therefore only appropriate for mapping to a transcriptome. Splice-aware aligners map reads over exon/intro junctions and are appropriate for aligning reads to a genome reference. This is an important point to consider when choosing an aligner and why it is necessary to use different programs for RNA- versus DNA-seq data.

Potential aligners include [HISAT2](http://daehwankimlab.github.io/hisat2/manual/), [STAR](https://github.com/alexdobin/STAR) and [Tophat2](https://ccb.jhu.edu/software/tophat/index.shtml). We will be using HISAT2. This is the next generation of spliced aligner from the same research group that developed TopHat and has considerable advantages in terms of speed and memory requirements. When choosing an aligner appropriate to your study, it is useful to refer to the numerous studies that compare performance statistics across programs (eg [PeerJ review](https://peerj.com/preprints/27283/), [PLoS One review](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0190152)). The blog [RNA-Seq Blog](https://www.rna-seqblog.com/review-of-rna-seq-data-analysis-tools/) is a valuable resource.  

---

## PRACTICAL - 1. Prepare reference genome

Prepare the reference genome before mapping RNA-seq reads with [HISAT2](http://daehwankimlab.github.io/hisat2/manual/)

* First load an interactive session on Bessemer.

        srun --pty bash -l

* Make a new folder in your align folder in fastdata. Remember to change $USER to your username.

        mkdir /fastdata/$USER/1.align/ref
        cd /fastdata/$USER/1.align/ref

* Copy the reference genome to this new folder

        cp /shared/genomicsdb2/shared/workshops/NGS_AdvSta_2024/NGS_data/Reference/Hmel2.fa /fastdata/$USER/1.align/ref
	
* Have a look at the reference genome file

		head Hmel2.fa

* Index the genome using HISAT2 

        hisat2-build -f Hmel2.fa Hmel2

This command indexes the reference genome fasta file and outputs eight files: .1.ht2  .2.ht2  .3.ht2  .4.ht2  .5.ht2  .6.ht2 ... These files constitute the index.It will take around five minutes to finish. While this is running, move onto section two to familiarise yourself with HISAT2 and prepare your command to map reads to the reference genome. When `hisat2-build` has finished you can then run HISAT2.
   
## READ - 2. Align RNA-seq reads to reference

Map RNA-seq reads to a reference genome using [HISAT2](http://daehwankimlab.github.io/hisat2/manual/). We need to supply the trimmed RNA-seq reads and indexed reference genome.
 
There are many different mapping parameters you can specify, see [here](http://daehwankimlab.github.io/hisat2/manual/). While it is often suffient to run HISAT2 with default settings, there are a number of parameters that should be considered:

**-q** Reads are FASTQ files. FASTQ files usually have extension .fq or .fastq. FASTQ is the default format.

**--dta/--downstream-transcriptome-assembly** Report alignments tailored for transcript assemblers including StringTie. With this option, HISAT2 requires longer anchor lengths for de novo discovery of splice sites. This leads to fewer alignments with short-anchors, which helps transcript assemblers improve significantly in computationa and memory usage. You must specify this option when running StringTie after HISAT2 to quantify gene expression.

**-S path** File to write SAM alignments to. By default, alignments are written to the "standard out" or "stdout" filehandle (i.e. the console).

**-p**  Use this many threads to align reads. The default is 1.

**--phred33** Input qualities are ASCII chars equal to the Phred quality plus 33. This is also called the "Phred+33" encoding, which is used by the very latest Illumina pipelines.

**--met-file path** Write hisat2 metrics to file <path>. Having alignment metric can be useful for debugging certain problems, especially performance issues. See also: --met. Default: metrics disabled.

The following two parameters should be specified if it is necessary to obtain very accurate and specific mapping data. A pair of reads that align with the expected relative mate orientation and range of distances between mates is said to align “concordantly”. If both mates have unique alignments, but the alignments do not match paired-end expectations (i.e. the mates aren’t in the expected relative orientation, or aren’t within the expected distance range, or both), the pair is said to align “discordantly”. Discordant alignments may be of particular interest, for instance, when seeking structural variants.

**--no-discordant** For paired reads, report only concordant mappings.

**--no-mixed** For paired reads, only report read alignments if both reads in a pair can be mapped. By default, when hisat2 cannot find a concordant or discordant alignment for a pair, it then tries to find alignments for the individual mates. This option disables that behavior.

HISAT2 outputs a SAM file with the alignment information for each read.

## PRACTICAL - 2. Align RNA-seq reads to reference

Let's map RNA-seq reads to the reference genome. This exercise will illustrate the effect of different parameters on the stringency of mapping. In particular, we will focus on the effect of `no-mixed`. In some cases, such as for SNP calling, we may want to ensure that mapping is particularly stringent and that we have high confidence a read is correctly aligned. We can ensure higher confidence in our mapping by specifying `no-mixed`, which means that a read is only reported if both pairs align to the reference. You can compare mapping statistics between two HISAT2 runs, one with `no-mixed` and one without.

You should submit these commands as jobs to Bessemer. Our data is paired end, strand-specific and has quality scores in phred33 format. We will focus on one individual (60A). 

* First, make some output folders.

        mkdir /fastdata/$USER/1.align/Trimmed_data
        mkdir /fastdata/$USER/1.align/HISAT2
        mkdir /fastdata/$USER/1.align/HISAT2/60A
 
* Copy two fastq files into your output folder. These are subsampled versions of the files you generated yourself to speed up computational time in this practical. You must be in interactive mode `srun --pty bash -l` to do this.

        cp /shared/genomicsdb2/shared/workshops/NGS_AdvSta_2024/NGS_data/Trimmed_files/60A_1.fq.gz /fastdata/$USER/1.align/Trimmed_data
        cp /shared/genomicsdb2/shared/workshops/NGS_AdvSta_2024/NGS_data/Trimmed_files/60A_2.fq.gz /fastdata/$USER/1.align/Trimmed_data

* Make an executable script where you can specify the job requirements. 
        
        cd /fastdata/$USER/1.align/HISAT2/60A
        nano HISAT2_60A.sh
        
* Specify the requirements for Bessemer and give the path to the Genomics Software Repository. Add this to your script. 
**Indentation is important! Make sure your code looks like below.**
**Change $USER to your username, otherwise the code will not run**

```bash
#!/bin/bash
#SBATCH --time=00:15:00
#SBATCH --mem=5G
#SBATCH --cpus-per-task=2
       
source /usr/local/extras/Genomics/.bashrc
```

* Underneath, add the following commands to map reads to the reference genome.

```bash
hisat2 \
	/fastdata/$USER/1.align/ref/Hmel2 \
	-1 /fastdata/$USER/1.align/Trimmed_data/60A_1.fq.gz \
	-2 /fastdata/$USER/1.align/Trimmed_data/60A_2.fq.gz \
	-p 2 \
	-q \
	--dta \
	-S /fastdata/$USER/1.align/HISAT2/60A/60A.sam \
	--met-file /fastdata/$USER/1.align/HISAT2/60A/60A.stats
```
        
* Finally, submit your job. Only submit this if `hisat2-build` has finished and you have an indexed reference genome. This should take around 10 minutes to run so move onto the next step.
        
        sbatch HISAT2_60A.sh
	
* Let's check the job status

		squeue --me

* `R` means the job is running, `PD` is pending (waiting for resource to be available), other codes including error codes can be found [here](https://slurm.schedmd.com/squeue.html#SECTION_JOB-STATE-CODES). When the job finishes running, it will disappear from the list.
        
* Next, repeat this process to run HISAT2 again but now specifying `--no-mixed`. You need to create a new working directory and executable script. It is essential that you specify the new output folder and working directory to ensure files aren't overwritten. You need to replace `$USER` with your username.

        mkdir /fastdata/$USER/1.align/HISAT2/60A_nomixed
        cd /fastdata/$USER/1.align/HISAT2/60A_nomixed
        nano HISAT2_60A_nomixed.sh

* Add this code to your script.
        
        #!/bin/bash
        #SBATCH --time=00:15:00
        #SBATCH --mem=5G
        #SBATCH --cpus-per-task=2
               
        source /usr/local/extras/Genomics/.bashrc
        
        hisat2 \
		/fastdata/$USER/1.align/ref/Hmel2 \
		-1 /fastdata/$USER/1.align/Trimmed_data/60A_1.fq.gz \
		-2 /fastdata/$USER/1.align/Trimmed_data/60A_2.fq.gz \
		-p 2 \
		-q \
		--dta \
		--no-mixed \
		-S /fastdata/$USER/1.align/HISAT2/60A_nomixed/60A.nomixed.sam \
		--met-file /fastdata/$USER/1.align/HISAT2/60A_nomixed/60A.nomixed.stats

* Run the script. 
        
       sbatch HISAT2_60A_nomixed.sh
       
* Check the job status of both scripts. Once they have finished, let's check they ran properly.

e.g when HISAT2_60A.sh has finished there should be a files slurm-XXXXXxx.out in the /fastdata/$USER/1.align/HISAT2/60A. Look at these files with `cat` or `less`. In the slurm-XXXXXxx.out you should see stats on the mapping and the % of reads aligned at the end. **If you do not see this, something has gone wrong and you should fix it before moving on**.

--

## PRACTICAL - 3. Assess mapping quality - HISAT2 Output

HISAT2 produces various statistics including i) reads processed, ii) number of reads mapped iii) number of pairs mapped. This information will be saved in the `slurm-jobnumber.out` output file when you ran HISAT2 on Bessemer. 

Open the `slurm-jobnumber.out` output file for your first HISAT2 run using `cat` or `less`.

What is the overall alignment rate?

Now open the `slurm-jobnumber.out` output file from HISAT2 with `no-mixed`.
	
What is the overall alignment rate?

What are the differences between the files? What are the consequence of specifying `no-mixed` for read mapping?
       
---

## PRACTICAL -  4. Visualise alignments

* **Interactive Genome Viewer**

We can view read alignments to the reference genome with [IGV](http://software.broadinstitute.org/software/igv/) on the Desktop. 
	
* First, sort the SAM file with [Samtools](http://www.htslib.org/doc/samtools-1.0.html). By default Samtools sorts files by coordinate. It is easier to handle BAM files as they are compressed so at the same time we convert to a BAM format.

		samtools view -Su 60A.sam | samtools sort -o 60A.sorted.bam

* You need to then index the BAM file with [Samtools](http://www.htslib.org/doc/samtools-1.0.html).

		samtools index 60A.sorted.bam

* Let's extract reads aligning to one scaffold (Hmel200115). We can do this with the following [Samtools](http://www.htslib.org/doc/samtools-1.0.html) command.
	
	        samtools view -b -h 60A.sorted.bam "Hmel200115" > 60A_Hmel200115.bam

* You need to index the new BAM file with [Samtools](http://www.htslib.org/doc/samtools-1.0.html)
	
	        samtools index 60A_Hmel200115.bam
	
* Copy the following BAM file, index and reference genome to your desktop.
	
	        60A_Hmel200115.bam
	
	        60A_Hmel200115.bam.bai
	
	        Hmel2.fa
		
* Now you need to download the appropriate [IGV](http://software.broadinstitute.org/software/igv/) viewer onto your desktop. Download [here](http://software.broadinstitute.org/software/igv/download)
                
* Load the reference genome into IGV
	
	        Genomes/Load Genome from File
	
* Load the BAM file into IGV

	        File/Load from File

* Using the coordinates above, locate the scaffold & exons in the IGV viewer. 

		Search Hmel200115 in the search box in IGV viewer. You can also use the scaffold dropdown.

* Do reads map to these exons? Can you see the intronic regions? Is there any evidence for many reads that haven't mapped correctly?

* What do the different read colours mean? Why are there coloured lines? Can you identify any SNPs in your data?
        
---

## READ - 5. Assess mapping quality - Info from BAM files

Mapped reads can be found in the BAM file. BAM and SAM formats are designed to contain the same information. The [SAM format](https://en.wikipedia.org/wiki/SAM_(file_format)) is more human readable, and easier to process by conventional text based processing programs. The BAM format provides binary versions of most of the same data, and is designed to compress reasonably well. You can visualise a BAM file using [Samtools](http://www.htslib.org/doc/samtools-1.0.html) and get general statistics about the BAM files eg the number of mapped reads using [Samtools](http://www.htslib.org/doc/samtools-1.0.html).

## PRACTICAL - 5. Assess mapping quality - Info from BAM files

Let's find the information on read mapping from the BAM file. We can do that using the Samtools command `flagstat`.

        cd /fastdata/bo1aewr/1.align/HISAT2/60A
	samtools flagstat 60A.sorted.bam
	
How many reads mapped?
        
Next, you might also want to identify reads with high quality alignments. This may be because you want to extract a set of reads for downstream analyses. This information is encoded in the BAM/SAM file as [flags](https://samtools.github.io/hts-specs/SAMv1.pdf). The following table gives an overview of the mandatory fields in the SAM format:

![alt text](https://github.com/alielw/APS-NGS-day1-PM/blob/master/SAM%20fields.jpg)

The fields that will be most useful are the FLAGS and MAPQ. The Broad Institute have a useful website to understand the [FLAGS](https://broadinstitute.github.io/picard/explain-flags.html). In addition, the MAPQ field indicates the mapping quality of the read.

Count how many reads have all of the following features; are paired, first in the pair, on the reverse strand and also have a read mapped in proper pair. This information is encoded in the FLAG section of the SAM/BAM file. We can do this in interactive mode for the BAM file you already generated (60A.sorted.bam).

* View the first few lines of the BAM file with samtools (`samtools view BAM`), Unix pipe (`|`) and Unix `head`. 

        samtools view 60A.sorted.bam | head

* Identify which column contains the FLAG field.

* Next, print the FLAG field for the first few lines using samtools (`samtools view BAM`), Unix pipe (`|`) and cut (`cut -f[Column number]`), Unix pipe (`|`) and Unix `head`. 

* Then, using the [Broad Institute website](https://broadinstitute.github.io/picard/explain-flags.html) identify the FLAG which specifies that a read is all of the following. The FLAG will be an integer. 
        
        paired
        is mapped in proper pair
        first in pair
        on the reverse strand

* Finally, count how many reads have this FLAG using samtools (`samtools view BAM`), Unix pipe (`|`) and cut (`cut -fColumn number`), Unix pipe (`|`) and grep (`grep -c "FLAG"`).

---

## READ - 6. Assemble transcripts

We can use [StringTie](https://ccb.jhu.edu/software/stringtie/index.shtml?t=manual) to assemble gene transcripts. StringTie is part of the  “classic” RNA-Seq workflow, which includes read mapping with HISAT2 followed by assembly with StringTie. It has replaced the Cufflinks program. Although often you may have a set of annotated genes in the reference genome, and therefore a reference gtf, this may be incomplete and some genes may not be annotated. StringTie will identify potential new transcripts.

StringTie takes as input a BAM file sorted by coordinates. A text file in SAM format which was produced by HISAT2 must be sorted and converted to BAM format using the samtools program. We have already done this earlier for the IGV exercise using the command `samtools view -Su sample.sam | samtools sort -o sample.sorted.bam`.

As an option, a reference annotation file in GTF/GFF3 format can be provided to StringTie. A GTF file contains information about assembled transcripts. In this case, StringTie will prefer to use these "known" genes from the annotation file, and for the ones that are expressed it will compute coverage, TPM and FPKM values. It will also produce additional transcripts to account for RNA-seq data that aren't covered by (or explained by) the annotation. 

* We will assemble transcripts using [StringTie](https://ccb.jhu.edu/software/stringtie/index.shtml?t=manual). There are many different mapping parameters you can specify, see [here](https://ccb.jhu.edu/software/stringtie/index.shtml?t=manual). While it is often suffient to run StringTie with default settings, there are a number of parameters that should be considered:
 
**-G <ref_ann.gff>**	Use the reference annotation file (in GTF or GFF3 format) to guide the assembly process. The output will include expressed reference transcripts as well as any novel transcripts that are assembled. This option is required by options -B, -b, -e, -C (see below).

**-t**	This parameter disables trimming at the ends of the assembled transcripts. This might be important for some downstream programs. By default StringTie adjusts the predicted transcript's start and/or stop coordinates based on sudden drops in coverage of the assembled transcript.

**-e**	Limits the processing of read alignments to only estimate and output the assembled transcripts matching the reference transcripts given with the -G option (requires -G, recommended for -B/-b). You might want to use this option if you only want to quantify expression of annotated genes. With this option, read bundles with no reference transcripts will be entirely skipped, which may provide a considerable speed boost when the given set of reference transcripts is limited to a set of target genes, for example.

**-o path** Sets the name of the output GTF file where StringTie will write the assembled transcripts. This can be specified as a full path, in which case directories will be created as needed. By default StringTie writes the GTF at standard output.

**-A path** Gene abundances will be reported (tab delimited format) in the output file with the given name.

## PRACTICAL - 6. Assemble transcripts

Assemble transcripts using [StringTie](https://ccb.jhu.edu/software/stringtie/index.shtml?t=manual).

* First, make new folders.

		mkdir /fastdata/$USER/2.assemble_transcripts
		mkdir /fastdata/$USER/2.assemble_transcripts/60A

* Copy over the Heliconius gff files into this folder.

		cp /shared/genomicsdb2/shared/workshops/NGS_AdvSta_2024/NGS_data/Reference/Hmel2.gff /fastdata/$USER/2.assemble_transcripts/60A

* Let's look at the gff file. This is a list of annotated genes and transcripts in the genome with information on where they are located.
		
		cd /fastdata/$USER/2.assemble_transcripts/60A
		head Hmel2.gff

* How many line are there in the gff file? You can count this with Unix (`wc -l`).

        wc -l Hmel2.gff
	
* How many genes are there in the gff file? You can count this with Unix (`grep`).

		grep -c "gene" Hmel2.gff

* Make an executable script where you can specify the job requirements for StringTie. 
        
        nano StringTie_60A.sh
        
* Specify the requirements for Bessemer and give the path to the Genomics Software Repository. Remember
**Indentation is important! Make sure your code looks like below.**
**Change $USER to your username, otherwise the code will not run**

        #!/bin/bash
        #SBATCH --time=00:15:00
        #SBATCH --mem=5G
        #SBATCH --cpus-per-task=3
          
        source /usr/local/extras/Genomics/.bashrc

* Next, add the following commands to assemble transcripts.

        stringtie \
		/fastdata/$USER/1.align/HISAT2/60A/60A.sorted.bam \
		-G /fastdata/$USER/2.assemble_transcripts/60A/Hmel2.gff \
		-p 3 \
		-e \
		-A /fastdata/$USER/2.assemble_transcripts/60A/60A.gene_abund \
		-o /fastdata/$USER/2.assemble_transcripts/60A/60A.gtf
		
* Finally, submit your job. 
        
        sbatch StringTie_60A.sh

* Check the job status. Once it has finished, check it ran properly. slurm-XXXXXxx.out should just have two lines about the Genomics Software Repository. The script should have generated two new files `60A.gtf` and `60A.gene_abund`. We will learn about these in the next session.


---
Return to the main course page: https://github.com/njnadeau/AdvDataAna-introNGS

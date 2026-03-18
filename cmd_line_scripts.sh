#Scripts are adapted from the following sources:
#https://github.com/z0on/2bRAD_denovo/blob/master/2bRAD_README.sh
#https://catchenlab.life.illinois.edu/stacks/manual/
#https://zzz.bwh.harvard.edu/plink/data.shtml

#Authors and Contact Info:  
#David J. Bradshaw II, Ph.D.- dbradshaw2015@fau.edu, dbradshaw3366@gmail.com  
#Laura E. King, MSc- lpescitelli@fau.edu, lepescitelli@gmail.com

##Create conda environments
conda create -n sra-tools -c bioconda sra-tools
conda create -n bowtie2-samtools -c bioconda samtools bowtie2 picard -y

##Download reads
#Got to BioProject () and click on SRA, scroll down to end of list and click Send to - Run Selector.
#In the Run Selector download the Accession List and the Metadata

#Activate the environemnt
conda activate sra-tools

#Create a variable to hold the SRR Accesion list
SRA_LIST="SRR_Acc_List.txt"

#Make a reads directory
mkdir reads

# Loop through each accession number in the list
while IFS= read -r accession; do
    echo "Processing $accession"
    prefetch $accession && fasterq-dump $accession --outdir reads
done < "$SRA_LIST"

#Check have right number
ls reads/*.fastq | wc -l 
#40/40

#Download 2bRAD GitHub
git clone https://github.com/z0on/2bRAD_denovo.git 


#Make a directory to store SRAs
mkdir sras

#Move SRAs to that directory
mv SRR*/*.sra sras/

#Remove empty directories
rmdir SRR*

##Download and unzip genome file from NCBI database

#Copy fna to working directory
 cp GCA_019788955.1_BYU_Aglo_1.0_ncbi_dataset/ncbi_dataset/data/GCA_019788955.1/GCA_019788955.1_BYU_Aglo_1.0_genomic.fna .

#Save fasta and dictionary files as variables
export GENOME_FASTA=GCA_019788955.1_BYU_Aglo_1.0_genomic.fna 
export GENOME_DICT=GCA_019788955.1_BYU_Aglo_1.0_genomic.dict

#indexing genome for bowtie2 mapper
conda activate bowtie2-samtools
bowtie2-build $GENOME_FASTA $GENOME_FASTA
samtools faidx $GENOME_FASTA
java -jar /home/djbradshaw2/miniconda3/envs/bowtie2-samtools/share/picard-3.4.0-0/picard.jar CreateSequenceDictionary R=$GENOME_FASTA  O=$GENOME_DICT

#Align sequences to genome to create SAM files
#Go to folder with post-trimmed sequences
cd /home/aquaomics/bonefish/filtered/reads/stacks/post_stacks
export GENOME_FASTA=/mnt/c/Users/dbrad/Documents/HBOI/Bonefish_Revised_Manuscript/Bonefish_PA/GCA_019788955.1_BYU_Aglo_1.0_genomic.fna

#Create set of scripts to run
2bRAD_denovo/2bRAD_bowtie2_launch.pl '\.fq$' $GENOME_FASTA > bt2
2bRAD_denovo/2bRAD_bowtie2_launch.pl '\.fastq$' $GENOME_FASTA > bt2

#Execute all commands written to bt2
bash bt2

#record alignments in spreadsheet and determine the average
#/mnt/c/Users/dbrad/Documents/HBOI/Bonefish_Revised_Manuscript/Bonefish_PA/reads/alignment_rates.ods
#Average = 77.91%

#Create a list of sam files
ls *.bt2.sam > sams

#Double check you have right number of files
cat sams | wc -l  # Should be 40

#Save genome as a variable
export GENOME_FASTA=/mnt/c/Users/dbrad/Documents/HBOI/Bonefish_Revised_Manuscript/Bonefish_PA/GCA_019788955.1_BYU_Aglo_1.0_genomic.fna

#Create a list of commands to compress, sort, and index SAM files to BAM files
cat sams | perl -pe 's/(\S+)\.sam/samtools import \$GENOME_FASTA $1\.sam $1\.unsorted\.bam && samtools sort -o $1\.sorted\.bam $1\.unsorted\.bam && java -Xmx5g -jar \$TACC_PICARD_DIR\/picard\.jar AddOrReplaceReadGroups INPUT=$1\.sorted\.bam OUTPUT=$1\.bam RGID=group1 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=$1 && samtools index $1\.bam/' >s2b

#Open s2b and Ctrl + H to find and replace all $TACC_PICARD_DIR/picard.jar with /home/aquaomics/miniconda3/envs/bowtie2-samtools/share/picard-2.18.29-0/picard.jar

#Execute all command sin s2b files
bash s2b

#Get rid of extra files
rm *sorted*

#Check to see if have right number of files
ls *bam | wc -l  # Should be 40

#Move all bams to alternative folder
mv *bam /mnt/c/Users/dbrad/Documents/HBOI/Bonefish_Revised_Manuscript/Bonefish_PA/bams/

#rename bam files to remove ".fq.bt2" to make sample name be compatible with Stacks
rename 's/\.bt2.\././' *.bam

#Activate stacks envirironment
conda activate stacks_2_6_0

#Run ggstacks on bams to get SNP Genotype Call calculations
gstacks -I /mnt/c/Users/dbrad/Documents/HBOI/Bonefish_Revised_Manuscript/Bonefish_PA/bams/ -M /mnt/c/Users/dbrad/Documents/HBOI/Bonefish_Revised_Manuscript/Bonefish_PA/scripts_summaries/popmap.txt -O /mnt/c/Users/dbrad/Documents/HBOI/Bonefish_Revised_Manuscript/Bonefish_PA/stacks_out_rerun -t 2
Logging to '/mnt/c/Users/dbrad/Documents/HBOI/Bonefish_Revised_Manuscript/Bonefish_PA/stacks_out_rerun/gstacks.log'.
Locus/sample distributions will be written to '/mnt/c/Users/dbrad/Documents/HBOI/Bonefish_Revised_Manuscript/Bonefish_PA/stacks_out_rerun/gstacks.log.distribs'.

#Make a folder to dump outputs to for initial populations run
mkdir pop_R_1_rerun

#Run populations 
populations -P /mnt/c/Users/dbrad/Documents/HBOI/Bonefish_Revised_Manuscript/Bonefish_PA/stacks_out_rerun/ -M /mnt/c/Users/dbrad/Documents/HBOI/Bonefish_Revised_Manuscript/Bonefish_PA/scripts_summaries/adults_popmap.txt -R 1 --min-maf 0.1 --write-random-snp --structure --vcf --genepop --fstats --plink --hwe -t 2 --verbose --out-path /mnt/c/Users/dbrad/Documents/HBOI/Bonefish_Revised_Manuscript/Bonefish_PA/pop_R_1_rerun

#Count the number of lines in the .map file to get SNPS (subtract 1 from lines to due to header)
wc -l populations.plink.map 
#8745

#Activate plink environment
conda activate plink_1_90b6_21

#Run plink to remove any SNPs with a 90% genotyping rate (10% missing genotypes)
plink --file populations.plink --missing --geno 0.1 --out plink_geno_0.1 --recode --allow-extra-chr

#Count lines to see if SNPs changed
wc -l plink_geno_0.1.map 
#8744

#Run plink to remove any SNPs with a minimum allele frquency below 0.05
plink --file plink_geno_0.1 --missing --maf 0.05 --out plink_geno0.1_maf_0.05 --recode --allow-extra-chr

#Count lines to see if SNPs changed
wc -l plink_geno0.1_maf_0.05.map
#8744

#Run plink to remove any SNPs that do not pass hwe filter
plink --file plink_geno0.1_maf_0.05 --missing --hardy --out plink_geno0.1_maf_0.05_hwe --recode --allow-extra-chr

plink --file plink_geno0.1_maf_0.05_hwe --missing --hwe 0.1 --out plink_geno0.1_maf_0.05_hwe_0.1 --recode --allow-extra-chr

#Count lines to see if SNPs changed
wc -l plink_geno0.1_maf_0.05_hwe_0.1.map
#7671

#Remove the column with the site_locus information, replace the _ with a tab and save the list as a txt file
#This is a Whitelist i.e. a list of SNPs to keep the next time we run populations
cut -f 2 plink_geno0.1_maf_0.05_hwe_0.1.map | sed 's/_/\t/' > Adult_SNPs_whitelist.txt

#Check to make sure have correct number of SNPs
wc -l Adult_SNPs_whitelist.txt 
#7671

#Activate stacks envirironment
conda activate stacks_2_6_0

#Go back one directory
cd ..

###"Populations" Analysi Version###

#Make a new directory for Whitelist guideded analysis of male and spawn samples
mkdir WL_M_F_S1_S2_R_1_hwe_0.1

#Rerun populations with a whitelist defined by variable SNPs in males but this time include all spawn samples
populations -P /mnt/c/Users/dbrad/Documents/HBOI/Bonefish_Revised_Manuscript/Bonefish_PA/stacks_out_rerun/ -M /mnt/c/Users/dbrad/Documents/HBOI/Bonefish_Revised_Manuscript/Bonefish_PA/scripts_summaries/popmap_F_M_S1_S2.txt -W /mnt/c/Users/dbrad/Documents/HBOI/Bonefish_Revised_Manuscript/Bonefish_PA/pop_R_1_rerun/Adult_SNPs_whitelist.txt --structure --vcf --genepop --fstats --plink -t 2 --out WL_M_F_S1_S2_R_1_hwe_0.1


###Plink Analysis Version####
#Make a new directory for Whitelist guideded analysis of male and spawn samples
mkdir WL_HBOI_R_1_hwe_0.1

#Rerun populations with a whitelist defined by variable SNPs in males but this time include all spawn samples
populations -P /mnt/c/Users/dbrad/Documents/HBOI/Bonefish_Revised_Manuscript/Bonefish_PA/stacks_out_rerun/ -M /mnt/c/Users/dbrad/Documents/HBOI/Bonefish_Revised_Manuscript/Bonefish_PA/scripts_summaries/popmap_hboi.txt -W /mnt/c/Users/dbrad/Documents/HBOI/Bonefish_Revised_Manuscript/Bonefish_PA/pop_R_1_rerun/Adult_SNPs_whitelist.txt --structure --vcf --genepop --fstats --plink -t 2 --out WL_HBOI_R_1_hwe_0.1

#Change directories into the newly populated folder
cd WL_HBOI_R_1_hwe_0.1

#Activate plink environment
conda activate plink_1_90b6_21

#Use plink to create an Identity-By-State Matrix to be inputed into R
plink --file populations.plink --cluster --matrix --allow-extra-chr

#Use plink to calculate mendelian errors - Test Script
plink --file populations.plink --missing --update-sex ../update_sex.txt --update-parents ../update_parents_S11_S9.txt --mendel --out plink_me_parents_S11_S9 --recode --allow-extra-chr

#Make a directory to hold all parentage combo dataframes
mkdir mendel_parentage

#Run R code to create the tables

#Check to make sure all of them were created
ls mendel_parentage/* | wc -l 
#207/207

#Use plink to calculate mendelian errors - Loop Script
for parentage in /mnt/c/Users/dbrad/Documents/HBOI/Bonefish_Revised_Manuscript/Bonefish_PA/mendel_parentage/*; do basename=$(echo $parentage | cut -d "/" -f 11 | cut -d "_" -f 3,4,5 | cut -d "." -f 1); full_basename="$plink_basename$basename"; plink --file populations.plink --missing --update-sex ../update_sex.txt --update-parents $parentage --mendel --out $full_basename --recode --allow-extra-chr; done

#Combine all imendel files
grep "" *.imendel > all_families_imendel.tsv

#Combine all fmendel files
grep "" *.fmendel > all_families_fmendel.tsv

###Resources, Questions, Errors
##Filter based on divergence from Hardy-Weinberg
#https://groups.google.com/g/stacks-users/c/gcLDhRFrkAg

##PLINK 1.9 mendel errors trios
#https://www.biostars.org/p/278229/

##How to add trio information to VCF or PED format (to compute Mendel error in Plink)
#https://www.biostars.org/p/240870/

##Plink2 Users Google Group
#https://groups.google.com/g/plink2-users/

##Plink guide
#https://www.cog-genomics.org/plink/1.9/

##Stacks Guide
#https://catchenlab.life.illinois.edu/stacks/

##Miller-Crews et al. GitHub
#https://github.com/imillercrews/ParentageAnalysis


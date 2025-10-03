# New_caledonia_RNAseq
## orthologous detection
```bash
# jlkang@hnu2024 Thu Sep 18 2025 10:39:59 ~/HK/New_caledonia/kraken/orthologue/orthofinder_input_pep/OrthoFinder/Results_Jun28/Orthogroups
less Orthogroups.GeneCount.tsv|perl -alne 'next if /^Orthogroup/;my @a=split;my $nb;for ($i=1;$i<@a-1;$i++){$nb++ if $a[$i]>=1};print if $nb==3;'|wc -l
# 33383 orthgroups with at least one transcript

# select these orthgroups for the blastp annotation
# my $uni_pro="~/software/database/uniprot";
# jlkang@hnu2024 Thu Sep 18 2025 10:54:42 ~/HK/New_caledonia/kraken/orthologue/orthofinder_input_pep/OrthoFinder/Results_Jun28
nohup perl blastp_uni.pl Orthogroups/Orthogroups.GeneCount.tsv > run_blastp.txt 2>&1 &
# [1] 914398

####################################
# get the mapped reads number matrix
# 1. get the corresponding nucleotide sequences of orf;
# 2. map the sequences to the nucleotide sequences (RSEM)
# jlkang@hnu2024 Thu Sep 18 2025 11:34:57 ~/HK/New_caledonia/kraken/orthologue
mkdir orthofinder_nuc
# rename the headers of each nuc orf fasta
# jlkang@hnu2024 Thu Sep 18 2025 11:40:35 ~/HK/New_caledonia/kraken/orthologue/Clun_assembly.fa.transdecoder_dir
less longest_orfs.cds|perl -alne 'if (/>/){$i++;$name="Clun_$i"}else{print ">$name\n$_"}' > ../orthofinder_nuc/Clun.fa

# jlkang@hnu2024 Thu Sep 18 2025 11:41:41 ~/HK/New_caledonia/kraken/orthologue/Daru_assembly.fa.transdecoder_dir
less longest_orfs.cds|perl -alne 'if (/>/){$i++;$name="Daru_$i"}else{print ">$name\n$_"}' > ../orthofinder_nuc/Daru.fa

# jlkang@hnu2024 Thu Sep 18 2025 11:42:30 ~/HK/New_caledonia/kraken/orthologue/Zlep_assembly.fa.transdecoder_dir
less longest_orfs.cds|perl -alne 'if (/>/){$i++;$name="Zlep_$i"}else{print ">$name\n$_"}' > ../orthofinder_nuc/Zlep.fa

# build the bowtie2 index
# nohup rsem-prepare-reference --bowtie2 Blenny.fa Blenny --bowtie2 &
# jlkang@hnu2024 Thu Sep 18 2025 11:48:30 ~/HK/New_caledonia/kraken/orthologue/orthofinder_nuc
nohup rsem-prepare-reference --bowtie2 Clun.fa Clun --bowtie2 > Clun_build_ref_bowtie2.log 2>&1 &
# [2] 1281372
nohup rsem-prepare-reference --bowtie2 Daru.fa Daru --bowtie2 > Daru_build_ref_bowtie2.log 2>&1 &
# [3] 1283291
nohup rsem-prepare-reference --bowtie2 Zlep.fa Zlep --bowtie2 > Zlep_build_ref_bowtie2.log 2>&1 &
# [4] 1285392

# Obtain the sample information
# jlkang@hnu2024 Thu Sep 18 2025 15:40:39 ~/HK/New_caledonia/kraken
ll Clun*_1.fq.gz|perl -alne 'my @a=split;(my $nm)=$a[-1]=~/(.*?)_1\.fq\.gz/;print "$nm"' > ~/HK/New_caledonia/kraken/orthologue/orthofinder_nuc/Clun_sample.txt
ll Daru*_1.fq.gz|perl -alne 'my @a=split;(my $nm)=$a[-1]=~/(.*?)_1\.fq\.gz/;print "$nm"' > ~/HK/New_caledonia/kraken/orthologue/orthofinder_nuc/Daru_sample.txt
ll Zlep*_1.fq.gz|perl -alne 'my @a=split;(my $nm)=$a[-1]=~/(.*?)_1\.fq\.gz/;print "$nm"' > ~/HK/New_caledonia/kraken/orthologue/orthofinder_nuc/Zlep_sample.txt
```

```run_rsem.pl
#!/usr/bin/perl
use strict;
use warnings;

my $spe=$ARGV[0];
my $sample=$ARGV[1];
my $dir="~/HK/New_caledonia/kraken";

my $result_dir=$spe."_RSEM_output";
mkdir $result_dir unless (-d $result_dir);
chdir "./$result_dir";

open FIL, "../$sample" or die "can not open ../$sample\n";
while (<FIL>) {
	chomp;
    my $re=$_.".genes.results";
    #print "$re\n";
    if (-e $re) {
        next;
    } else {
        my $R1=$_."_1.fq.gz";
        my $R2=$_."_2.fq.gz";
        my $cmd="rsem-calculate-expression -p 30 --bowtie2 --paired-end $dir/$R1 $dir/$R2 ../$spe $_ --no-bam-output";
       #print "$cmd\n";
		system($cmd);
		chdir "../";
    }
}
```


```bash
# jlkang@hnu2024 Thu Sep 18 2025 16:05:53 ~/HK/New_caledonia/kraken/orthologue/orthofinder_nuc
nohup perl run_rsem.pl Clun Clun_sample.txt > Clun_rsem.log 2>&1 &
# [1] 845198
nohup perl run_rsem.pl Daru Daru_sample.txt > Daru_rsem.log 2>&1 &
# [2] 846962
nohup perl run_rsem.pl Zlep Zlep_sample.txt > Zlep_rsem.log 2>&1 &
# [3] 848151

# Summary the reads number
# Clun
# jlkang@hnu2024 Fri Sep 19 2025 12:17:59 ~/HK/New_caledonia/kraken/orthologue/orthofinder_nuc/Clun_RSEM_output
ll *genes.results|perl -alne '$cmd="../merge_RSEM_frag_counts_single_table.pl ";$cmd1.=$F[-1]." ";END{$cmd1=~s/\s+$//;$cmd2=$cmd.$cmd1;`$cmd2 > total_Clun.gene.matrix`}'
less total_Clun.gene.matrix|perl -alne 's/\.genes\.results//g;print' > total_Clun.gene.matrix.1; mv total_Clun.gene.matrix.1 total_Clun.gene.matrix
# Daru
# jlkang@hnu2024 Fri Sep 19 2025 13:27:39 ~/HK/New_caledonia/kraken/orthologue/orthofinder_nuc/Daru_RSEM_output
ll *genes.results|perl -alne '$cmd="../merge_RSEM_frag_counts_single_table.pl ";$cmd1.=$F[-1]." ";END{$cmd1=~s/\s+$//;$cmd2=$cmd.$cmd1;`$cmd2 > total_Daru.gene.matrix`}'
less total_Daru.gene.matrix|perl -alne 's/\.genes\.results//g;print' > total_Daru.gene.matrix.1; mv total_Daru.gene.matrix.1 total_Daru.gene.matrix
# Zlep
# jlkang@hnu2024 Fri Sep 19 2025 13:30:00 ~/HK/New_caledonia/kraken/orthologue/orthofinder_nuc/Zlep_RSEM_output
ll *genes.results|perl -alne '$cmd="../merge_RSEM_frag_counts_single_table.pl ";$cmd1.=$F[-1]." ";END{$cmd1=~s/\s+$//;$cmd2=$cmd.$cmd1;`$cmd2 > total_Zlep.gene.matrix`}'
less total_Zlep.gene.matrix|perl -alne 's/\.genes\.results//g;print' > total_Zlep.gene.matrix.1; mv total_Zlep.gene.matrix.1 total_Zlep.gene.matrix

# jlkang@hnu2024 Fri Sep 19 2025 13:32:08 ~/HK/New_caledonia/kraken/orthologue/orthofinder_input_pep/OrthoFinder/Results_Jun28
mkdir reads_matrix; cd reads_matrix
# jlkang@hnu2024 Fri Sep 19 2025 13:33:08 ~/HK/New_caledonia/kraken/orthologue/orthofinder_input_pep/OrthoFinder/Results_Jun28/reads_matrix
cp ~/HK/New_caledonia/kraken/orthologue/orthofinder_nuc/Clun_RSEM_output/total_Clun.gene.matrix ./
cp ~/HK/New_caledonia/kraken/orthologue/orthofinder_nuc/Daru_RSEM_output/total_Daru.gene.matrix ./
cp ~/HK/New_caledonia/kraken/orthologue/orthofinder_nuc/Zlep_RSEM_output/total_Zlep.gene.matrix ./
# jlkang@hnu2024 Fri Sep 19 2025 13:37:25 ~/HK/New_caledonia/kraken/orthologue/orthofinder_input_pep/OrthoFinder/Results_Jun28
# change "my @specs" and the species number to "3" ($info1_1{$key}->{num_spec}==2) in get_best_blast_orthogroup.pl
nohup perl get_best_blast_orthogroup.pl -blast_result=blastp_result -reads_matrix=reads_matrix > get_best_blast_orthogroup.log 2>&1 &
# [1] 2296260
# jlkang@hnu2024 Fri Sep 19 2025 16:48:38 ~/HK/New_caledonia/kraken/orthologue/orthofinder_input_pep/OrthoFinder/Results_Jun28/final_blast_orth_group
ll *blastp_result|wc -l # 20936 orthogroups

# Create the references
# jlkang@hnu2024 Fri Sep 19 2025 16:59:28 ~/HK/New_caledonia/kraken/orthologue/orthofinder_input_pep/OrthoFinder/Results_Jun28
mkdir orthofinder_nuc; cd orthofinder_nuc
# jlkang@hnu2024 Fri Sep 19 2025 17:00:36 ~/HK/New_caledonia/kraken/orthologue/orthofinder_input_pep/OrthoFinder/Results_Jun28/orthofinder_nuc
cp ~/HK/New_caledonia/kraken/orthologue/orthofinder_nuc/Clun.fa ./Clun_nuc.fa
cp ~/HK/New_caledonia/kraken/orthologue/orthofinder_nuc/Daru.fa ./Daru_nuc.fa
cp ~/HK/New_caledonia/kraken/orthologue/orthofinder_nuc/Zlep.fa ./Zlep_nuc.fa

# jlkang@hnu2024 Fri Sep 19 2025 17:02:13 ~/HK/New_caledonia/kraken/orthologue/orthofinder_input_pep/OrthoFinder/Results_Jun28
perl get_sequences_ref.pl -input=final_blast_orth_group -nuc=orthofinder_nuc -output=final_reference
cp -r final_reference/ ~/HK/New_caledonia/kraken/
# RSEM again
# jlkang@hnu2024 Fri Sep 19 2025 17:08:35 ~/HK/New_caledonia/kraken/final_reference
nohup rsem-prepare-reference --bowtie2 Clun.fa Clun --bowtie2 > Clun_build_ref_bowtie2.log 2>&1 &
nohup rsem-prepare-reference --bowtie2 Daru.fa Daru --bowtie2 > Daru_build_ref_bowtie2.log 2>&1 &
nohup rsem-prepare-reference --bowtie2 Zlep.fa Zlep --bowtie2 > Zlep_build_ref_bowtie2.log 2>&1 &

# jlkang@hnu2024 Fri Sep 19 2025 17:11:12 ~/HK/New_caledonia/kraken/final_reference
cp ~/HK/New_caledonia/kraken/orthologue/orthofinder_nuc/*.pl ./
cp ~/HK/New_caledonia/kraken/orthologue/orthofinder_nuc/*_sample.txt ./

# jlkang@hnu2024 Fri Sep 19 2025 17:14:03 ~/HK/New_caledonia/kraken/final_reference
nohup perl run_rsem.pl Clun Clun_sample.txt > Clun_rsem.log 2>&1 &
# [1] 6705
nohup perl run_rsem.pl Daru Daru_sample.txt > Daru_rsem.log 2>&1 &
# [2] 6777
nohup perl run_rsem.pl Zlep Zlep_sample.txt > Zlep_rsem.log 2>&1 &
# [3] 6856

# Summary the reads number
# Clun
# jlkang@hnu2024 Sat Sep 20 2025 08:18:32 ~/HK/New_caledonia/kraken/final_reference/Clun_RSEM_output
ll *genes.results|perl -alne '$cmd="../merge_RSEM_frag_counts_single_table.pl ";$cmd1.=$F[-1]." ";END{$cmd1=~s/\s+$//;$cmd2=$cmd.$cmd1;`$cmd2 > total_Clun.gene.matrix`}'
less total_Clun.gene.matrix|perl -alne 's/\.genes\.results//g;print' > total_Clun.gene.matrix.1; mv total_Clun.gene.matrix.1 total_Clun.gene.matrix

# Daru
# jlkang@hnu2024 Sat Sep 20 2025 08:20:28 ~/HK/New_caledonia/kraken/final_reference/Daru_RSEM_output
ll *genes.results|perl -alne '$cmd="../merge_RSEM_frag_counts_single_table.pl ";$cmd1.=$F[-1]." ";END{$cmd1=~s/\s+$//;$cmd2=$cmd.$cmd1;`$cmd2 > total_Daru.gene.matrix`}'
less total_Daru.gene.matrix|perl -alne 's/\.genes\.results//g;print' > total_Daru.gene.matrix.1; mv total_Daru.gene.matrix.1 total_Daru.gene.matrix

# Zlep
# jlkang@hnu2024 Sat Sep 20 2025 08:21:43 ~/HK/New_caledonia/kraken/final_reference/Zlep_RSEM_output
ll *genes.results|perl -alne '$cmd="../merge_RSEM_frag_counts_single_table.pl ";$cmd1.=$F[-1]." ";END{$cmd1=~s/\s+$//;$cmd2=$cmd.$cmd1;`$cmd2 > total_Zlep.gene.matrix`}'
less total_Zlep.gene.matrix|perl -alne 's/\.genes\.results//g;print' > total_Zlep.gene.matrix.1; mv total_Zlep.gene.matrix.1 total_Zlep.gene.matrix

# summary the information of sequenced samples
# kangjingliang@KangdeMacBook-Pro-2 六  9 20 2025 16:50:49 ~/Documents/2025/New_caledonia
perl temp1.pl > sample_seq_information.txt
```

## DEGs detection
```bash
# Clun
# Obtain the reads number matrix according to the sample order in coldata_Clun_brain.txt
# kangjingliang@KangdeMacBook-Pro-2 六  9 20 2025 17:19:47 ~/Documents/2025/New_caledonia/Clun
extract_reads_nb --matrix ../Reads_nb_matrix/total_Clun.gene.matrix --samples coldata_Clun_brain.txt > Clun_brain_reads_matrix.xls
# int the reads number
less Clun_brain_reads_matrix.xls|perl -alne 'if (/^\s+/){print}else{my $info;for($i=1;$i<@F;$i++){$info.=int($F[$i])."\t"};$info=~s/\s+$//;print"$F[0]\t$info"}' > Clun_brain_reads_matrix.xls.1; mv Clun_brain_reads_matrix.xls.1 Clun_brain_reads_matrix.xls

# Daru
# kangjingliang@KangdeMacBook-Pro-2 六  9 20 2025 17:44:39 ~/Documents/2025/New_caledonia/Daru
# brain
extract_reads_nb --matrix ../Reads_nb_matrix/total_Daru.gene.matrix --samples coldata_Daru_brain.txt > Daru_brain_reads_matrix.xls
less Daru_brain_reads_matrix.xls|perl -alne 'if (/^\s+/){print}else{my $info;for($i=1;$i<@F;$i++){$info.=int($F[$i])."\t"};$info=~s/\s+$//;print"$F[0]\t$info"}' > Daru_brain_reads_matrix.xls.1; mv Daru_brain_reads_matrix.xls.1 Daru_brain_reads_matrix.xls
# gill
extract_reads_nb --matrix ../Reads_nb_matrix/total_Daru.gene.matrix --samples coldata_Daru_gill.txt > Daru_gill_reads_matrix.xls
less Daru_gill_reads_matrix.xls|perl -alne 'if (/^\s+/){print}else{my $info;for($i=1;$i<@F;$i++){$info.=int($F[$i])."\t"};$info=~s/\s+$//;print"$F[0]\t$info"}' > Daru_gill_reads_matrix.xls.1; mv Daru_gill_reads_matrix.xls.1 Daru_gill_reads_matrix.xls

# Zlep
# kangjingliang@KangdeMacBook-Pro-2 六  9 20 2025 17:54:15 ~/Documents/2025/New_caledonia/Zlep
# brain
extract_reads_nb --matrix ../Reads_nb_matrix/total_Zlep.gene.matrix --samples coldata_Zlep_brain.txt > Zlep_brain_reads_matrix.xls
less Zlep_brain_reads_matrix.xls|perl -alne 'if (/^\s+/){print}else{my $info;for($i=1;$i<@F;$i++){$info.=int($F[$i])."\t"};$info=~s/\s+$//;print"$F[0]\t$info"}' > Zlep_brain_reads_matrix.xls.1; mv Zlep_brain_reads_matrix.xls.1 Zlep_brain_reads_matrix.xls
# gill
extract_reads_nb --matrix ../Reads_nb_matrix/total_Zlep.gene.matrix --samples coldata_Zlep_gill.txt > Zlep_gill_reads_matrix.xls
less Zlep_gill_reads_matrix.xls|perl -alne 'if (/^\s+/){print}else{my $info;for($i=1;$i<@F;$i++){$info.=int($F[$i])."\t"};$info=~s/\s+$//;print"$F[0]\t$info"}' > Zlep_gill_reads_matrix.xls.1; mv Zlep_gill_reads_matrix.xls.1 Zlep_gill_reads_matrix.xls

# PCA
# kangjingliang@KangdeMacBook-Pro-2 一  9 22 2025 10:36:33 ~/Documents/2025/New_caledonia/Clun
mpca --matrix Clun_brain_reads_matrix.xls --samples coldata_Clun_brain.txt --column site --title Clun --label --prefix Clun_all_gene
# kangjingliang@KangdeMacBook-Pro-2 一  9 22 2025 10:33:55 ~/Documents/2025/New_caledonia/Daru
mpca --matrix Daru_brain_reads_matrix.xls --samples coldata_Daru_brain.txt --column site --title Daru_brain --label --prefix Daru_brain_all_gene
mpca --matrix Daru_gill_reads_matrix.xls --samples coldata_Daru_gill.txt --column site --title Daru_gill --label --prefix Daru_gill_all_gene
# kangjingliang@KangdeMacBook-Pro-2 一  9 22 2025 10:37:09 ~/Documents/2025/New_caledonia/Zlep
mpca --matrix Zlep_brain_reads_matrix.xls --samples coldata_Zlep_brain.txt --column site --title Zlep_brain --label --prefix Zlep_brain_all_gene
mpca --matrix Zlep_gill_reads_matrix.xls --samples coldata_Zlep_gill.txt --column site --title Zlep_gill --label --prefix Zlep_gill_all_gene

# DEseq
# Obtain the orthogroup annotation
# ~/software/database/swiss_pro_info.txt
# jlkang@hnu2024 Mon Sep 22 2025 16:06:28 ~/HK/New_caledonia/kraken/orthologue/orthofinder_input_pep/OrthoFinder/Results_Jun28
perl Extract_anno.pl > unprot_name_description_orthgroup.txt

# Clun
# kangjingliang@KangdeMacBook-Pro-2 一  9 22 2025 16:14:15 ~/Documents/2025/New_caledonia/Clun
DESeq --matrix Clun_brain_reads_matrix.xls --samples coldata_Clun_brain.txt --column site --prefix Clun_brain
# 7 DEGs
extract_anno --genes Clun_brain_Bourake_Control.DEGs.txt --anno ../unprot_name_description_orthgroup.txt --col 1 >Clun_brain_Bourake_Control.DEGs.ano.txt

# Daru
# brain
# kangjingliang@KangdeMacBook-Pro-2 一  9 22 2025 16:18:27 ~/Documents/2025/New_caledonia/Daru
DESeq --matrix Daru_brain_reads_matrix.xls --samples coldata_Daru_brain.txt --column site --prefix Daru_brain
# 779 DEGs
extract_anno --genes Daru_brain_Bourake_Control.DEGs.txt --anno ../unprot_name_description_orthgroup.txt --col 1 >Daru_brain_Bourake_Control.DEGs.ano.txt
# gill
DESeq --matrix Daru_gill_reads_matrix.xls --samples coldata_Daru_gill.txt --column site --prefix Daru_gill
# 2000 DEGs
extract_anno --genes Daru_gill_Bourake_Control.DEGs.txt --anno ../unprot_name_description_orthgroup.txt --col 1 >Daru_gill_Bourake_Control.DEGs.ano.txt

# Zlep
DESeq --matrix Zlep_brain_reads_matrix.xls --samples coldata_Zlep_brain.txt --column site --prefix Zlep_brain
# 526 DEGs
extract_anno --genes Zlep_brain_Bourake_Control.DEGs.txt --anno ../unprot_name_description_orthgroup.txt --col 1 >Zlep_brain_Bourake_Control.DEGs.ano.txt
# gill
DESeq --matrix Zlep_gill_reads_matrix.xls --samples coldata_Zlep_gill.txt --column site --prefix Zlep_gill
# 2443 DEGs
extract_anno --genes Zlep_gill_Bourake_Control.DEGs.txt --anno ../unprot_name_description_orthgroup.txt --col 1 >Zlep_gill_Bourake_Control.DEGs.ano.txt

# PCA including all individuals
# kangjingliang@KangdeMacBook-Pro-2 五  9 26 2025 13:31:14 ~/Documents/2025/New_caledonia/Reads_nb_matrix
perl temp1.pl > All_sample_reads_nb.txt
# kangjingliang@KangdeMacBook-Pro-2 五  9 26 2025 13:49:45 ~/Documents/2025/New_caledonia/Reads_nb_matrix
less All_sample_reads_nb.txt|head -n 1|perl -alne '@a=split /\t/;for (my $i=1; $i<@a; $i++){print $a[$i]}' > sample_info.txt
# kangjingliang@KangdeMacBook-Pro-2 五  9 26 2025 14:25:50 ~/Documents/2025/New_caledonia/Reads_nb_matrix
perl temp2.pl > sample_info3.txt
# use sample_info3.txt
```

## Enrichment analysis
```bash
# obtain the sequences of mapped uniprot ID for the functional enrichment
# jlkang@hnu2024 Tue Sep 23 2025 08:28:20 ~/HK/New_caledonia/kraken/orthologue/orthofinder_input_pep/OrthoFinder/Results_Jun28
perl create_orth_seq.pl > orth_seq.fasta
# Omicsbox
```

## The common and species-specific DEGs
### The common DEGs between tissues of different species
```bash
# kangjingliang@KangdeMacBook-Pro-2 日  9 28 2025 21:57:54 ~/Documents/2025/New_caledonia/Reads_nb_matrix
# Common_DEGs_gill_DaruZlep.txt: 370
# Common_DEGs_brain_DaruZlep.txt: 90
perl temp3.pl Common_DEGs_brain_DaruZlep.txt > Common_DEGs_brain_DaruZlep_ano.txt
perl temp3.pl Common_DEGs_gill_DaruZlep.txt > Common_DEGs_gill_DaruZlep_ano.txt

# put all common DEGs between tissues for functional enrichment
# kangjingliang@KangdeMacBook-Pro-2 一  9 29 2025 11:31:48 ~/Documents/2025/New_caledonia/Reads_nb_matrix/Common_enrich_gill
less Common_gill_fisher.txt|perl -alne '@a=split /\t/;print if ($a[3] eq "BIOLOGICAL_PROCESS" && $a[4] <=0.05)'
# 36 significantly enriched functions: circadian rhythm, glucocorticoid

# kangjingliang@KangdeMacBook-Pro-2 一  9 29 2025 11:41:40 ~/Documents/2025/New_caledonia/Reads_nb_matrix/Common_enrich_brain
less Common_brain_fisher.txt|perl -alne '@a=split /\t/;print if ($a[3] eq "BIOLOGICAL_PROCESS" && $a[4] <=0.05)'
# 30 significantly enriched functions: circadian rhythm, glucocorticoid, photoperiodism, blue light signaling pathway, cellular response to blue light, entrainment of circadian clock by photoperiod
```
### The species-specific DEGs in each species
```bash
# Daru
# kangjingliang@KangdeMacBook-Pro-2 一  9 29 2025 17:24:24 ~/Documents/2025/New_caledonia/Daru/Enrichment_brain
perl Add_genes.pl Daru_brain_fisher_GOseq.txt Daru_brain_fisher.txt > DaruBrain_enrichment.txt
# kangjingliang@KangdeMacBook-Pro-2 一  9 29 2025 17:26:22 ~/Documents/2025/New_caledonia/Daru/Enrichment_gill
perl Add_genes.pl Daru_gill_fisher_GOseq.txt Daru_gill_fisher.txt > DaruGill_enrichment.txt

# Zlep
# kangjingliang@KangdeMacBook-Pro-2 一  9 29 2025 17:29:18 ~/Documents/2025/New_caledonia/Zlep/Enrichment_brain
perl Add_genes.pl Zlep_brain_fisher_GOseq.txt Zlep_brain_fisher.txt > ZlepBrain_enrichment.txt
# kangjingliang@KangdeMacBook-Pro-2 一  9 29 2025 17:30:16 ~/Documents/2025/New_caledonia/Zlep/Enrichment_gill
perl Add_genes.pl Zlep_gill_fisher_GOseq.txt Zlep_gill_fisher.txt > ZlepGill_enrichment.txt

# kangjingliang@KangdeMacBook-Pro-2 一  9 29 2025 17:46:20 ~/Documents/2025/New_caledonia/Daru/Enrichment_brain
extract_gene_functions -i DaruBrain_enrichment.txt -a ../../unprot_name_description_orthgroup.txt --gene_column 1 --func_column 3 --functions light_func.txt --output light_func_DEGs

# common DEGs
# kangjingliang@KangdeMBP-2 三 10 01 2025 10:17:20 ~/Documents/2025/New_caledonia/Reads_nb_matrix/Common_enrich_brain
perl Add_genes.pl Common_brain_fisher_GOseq.txt Common_brain_fisher.txt > CommonBrain_enrichment.txt
extract_gene_functions -i CommonBrain_enrichment.txt -a ../../unprot_name_description_orthgroup.txt --gene_column 1 --func_column 3 --functions circadian_functions.txt --output circadian_func_DEGs

# kangjingliang@KangdeMBP-2 三 10 01 2025 10:18:45 ~/Documents/2025/New_caledonia/Reads_nb_matrix/Common_enrich_gill
perl Add_genes.pl Common_gill_fisher_GOseq.txt Common_gill_fisher.txt > CommonGill_enrichment.txt
extract_gene_functions -i CommonGill_enrichment.txt -a ../../unprot_name_description_orthgroup.txt --gene_column 1 --func_column 3 --functions circadian_functions.txt --output circadian_func_DEGs
```


```
# download and install Trinity
# jlkang@hnu2024 Fri Sep 19 2025 08:34:45 ~/software
wget https://github.com/trinityrnaseq/trinityrnaseq/releases/download/Trinity-v2.15.2/trinityrnaseq-v2.15.2.FULL.tar.gz
tar -zxvf trinityrnaseq-v2.15.2.FULL.tar.gz
cd trinityrnaseq-v2.15.2
# jlkang@hnu2024 Fri Sep 19 2025 08:48:30 ~/software/trinityrnaseq-v2.15.2
make
# error
sudo apt-get install autoconf automake libtool
# still have error
# try to use it in docker
# jlkang@hnu2024 Fri Sep 19 2025 08:58:03 ~/software/trinityrnaseq-v2.15.2/trinity-plugins/bamsifter
sudo docker pull trinityrnaseq/trinityrnaseq
```

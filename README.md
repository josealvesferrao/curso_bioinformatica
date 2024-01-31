# curso_bioinformatica


## Unix tutorial / Command Line 
Open a web browser and find this github repository https://github.com/krother/bash_tutorial

go to your home directory 
```
cd ~
```
```
git clone https://github.com/krother/bash_tutorial.git
```

Additional tutorial: https://ubuntu.com/tutorials/command-line-for-beginners#1-overview


## FASTQC

```
fastqc -t 1 -o /mnt/sdb/curso_bioinformatica/output/grupo1/fastqc /mnt/sdb/curso_bioinformatica/raw_data/sample1/selected_fastq/sample1_chr2-171000000-172000000_R1.fastq.gz sample1_chr2-171000000-172000000_R2.fastq.gz
```

##  Alignment

```
bwa mem -t 1 -R '@RG\tID:sample1\tPU:NB551361.1.sample1\tSM:sample1\tPL:ILLUMINA\tLB:sample1\tCN:curso' /mnt/sdb/curso_bioinformatica/raw_data/fasta_ref_genome/hg38/Homo_sapiens.GRCh38.dna.primary_assembly.par_y_n_masked.chr2.fasta /mnt/sdb/curso_bioinformatica/raw_data/sample1/selected_fastq/sample1_chr2-171000000-172000000_R1.fastq.gz /mnt/sdb/curso_bioinformatica/raw_data/sample1/selected_fastq/sample1_chr2-171000000-172000000_R2.fastq.gz > sample1_chr2-171000000-172000000.sam
```

## Variant calling

```
gatk HaplotypeCaller --reference /mnt/sdb/curso_bioinformatica/raw_data/fasta_ref_genome/hg38/Homo_sapiens.GRCh38.dna.primary_assembly.par_y_n_masked.chr2.fasta --input /home/grupo1/resultados_selected/recalibrate/sample1_chr2-171000000-172000000.sorted.marked.base_recalibration.bam --output /home/grupo1/resultados_selected/sample1_chr2-171000000-172000000_first.vcf.gz --intervals /mnt/sdb/curso_bioinformatica/raw_data/targets_TSO/trusight_one_targets_4columns.chr2-171000000-172000000.bed --interval-padding 100
```

## Variant filtering

```
gatk VariantFiltration --variant /home/formador/resultados_selected/variant_calling/sample1_chr2-171000000-172000000_first.vcf.gz --output /home/formador/resultados_selected/variant_calling/sample1_chr2-171000000-172000000_second_filtered.vcf.gz --filter-name QUAL30 --filter-expression "QUAL < 30.0" --genotype-filter-name DP10 --genotype-filter-expression "DP < 10"
```

## Variant selection
```

```

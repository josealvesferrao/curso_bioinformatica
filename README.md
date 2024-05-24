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

Sample1 (replace x by your group number)
```
fastqc -t 1 -o /mnt/sdb/curso_bioinformatica/output/grupox/fastqc /mnt/sdb/curso_bioinformatica/raw_data/sample1/selected_fastq/sample1_chr2-171000000-172000000_R1.fastq.gz /mnt/sdb/curso_bioinformatica/raw_data/sample1/selected_fastq/sample1_chr2-171000000-172000000_R2.fastq.gz
```

Sample2 (replace x by your group number)
```
fastqc -t 1 -o /mnt/sdb/curso_bioinformatica/output/grupox/fastqc /mnt/sdb/curso_bioinformatica/raw_data/sample2/selected_fastq/sample2_chr6-121000000-122000000_R1.fastq.gz /mnt/sdb/curso_bioinformatica/raw_data/sample2/selected_fastq/sample2_chr6-121000000-122000000_R2.fastq.gz
```

Sample3 (replace x by your group number)
```
fastqc -t 1 -o /mnt/sdb/curso_bioinformatica/output/grupo1/fastqc /mnt/sdb/curso_bioinformatica/raw_data/sample3/selected_fastq/sample3_chrX-129000000-131000000_R1.fastq.gz /mnt/sdb/curso_bioinformatica/raw_data/sample3/selected_fastq/sample3_chrX-129000000-131000000_R2.fastq.gz
```

##  Alignment

Sample1 (replace x by your group number)
```
bwa mem -t 1 -R '@RG\tID:sample1\tPU:NB551361.1.sample1\tSM:sample1\tPL:ILLUMINA\tLB:sample1\tCN:curso' /mnt/sdb/curso_bioinformatica/raw_data/fasta_ref_genome/hg38/Homo_sapiens.GRCh38.dna.primary_assembly.par_y_n_masked.chr2.fasta -M /mnt/sdb/curso_bioinformatica/raw_data/sample1/selected_fastq/sample1_chr2-171000000-172000000_R1.fastq.gz /mnt/sdb/curso_bioinformatica/raw_data/sample1/selected_fastq/sample1_chr2-171000000-172000000_R2.fastq.gz -o /mnt/sdb/curso_bioinformatica/output/grupox/sample1_chr2-171000000-172000000.sam
```

## Sort SAM

```
samtools sort –o sample_sorted.sam sample.sam
```

## SAM to BAM

```
samtools view –b sample_sorted.sam > sample_sorted.bam
```

## Index BAM

```
samtools index sample_sorted.bam
```

## MarkDuplicates

```
gatk MarkDuplicates -I sample.bam -O sample_marked.bam -M sample_marked_metrics.txt
```

## BaseRecalibrator

```
gatk BaseRecalibrator -I sample_marked.bam -R genome.fasta --known-sites known_snps.vcf --known-sites known_indels.vcf [ --intervals target_positions.bed --interval-padding 100 ] -O sample_marked_baserecalibrator_report.txt
```

## ApplyBaseRecalibrator

```
gatk ApplyBQSR -I sample_marked.bam -R genome.fasta --bqsr-recal-file sample_marked_baserecalibrator_report.txt -O sample_marked_baserecalibrator.bam
```

## Qualimap

```
qualimap bamqc -bam sample_marked_recalibrated.bam -outdir . -outfile sample_marked_recalibrated_qualimap.pdf -outformat PDF [ -gff targets_position.bed ]
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
gatk SelectVariants --variant /home/formador/resultados_selected/variant_calling/sample1_chr2-171000000-172000000_second_filtered.vcf.gz --output /home/formador/resultados_selected/variant_calling/sample1_chr2-171000000-172000000_second_filtered_selected.vcf.gz --exclude-filtered true --intervals /mnt/sdb/curso_bioinformatica/raw_data/targets_TSO/trusight_one_targets_4columns.chr2-171000000-172000000.bed
```

## Link VEP Ensembl Web interface https://www.ensembl.org/Tools/VEP

## VEP annotation - VCF output (conda activate curso_amb_vep)

```
vep -i /mnt/sdb/curso_bioinformatica/raw_data/sample1/full_vcf/sample1_single_call_standard_filtered_selected_targets.final.vcf.gz --cache --dir_cache /mnt/sdb/curso_bioinformatica/vep --dir_plugins /mnt/sdb/curso_bioinformatica/vep/Plugins --species homo_sapiens --assembly GRCh38 --format vcf --vcf --output_file /mnt/sdb/curso_bioinformatica/output/formador/sample1_single_call_standard_filtered_selected_targets.final_annotated.vcf.gz --force_overwrite --stats_file /mnt/sdb/curso_bioinformatica/output/formador/sample1_single_call_standard_filtered_selected_targets.final_annotated.html --fork 2 --fasta /mnt/sdb/curso_bioinformatica/raw_data/fasta_ref_genome/hg38/Homo_sapiens.GRCh38.dna.primary_assembly.par_y_n_masked.fa --everything --merged --nearest symbol --allele_number --show_ref_allele --hgvsg --transcript_version --mane_select --check_existing --exclude_null_alleles --clin_sig_allele 1 --exclude_predicted --pick_allele_gene
```

## VEP annotation - TAB output (conda activate curso_amb_vep)

```
vep -i /mnt/sdb/curso_bioinformatica/raw_data/sample1/full_vcf/sample1_single_call_standard_filtered_selected_targets.final.vcf.gz --cache --dir_cache /mnt/sdb/curso_bioinformatica/vep --dir_plugins /mnt/sdb/curso_bioinformatica/vep/Plugins --species homo_sapiens --assembly GRCh38 --format vcf --tab --output_file /mnt/sdb/curso_bioinformatica/output/formador/sample1_single_call_standard_filtered_selected_targets.final_annotated.tab --force_overwrite --stats_file /mnt/sdb/curso_bioinformatica/output/formador/sample1_single_call_standard_filtered_selected_targets.final_annotated_new.html --fork 2 --fasta /mnt/sdb/curso_bioinformatica/raw_data/fasta_ref_genome/hg38/Homo_sapiens.GRCh38.dna.primary_assembly.par_y_n_masked.fa --everything --merged --nearest symbol --allele_number --show_ref_allele --hgvsg --transcript_version --mane_select --check_existing --exclude_null_alleles --clin_sig_allele 1 --exclude_predicted --pick_allele_gene
```

## Link Exomiser Web interface https://exomiser.monarchinitiative.org/exomiser/

## Exomiser (conda activate base)

```
java -Xmx4g -jar /mnt/sdb/curso_bioinformatica/exomiser/exomiser-cli-13.1.0/exomiser-cli-13.1.0.jar --assembly hg38 --analysis /mnt/sdb/curso_bioinformatica/raw_data/sample1/exomiser_files/sample1_full_clinical_exome_analysis_file.yml --output-prefix /mnt/sdb/curso_bioinformatica/output/formador/exomiser_output/sample1_Exomiser --spring.config.location=/mnt/sdb/curso_bioinformatica/exomiser/exomiser-cli-13.1.0/application.properties
```

## [Optional] Virtual gene panel / Filter by HPO 150 most associated genes and AF cutoff. Last argument is AF (eg. 0.01 is 1%), try other AF values.

```
python /mnt/sdb/curso_bioinformatica/raw_data/scripts/Script_HPOtoGENE_vs_originalVEP-VCF_genome19_JAF_17-10-2022.py sample1 /mnt/sdb/curso_bioinformatica/raw_data/sample1/vep_files/sample1_vep_annotation_no_comments.tab /mnt/sdb/curso_bioinformatica/raw_data/sample1/phen2gene_list/sample1_phen2gene_list.txt /mnt/sdb/curso_bioinformatica/raw_data/scripts/gene_list_trusight_one.txt 0.01
```

## Useful links
https://hpo.jax.org/app/ - Human Phenotype Ontology (HPO)

https://bioportal.bioontology.org/ontologies/HP/?p=classes&conceptid=root - BioPortal has many ontologies, including HPO

https://doc2hpo.wglab.org/ - To automatically convert clinical description/phenotype to HPO terms

https://phen2gene.wglab.org/ - Get a list of genes associated with a set of HPO terms

https://www.internationalgenome.org/data-portal/sample - Here you may find and download fastq files of exome sequencing. Choose one sample and select "sequence" and "exome". They also have alignment (SAM/BAM/CRAM) and variant (VCF) files

https://www.ensembl.org/index.html - You may look for information regarding a certain gene in Ensembl. Choose Human and search for a gene, press "Go"

https://hsf.genomnis.com/login - You can use HSF for splicing analysis (predict the impact of variants)

https://www.ncbi.nlm.nih.gov/snp/

https://www.ncbi.nlm.nih.gov/clinvar/

https://mastermind.genomenon.com/ - You may use Mastermind to search for a specific variant on the literature 

https://www.uniprot.org/

https://gtexportal.org/home/, https://www.bgee.org/ - You may check gene expression levels (in different species or tissues)

https://genome-euro.ucsc.edu/cgi-bin/hgGateway?redirect=manual&source=genome.ucsc.edu - You may search for a gene or specific coordinate, and select a set of tracks to see a lot of information regarding that genomic region

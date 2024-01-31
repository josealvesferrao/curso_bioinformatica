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



# ifpan-nalepa-smog

1. QC control of Fastq files:
``` 
RNA_FILES=`ls */*/*/*gz`
echo $RNA_FILES | tr -s ' \011' '\012' | head -n 20 | tail -n 10 | xargs -i docker run -d --rm -v $PWD:/data pegi3s/fastqc /data/{}
```

2. MultiQC to generate a full report:


3. Mus musculus reference genomes was downloaded from Ensmbl (release 102), and indexed with Hisat2:

4. Alignement for a single sample:
```
ls */*/*/*gz | cut -d "/" -f 4 | cut -d "_" -f 1 | uniq >> sample-list.txt

cat sample-list.txt | head -2 | xargs -I {} docker run --rm -v $PWD:/data zlskidmore/hisat2:2.1.0 hisat2 -x /data/index/genome -1 /data/X201SC20111890-Z01-F001/raw_data/{}/{}_1.fq.gz -2 /data/X201SC20111890-Z01-F001/raw_data/{}/{}_2.fq.gz -S /data/{}.sam --summary-file /data/{}.txt --dta-cufflinks

```

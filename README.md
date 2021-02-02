# ifpan-nalepa-smog

1. QC control of Fastq files:
``` 
RNA_FILES=`ls */*/*/*gz`
echo $RNA_FILES | tr -s ' \011' '\012' | head -n 20 | tail -n 10 | xargs -i docker run -d --rm -v $PWD:/data pegi3s/fastqc /data/{}
```

2. MultiQC to generate a full report:


3. Mus musculus reference genomes was downloaded from Ensmbl (release 102), and indexed with Hisat2:

4. Alignement for a single sample:


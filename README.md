# ifpan-nalepa-smog

1. QC control of Fastq files:
``` 
RNA_FILES=`ls */*/*/*gz`
echo $RNA_FILES | tr -s ' \011' '\012' | head -n 20 | tail -n 10 | xargs -i docker run -d --rm -v $PWD:/data pegi3s/fastqc /data/{}
```

2. MultiQC to generate a full report:

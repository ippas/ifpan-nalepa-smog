# ifpan-nalepa-smog

1. QC control of Fastq files:
``` 
RNA_FILES=`ls */*/*/*gz`
echo $RNA_FILES | tr -s ' \011' '\012' | head -n 20 | tail -n 10 | xargs -i docker run -d --rm -v $PWD:/data pegi3s/fastqc /data/{}
```

2. MultiQC to generate a full report:
```
docker run --rm -v $PWD:/data ewels/multiqc:latest multiqc /data -o /data
```

3. Mus musculus reference genomes was downloaded from Ensmbl (release 102), and indexed with Hisat2:

4. Alignement for a single sample:
```
ls */*/*/*gz | cut -d "/" -f 4 | cut -d "_" -f 1 | uniq >> sample-list.txt

cat sample-list.txt | head -2 | xargs -I {} docker run --rm -v $PWD:/data zlskidmore/hisat2:2.1.0 hisat2 -x /data/index/genome -1 /data/X201SC20111890-Z01-F001/raw_data/{}/{}_1.fq.gz -2 /data/X201SC20111890-Z01-F001/raw_data/{}/{}_2.fq.gz -S /data/{}.sam --summary-file /data/{}.txt --dta-cufflinks

```
5. Convert sam to bam
```
cat sample-list.txt | head -6 | xargs -I {} docker run -d --rm -v $PWD:/data zlskidmore/samtools:1.9 samtools sort -O bam -T /data/$FILE.sort -o /data/{}.bam /data/{}.sam
```
6. Run cuffcuant:

```
cat sample-list.txt | head -6 | xargs -I {} docker run -d --rm -v $PWD:/data octavianus90/cufflinks_final:latest cuffquant -o /data/{}.cuant /data/Mus_musculus.GRCm38.102.gtf /data/{}.bam
```


### Fragments, mateusz analysis

Preparation of marker genes:

The marker genes were [taken from](http://mousebrain.org/celltypes/?fbclid=IwAR2uLbp0fYm2Eaet7l_vz9OYeoTIV_qByP6eEddBvwIx6-55GKGnHu5TaiQ), and extract need information:

```
cat gene-markes | 
    cut -f3,5 | 
    tail +2  | 
    grep -P 'Excitatory|Inhibitory|Microglia|astrocytes' | 
    awk -F"\t" 'BEGIN {OFS="\t"} {gsub(" ", "-", $1); print}'| 
    sed 's/ /\t/g' |  
    awk -F"\t" -v OFS="\t" '{print $1, $2"\n"$1, $3"\n"$1, $4"\n"$1, $5"\n"$1, $6"\n"$1, $7"\n"$1, $8}' | 
    awk '$2!=""' 
    
    cat gene-markes | cut -f3,5,7 | tail +2 | grep -P 'Excitatory|Inhibitory|Microglia|astrocytes' | awk -F"\t" -v OFS="\t" '{print $1, $3, $2}' | awk -F"\t" 'BEGIN {OFS="\t"} {gsub(" ", ";", $1); print}' | awk -F"\t" 'BEGIN {OFS="\t"} {gsub(" ", ";", $2); print}' | sed 's/ /\t/g' | awk -F"\t" -v OFS="\t" '{print $1, $2, $3"\n"$1, $2, $4"\n"$1, $2, $5"\n"$1, $2, $6"\n"$1, $2,  $7"\n"$1, $2, $8"\n"$1, $2,  $9}' | awk '$3!=""' | grep -P "Cortex|CNS|Hypothalamus" | sed 's/;/ /g'

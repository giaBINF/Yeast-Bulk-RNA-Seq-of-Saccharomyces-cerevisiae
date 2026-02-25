
## 0. Downloading SSRs and References
*SRA tools tutorials: [prefetch](https://github.com/ncbi/sra-tools/wiki/08.-prefetch-and-fasterq-dump) & [fasterq-dump](https://github.com/ncbi/sra-tools/wiki/HowTo:-fasterq-dump)*
```
# Extracting FASTQ files from SRA-accessions
for srr in SRR10551665 SRR10551664 SRR10551663 \
           SRR10551662 SRR10551661 SRR10551660 \
           SRR10551659 SRR10551658 SRR10551657
do
    prefetch $srr
    fasterq-dump $srr --split-files -e 4
    gzip ${srr}_*.fastq
done

# Download reference transcriptome 
curl -L \
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.gtf.gz \
-o data/references/S_cerevisiae_R64_annotation.gtf.gz

# Download reference annotation 
curl -L \
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.gtf.gz \
-o data/references/S_cerevisiae_R64_annotation.gtf.gz
```

## 1. Quality Control 
```
# Running FastQC using Docker
docker run --rm \
  -v "$PWD":/work \
  -w /work \
  quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0 \
  fastqc data/raw/*_1.fastq.gz -o results/fastqc -t 4

# Running MultiQC using Docker
docker run --rm \
  -v "$PWD":/work \
  -w /work \
  quay.io/biocontainers/multiqc:1.21--pyhdfd78af_0 \
  multiqc results/fastqc -o results
```

## 2. Quantification with Salmon
*Indexing and Quantification based on [this tutorial](https://www.hadriengourle.com/tutorials/rna/)*
```
# Indexing transcriptome 
docker run --rm \
  -v "$PWD":/work \
  -w /work \
  combinelab/salmon:latest \
  salmon index \
    -t data/references/S_cerevisiae_R64_transcripts.fa \
    -i data/references/salmon_index \
    -p 8

# Quantifying reads suing Salmon
for i in data/raw/*_1.fastq.gz
do
   prefix=$(basename "$i" _1.fastq.gz)

   echo "Processing ${prefix}"

   docker run --rm \
     -v "$PWD":/work \
     -w /work \
     combinelab/salmon:latest \
     salmon quant \
        -i data/references/salmon_index \
        --libType A \
        -r "$i" \
        --validateMappings \
        -p 8 \
        -o results/salmon/${prefix}_quant

done

# Creating sample table for tximport in R
find results/salmon -name "quant.sf" \
| awk -F'/' '{s=$(NF-1); sub(/_quant$/, "", s); print s "\t" $0}' \
> results/salmon_import.tsv  
```



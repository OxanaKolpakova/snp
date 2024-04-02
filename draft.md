# SNP
Try to make nextflow for variant calling

```
java -jar picard.jar CreateSequenceDictionary R=/home/alexandr/Documents/SNP/data/MT_short.fna O=/home/alexandr/Documents/SNP/data/MT_short.dict

    gatk HaplotypeCaller \
         -R /home/alexandr/Documents/SNP/data/MT_short.fna \
         -I /home/alexandr/Documents/SNP/results/bwamem/sample.bam \
         -O /home/alexandr/Documents/SNP/results/sample.vcf
```

```
mkdir -p $GOPATH/src/github.com/sylabs && \
    cd $GOPATH/src/github.com/sylabs && \
    wget https://github.com/sylabs/singularity/releases/download/v4.1.2/singularity-ce-4.1.2.tar.gz && \
    tar -xzf singularity-${VERSION}.tar.gz && \
    cd ./singularity && \
    ./mconfig
```

```
singularity exec library://funsinaime/gatk/withlc2 gatk HaplotypeCaller \
    -R /path/to/your/reference.fasta \
    -I /path/to/your/input.bam \
    -O output.vcf
```

```
samtools sort SRR14298060.bam -o SRR14298060.sorted.bam
samtools index ${bamFile.baseName}.sorted.bam
samtools mpileup -uf $reference ${bamFile.baseName}.sorted.bam > ${bamFile.baseName}.pileup
```
# SNP Pipeline

SNP Pipeline is a Nextflow-based pipeline for Single Nucleotide Polymorphism (SNP) calling from genome sequencing data. This pipeline provides streamlined processing from sequencing to annotation, assisting researchers in analyzing genetic variations with minimal effort.

## Features

- Alignment using BWA mem for accurate mapping of sequences to the reference genome.
- Quality control using fastp to remove adapters and low-quality reads.
- Variant calling using mpileup for SNP and indel detection.
- Variant annotation using Variant Effect Predictor (VEP) to identify functional consequences of variations.
- Support for containerization using Singularity, ensuring ease of deployment and management of the execution environment.

## DAG (Directed Acyclic Graph)


## Requirements

- Installed Nextflow (https://www.nextflow.io/docs/latest/getstarted.html)
- Installed Singularity (https://docs.sylabs.io/guides/3.0/user-guide/installation.html)
- Resources for analysis: computational power, access to sequencing data.
- Singularity for containerization. Make sure Singularity is installed and accessible in your environment.

## Installation and Execution

1. Install Nextflow by following the instructions on the official website. (https://www.nextflow.io/docs/latest/getstarted.html)

2. Clone the SNP Pipeline repository to your local machine.
```
git clone https://github.com/glebus-sasha/SNP.git
```
3. Navigate to the directory containing the pipeline and execute it using the command `nextflow run SNP.nf`.

## Usage ![DAG](Examplefile:///home/alexandr/%D0%98%D0%B7%D0%BE%D0%B1%D1%80%D0%B0%D0%B6%D0%B5%D0%BD%D0%B8%D1%8F/%D0%A1%D0%BD%D0%B8%D0%BC%D0%BA%D0%B8%20%D1%8D%D0%BA%D1%80%D0%B0%D0%BD%D0%B0/%D0%A1%D0%BD%D0%B8%D0%BC%D0%BE%D0%BA%20%D1%8D%D0%BA%D1%80%D0%B0%D0%BD%D0%B0%20%D0%BE%D1%82%202024-04-02%2011-25-22.png)


Example command to run the pipeline:

```
nextflow run SNP.nf
```

## Contributors

- Glebus Aleksandr ([@glebus-sasha](https://github.com/glebus-sasha/))
- Oxana Kolpakova ([@OxanaKolpakova](https://github.com/OxanaKolpakova))

## License

SNP Pipeline is distributed under the MIT license. See the LICENSE file for details.

# SNP Pipeline

SNP Pipeline is a Nextflow-based pipeline for Single Nucleotide Polymorphism (SNP) calling from genome sequencing data. This pipeline provides streamlined processing from sequencing to annotation, assisting researchers in analyzing genetic variations with minimal effort.

## Features

- Alignment using BWA mem for accurate mapping of sequences to the reference genome.
- Quality control using fastp to remove adapters and low-quality reads.
- Variant calling using mpileup for SNP and indel detection.
- Variant annotation using Variant Effect Predictor (VEP) to identify functional consequences of variations.
- Support for containerization using Singularity, ensuring ease of deployment and management of the execution environment.

## DAG (Directed Acyclic Graph)

```html

<!--
  ~ Copyright 2013-2023, Seqera Labs
  ~
  ~ Licensed under the Apache License, Version 2.0 (the "License");
  ~ you may not use this file except in compliance with the License.
  ~ You may obtain a copy of the License at
  ~
  ~     http://www.apache.org/licenses/LICENSE-2.0
  ~
  ~ Unless required by applicable law or agreed to in writing, software
  ~ distributed under the License is distributed on an "AS IS" BASIS,
  ~ WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  ~ See the License for the specific language governing permissions and
  ~ limitations under the License.
  -->
<html>
<head>
<meta name="viewport" content="width=device-width, user-scalable=no, initial-scale=1, maximum-scale=1">
</head>
<body>
<pre class="mermaid" style="text-align: center;">
flowchart TB
    subgraph " "
    v0["Channel.fromFilePairs"]
    v1["reference"]
    v5["reference"]
    v8["reference"]
    end
    v2([REFINDEX])
    v3([QCONTROL])
    subgraph " "
    v4[" "]
    v11[" "]
    end
    v6([ALIGN])
    v7([PREPARE])
    v9([CALLSNP])
    v10([ANNOTATE])
    v0 --> v3
    v1 --> v2
    v2 --> v6
    v3 --> v6
    v3 --> v4
    v5 --> v6
    v6 --> v7
    v7 --> v9
    v8 --> v9
    v9 --> v10
    v10 --> v11

</pre>
<script type="module">
  import mermaid from 'https://cdn.jsdelivr.net/npm/mermaid@10/dist/mermaid.esm.min.mjs';
  mermaid.initialize({ startOnLoad: true });
</script>
</body>
</html>
```

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

## Usage Example

Example command to run the pipeline:

```
nextflow run SNP.nf
```

## Contributors

- Glebus Aleksandr ([@glebus-sasha](https://github.com/glebus-sasha/))
- Oxana Kolpakova ([@OxanaKolpakova](https://github.com/OxanaKolpakova))

## License

SNP Pipeline is distributed under the MIT license. See the LICENSE file for details.
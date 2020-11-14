# reprohackathon-bioinfo

The goal is to reproduce parts of the analysis described in these papers :<br>
  - https://pubmed.ncbi.nlm.nih.gov/23313955/<br>
  - https://pubmed.ncbi.nlm.nih.gov/23861464/<br>

They performed RNA-Seq in samples from patients with uveal melanoma. Some samples are mutated in SF3B1. <br>

Data from the studies are avalaible at : http://www.ncbi.nlm.nih.gov/sra?term=SRA062359

At least 8 cpus and 32GB of RAM are required to run the pipeline.

To run the pipeline on IFB Cloub on Appliance BioPipes :
```bash
conda activate
nextflow run main.nf -c nextflow.config -resume {-N emailAddress(for email report) -bg (for background execution) }
```

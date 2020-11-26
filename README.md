# reprohackathon-bioinfo

The goal is to reproduce parts of the analysis described in these papers :<br>
  - https://pubmed.ncbi.nlm.nih.gov/23313955/<br>
  - https://pubmed.ncbi.nlm.nih.gov/23861464/<br>
  
download at the site down below on "supplementary tables" to learn about the number of genes identified as presenting differential exon usage <br> 
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5321577/?fbclid=IwAR0N5sBxvbJ4imv3H0caQ9noRJc0qNsdUx0lCBThdPnL21DhyGFfpVwlq3s


They performed RNA-Seq in samples from patients with uveal melanoma. Some samples are mutated in SF3B1. <br>

Data from the studies are avalaible at : http://www.ncbi.nlm.nih.gov/sra?term=SRA062359

At least 8 cpus and 40GB of RAM are required to run the pipeline.

To run the pipeline on IFB Cloub on Appliance BioPipes :
```bash
conda activate
nextflow run main.nf -c nextflow.config -resume {-N emailAddress(for email report) -bg (for background execution) }
```

Results of the pipeline can be found in /files/results/results.txt

# canine-rnaseq-benchmark
Benchmarking DESeq2, edgeR, and limma-voom using simulated and canine lymphoma RNA-seq data

## Project goal
This project benchmarks DESeq2, edgeR, and limma-voom for identifying differentially expressed genes from RNA-seq count data. The methods will first be evaluated using simulated data where the truly differentially expressed genes are known, and then applied to a real canine lymphoma dataset. The project connects statistical method evaluation with biologically relevant questions about canine lymphoma.

## Tools
This project compares three R/Bioconductor methods for RNA-seq differential-expression analysis:

- DESeq2
- edgeR
- limma-voom
  
## Data
The real-data analysis will use GSE196186, a canine whole-blood RNA-seq count dataset containing seven lymphoma samples and seven non-lymphoma samples. Simulated datasets will be generated with known differential-expression truth for evaluating accuracy.

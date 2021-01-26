# amplicon_cov

This script computes per amplicon relative coverage for a batch of samples. It outputs a .csv with samples in the rows, first column is sample name and 98 other columns are amplicons.  
It needs a bedfile of primers, a .tsv to list the samples (like `/cluster/project/pangolin/working/samples.tsv`) and a directory where to find the aligned samples (like `/cluster/project/pangolin/working/samples/`).  
Optionally outputs plots.

## Example of usage:

If we want to analyze coverage of samples in 'samples20210122_HY53JDRXX.tsv', with primers info in 'articV3primers.bed', while plotting and with verbose, outputting into '20210122out' (I think you need to create directory first): 

```python ./amplicon_covs.py -pv -s samples20210122_HY53JDRXX.tsv -r articV3primers.bed -o 20210122out/```

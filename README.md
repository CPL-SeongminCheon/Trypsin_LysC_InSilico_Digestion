# Trypsin_LysC_InSilico_Digestion
In Silico Trypsin/Lys-C Clevage Site Prediction

---

### Requirement
Python 3.7
Install Biopython


### Installation requirement on Conda environment
```
conda install -c conda-forge biopython
```
[Biopython Install guide with Anaconda](https://anaconda.org/conda-forge/biopython)
  


### Usage
```
python Trypsin-lysC_clevage_predict.py -i <input fasta file> -o <input fasta file> -miss <Allow missclevage> -min <Minimum length of digested peptide>
```
---

### Example usage
```
python Trypsin-lysC_clevage_predict.py -i example.fa -o output.fa -miss 2 -min 6
```

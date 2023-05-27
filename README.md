# EDWARD (CSE 185 Project)

A data visualization tool that can reveal high-level relationships among classes in a dataset by performing PCA and t-SNE dimensionality reduction given an input VCF file of SNPs or an RNA-sequencing count matrix.

# Installation Instruction

Package requires `numpy`, `cyvcf2`, `pandas` all specified in requirements.txt. 
```
pip install -r requirements.txt
```

Once required libraries are installed, you can install mypileup with the following command:

```
python setup.py install
```

for testing run './untitled.py -t v -i [input_file] --pca --num_PCs'
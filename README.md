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
*note: this package can only be installed on linux (including MacOS) devices 

# Basic Usage

The basic usage of `edward` is
```
edward [required] --type <v or c> --pca/umap/tsne --input <input> [options] --prefix --leiden --louivain --number-of-pcs 
```

To run `edward` on a small test example (using files in this repo):

```
edward -t v -i example_data/sample.vcf.gz --pca -n 2
```
* For now `output.html` will be generated

# Documentation 

```
  INPUTS AND OUTPUTS:
    --type -t: 'v' for vcf, 'c' for RNA count matrix;
    --prefix -p: prefix for tsv and html outputs [default=None]
    DIMENSIONALITY TOOLS:
    --pca: run pca algorithm
    --umap: run umap algorithm
    --tsne: run tsne algorithm 
    CLUSTERING (optional):
    --leiden: [default=None]
    --louvain: [default=None]
    EXTRA OPTIONS: 
    --number-of-pcs -n: number of pcs for PCA [default=5]
```

# Contributors

This project is for Dr. Melissa Gymreks CSE 185 class with instructional staff Ryan Eveloff, Luisa Amaral, and Himanshu Nln. <p>

The authors are Minyoung Ahn, Joelle Faybishenko, and Minh-Son Tran
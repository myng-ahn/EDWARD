import argparse
import sys, getopt
from cyvcf2 import VCF
import numpy as np 
'''
pca input
- data: numpy ndarray
- ncomponents: int
'''

def main():
    pass

def process_pca(vcf_path):
    for variant in VCF('some.vcf.gz'):

        variant.gt_types # numpy array
        variant.gt_ref_depths, variant.gt_alt_depths # numpy arrays
        variant.gt_phases, variant.gt_quals # numpy arrays
        variant.gt_bases # numpy array
        variant.CHROM, variant.start, variant.end, variant.ID, \
        variant.REF, variant.ALT, variant.FILTER, variant.QUAL
        variant.INFO.get('DP') # int
        variant.INFO.get('FS') # float
        variant.INFO.get('AC') # float
        a = variant.gt_phred_ll_homref # numpy array
        b = variant.gt_phred_ll_het # numpy array
        c = variant.gt_phred_ll_homalt # numpy array

        str(variant)

if __name__ == "__main__":
    
    USAGE_MSG = "USAGE: python edward.py -h [required] --type <v or c> --pca/umap/tsne <input> [options] --prefix --leiden --louivain --number-of-pcs "

    HELP='''
    USAGE: python edward.py -h [required] --type <v or c> --pca/umap/tsne <input> [options] --prefix --leiden --louivain --number-of-pcs
    
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
    '''



    # Create the parser
    parser = argparse.ArgumentParser(description='Description of your program')

    # Define the arguments
    parser.add_argument('-t', '--type', required=True, choices=['v', 'c'], help="'v' for vcf, 'c' for RNA count matrix")
    parser.add_argument('-i', '--input', required=True, help="input file as specified by type")
    parser.add_argument('-p', '--prefix', default=None, help="prefix for tsv and html outputs")
    parser.add_argument('--pca', required=True, action='store_true', help="run pca algorithm")
    parser.add_argument('--umap', action='store_true', help="run umap algorithm")
    parser.add_argument('--tsne', action='store_true', help="run tsne algorithm")
    parser.add_argument('--leiden', default=None, help="leiden argument")
    parser.add_argument('--louvain', default=None, help="louvain argument")
    parser.add_argument('-n', '--number-of-pcs', default=5, type=int, help="number of pcs for PCA")

    # Parse the arguments
    args = parser.parse_args()

    # Access the parsed arguments
    type_arg = args.type
    input_arg = args.input
    prefix_arg = args.prefix
    pca_arg = args.pca
    umap_arg = args.umap
    tsne_arg = args.tsne
    leiden_arg = args.leiden
    louvain_arg = args.louvain
    number_of_pcs_arg = args.number_of_pcs

    # Print the parsed arguments
    print(f'Type: {type_arg}')
    print(f'Prefix: {prefix_arg}')
    print(f'PCA: {pca_arg}')
    print(f'UMAP: {umap_arg}')
    print(f'TSNE: {tsne_arg}')
    print(f'Leiden: {leiden_arg}')
    print(f'Louvain: {louvain_arg}')
    print(f'Number of PCs: {number_of_pcs_arg}')


    

    pca_matrix = process_pca(input_arg)
    

    

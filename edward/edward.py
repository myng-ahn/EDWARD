#!/usr/bin/env python

import argparse
import sys
import numpy as np
import os
from cyvcf2 import VCF
import pca as pca
from html import *
#from edward import __version__


'''
pca input
- data: numpy ndarray
- ncomponents: int
'''

def main():
    pass

# for testing run './untitled.py -t v -i [input_file] --pca --num_PCs'
def process_pca(vcf_path, number_of_pcs_arg):
    '''
    TODO: documentation
    '''
    gt_array = []
    for variant in VCF(vcf_path):
        gt = [sum(genotypes[:-1]) for genotypes in variant.genotypes]
        # variant.genotypes is list containing [0, 0, True], [1, 1, True], [0, 1, False], ... 
        # where the first 2 elements are the genotype and the third element is a boolean indicating
        # whether that genotype is phased or not (for now I'm ignoring this but it may be important
        # to consider later. 
        # 
        # NOTE: some of the genotypes only have one allele and some have negative alleles 
        # whatever that even means. For now I'm just gonna keep going but something to investigate.
        gt_array.append(gt)
    gt_array = np.array(gt_array) # convert to numpy ndarray
    pca_transformed = pca(gt_array, number_of_pcs_arg) # perform pca
    
    print(pca_transformed) # for testing. should comment out
    # TODO: plot pca_transformed and output to html

def main():
    
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
    parser = argparse.ArgumentParser(
        prog='edward',
        description='exceptional-delicate-wonderful-attempt(at)-reducing-dimensions')

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
    #parser.add_argument("--version", help="Print the version and quit", \
	#	action="version", version = '{version}'.format(version=__version__))

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

    
    
    pca_matrix = process_pca(input_arg, number_of_pcs_arg)
    pca_output = pca(pca_matrix, number_of_pcs_arg)
    ## todo: add min-sohn html function call here

    sys.exit(0)


if __name__ == "__main__":
    main()
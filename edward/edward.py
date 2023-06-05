#!/usr/bin/env python

import argparse
import sys
import numpy as np
import pandas as pd 
import os
from cyvcf2 import VCF
#from pca import pca
#from html import *
#from write_html import write_html
import matplotlib.pyplot as plt
from edward import __version__
from . import pca as pca
from . import write_html as write_html


'''
pca input
- data: numpy ndarray
- ncomponents: int
'''

def main():
    pass

# for testing run './untitled.py -t v -i [input_file] --pca --num_PCs'
def process_vcf(vcf_path):
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
    return gt_array
    pca_proj, sorted_eigvals, sorted_eigvecs = pca.pca(gt_array, number_of_pcs_arg) # perform pca
    return pca_proj, sorted_eigvals, sorted_eigvecs
    print(pca_transformed) # for testing. should comment out
    # TODO: plot pca_transformed and output to html

def process_count(matrix_path):
    '''
    TODO: documentation
    '''
    # rows = samples, columns = genes
    gt_array = []
    gt = pd.read_csv(matrix_path)  
    gt_array = np.array(gt.values) # convert to numpy ndarray
    return gt_array
    pca_proj, sorted_eigvals, sorted_eigvecs = pca.pca(gt_array, number_of_pcs_arg) # perform pca
    return pca_proj, sorted_eigvals, sorted_eigvecs
    print(pca_transformed) # for testing. should comment out


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
    parser.add_argument('--pca', action='store_true', help="run pca algorithm")
    parser.add_argument('--umap', action='store_true', help="run umap algorithm")
    parser.add_argument('--tsne', action='store_true', help="run tsne algorithm")
    parser.add_argument('--leiden', action='store_true', help="leiden argument")
    parser.add_argument('--louvain', action='store_true', help="louvain argument")
    parser.add_argument('-n', '--number-pcs', default=5, type=int, help="number of pcs for PCA")
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
    number_of_pcs_arg = args.number_pcs

    # Print the parsed arguments
    print(f'Type: {type_arg}')
    print(f'Prefix: {prefix_arg}')
    print(f'PCA: {pca_arg}')
    print(f'UMAP: {umap_arg}')
    print(f'TSNE: {tsne_arg}')
    print(f'Leiden: {leiden_arg}')
    print(f'Louvain: {louvain_arg}')
    print(f'Number of PCs: {number_of_pcs_arg}')

    
    ### 
    # check out arguments
    #
    #
    ###

    if not (pca_arg or umap_arg or tsne_arg): 
        print("PCA, UMAP, or TSNE not declared")
        sys.exit(1)

    if type_arg =='v': 
        array = process_vcf(input_arg)
    else:
        array = process_count(input_arg)

    if pca_arg:
        pca_output, pca_eigvals, pca_eigvecs = pca.pca(array, number_of_pcs_arg) # perform pca
        # make a relevant figure??? 
        fig = plt.figure()
        write_html.write_html(input_arg, -1, -1, pca_figs=[fig], pca_eigvals=pca_eigvals, pca_eigvecs=pca_eigvecs)
    
    ##testing dummy figures
    #fig = plt.figure()
    #write_html.write_html(input_arg, -1, -1, pca_figs=[fig], pca_eigvals=pca_eigvals, pca_eigvecs=pca_eigvecs)

    sys.exit(0)


if __name__ == "__main__":
    main()
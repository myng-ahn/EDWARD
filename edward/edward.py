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
import matplotlib.colors as mcolors
from edward import __version__
from . import dim_reduce as dim_reduce
from . import write_html as write_html
from . import cluster as cluster

'''
pca input
- data: numpy ndarray
- ncomponents: int
'''

COLORS = list(mcolors.XKCD_COLORS.values())

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
    row, col = gt_array.shape
    return gt_array, row, col
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
    row, col = gt_array.shape
    return gt_array, row, col
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
    parser.add_argument('-p', '--prefix', default='', help="prefix for tsv and html outputs")
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
        array, samples, observations = process_vcf(input_arg)
    else:
        array, samples, observations = process_count(input_arg)


 
    pca_figs = None
    pca_eigvals = None
    pca_eigvecs = None
    tsne_fig = None 
    umap_fig = None

    if pca_arg:
        pca_output, pca_eigvals, pca_eigvecs = dim_reduce.pca(array, number_of_pcs_arg) # perform pca
        k = np.sqrt(pca_output.shape[0])
        c = ['#52b2BF'] * pca_output.shape[0]
        if leiden_arg:
            idents = cluster.leiden(pca_output, k)
            c = [COLORS[idx] for idx in idents]
        elif louvain_arg:
            idents = cluster.louvain(pca_output, k)
            c = [COLORS[idx] for idx in idents]

        pca_figs=[]
        for i in range(number_of_pcs_arg):
            for j in range(number_of_pcs_arg):
                if j <= i:
                    continue
                fig = plt.figure()
                ax = fig.add_subplot(111)
                ax.scatter(pca_output[:, i], pca_output[:, j], c=c)
                ax.set_xlabel('PC{}'.format(i+1))
                ax.set_ylabel('PC{}'.format(j+1))
                pca_figs.append(fig)
        
        path = prefix_arg + '_pca_eigenval.txt'
        with open(path, "w") as file:
            for val in pca_eigvals:
                file.write(str(val) + "\n")
        path = prefix_arg + '_pca_eigenvec.txt'
        with open(path, "w") as file:
            for row in pca_eigvecs:
                row_str = " ".join(str(element) for element in row)
                file.write(row_str + "\n")
        
    if umap_arg:
        umap_transformed = dim_reduce.umap(array)
        k = np.sqrt(umap_transformed.shape[0])
        c = ['#52b2BF'] * umap_transformed.shape[0]
        if leiden_arg:
            idents = cluster.leiden(umap_transformed, k)
            c = [COLORS[idx] for idx in idents]
        elif louvain_arg:
            idents = cluster.louvain(umap_transformed, k)
            c = [COLORS[idx] for idx in idents]

        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.scatter(umap_transformed[:, 0], umap_transformed[:, 1], c=c)
        ax.set_xlabel('UMAP_1')
        ax.set_ylabel('UMAP_2')
        umap_fig = fig

    if tsne_arg:
        perplexity = np.sqrt(array.shape[0])
        tsne_transformed = dim_reduce.tsne(array, perplexity)
        k = np.sqrt(tsne_transformed.shape[0])
        c = ['#52b2BF'] * tsne_transformed.shape[0]
        if leiden_arg:
            idents = cluster.leiden(tsne_transformed, k)
            c = [COLORS[idx] for idx in idents]
        elif louvain_arg:
            idents = cluster.louvain(tsne_transformed, k)
            c = [COLORS[idx] for idx in idents]

        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.scatter(tsne_transformed[:, 0], tsne_transformed[:, 1], c=c)
        ax.set_xlabel('TSNE_1')
        ax.set_ylabel('TSNE_2')
        tsne_fig = fig


    write_html.write_html(input_arg, 
                          samples, observations, 
                          pca_figs=pca_figs, pca_eigvals=pca_eigvals, pca_eigvecs=pca_eigvecs, 
                          umap_fig=umap_fig,
                          tsne_fig=tsne_fig,
                          output=prefix_arg)


    sys.exit(0)


if __name__ == "__main__":
    main()
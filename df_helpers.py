import pandas as pd
import numpy as np
import os

dirs_list = ['1. Batch Counts', '2. DESeq2 Input', '3. DESeq2 Output', '4. Significant DEGs', 
        '5. PCA Plots', '6. GSEA Plots', '6. msigdb gene sets', '7. Heatmaps','8. Volcano Plots','9. TF Analysis']
batchcounts, deseq2input, deseq2output, sigdegs, pcaplots, gsea, msigdb, heatmaps, volcano, tfs = dirs_list

class dirs: 
    def __init__(self):
        self.batchcounts = '1. Batch Counts'
        self.deseq2input = '2. DESeq2 Input'
        self.deseq2output = '3. DESeq2 Output'
        self.sigdegs = '4. Significant DEGs'
        self.pcaplots = '5. PCA Plots'
        self.gsea = '6. GSEA Plots'
        self.msigdb = '6. msigdb gene sets'
        self.heatmaps = '7. Heatmaps'
        self.volcano = '8. Volcano Plots'
        self.tfs = '9. TF Analysis'
        
        for dir in dirs_list:
            os.makedirs(dir, exist_ok=True)

def read_sigs(updown='Up',contrast='D35'):
    if updown != None: sigs_df = pd.read_excel(f'{sigdegs}/{updown}regulated DEGs.xlsx',sheet_name=f'{contrast}_{updown.lower()}sigs')
    else: sigs_df = pd.read_excel(f'{sigdegs}/Significant DEG.xlsx',sheet_name=f'{contrast}_sigs')
    return sigs_df

def read_normcounts():
    normcounts_df = pd.read_excel(f'{deseq2output}/norm_counts.xlsx',header=0)
    normcounts_df.dropna(inplace=True)
    normcounts_df.rename(columns={'index':'Symbol'}, inplace = True)
    normcounts_df.drop(columns='Unnamed: 0',inplace=True)
    normcounts_df.set_index('Symbol', inplace = True)
    normcounts_df.drop_duplicates(keep='first', inplace=True)
    return normcounts_df

def read_metadata():
    metadata_df = pd.read_excel(f'{deseq2input}/metadata.xlsx')
    metadata_df.set_index('Sample', inplace = True) 
    return metadata_df

def read_rawcounts():
    rawcounts_df = pd.read_excel(f'{deseq2input}/raw_counts.xlsx')
    rawcounts_df.dropna(inplace=True)
    rawcounts_df.rename(columns={'Unnamed: 0':'Symbol'}, inplace = True)
    rawcounts_df.set_index('Symbol', inplace = True)
    rawcounts_df.drop_duplicates(keep='first', inplace=True)
    return rawcounts_df

def read_results(updown='Up',contrast='DCTRL0'):
    if updown != None: sigs_df = pd.read_excel(f'{deseq2output}/DESeq2 {updown}regulated Results.xlsx',sheet_name=f'{contrast}_{updown.lower()}results')
    else: 
        sigs_df = pd.read_excel(f'{deseq2output}/DESeq2 Results.xlsx',sheet_name=f'{contrast}_results')
    sigs_df.rename(columns={'Unnamed: 0':'Symbol'}, inplace = True)
    sigs_df.set_index('Symbol', inplace = True)
    sigs_df.drop_duplicates(keep='first', inplace=True)
    return sigs_df

def get_cytoKEGGlist():
    gene_set = {}
    for gene_set_collection in os.listdir(msigdb):
        with open(f'{msigdb}/{gene_set_collection}') as gmt:
            h = gmt.read()
        gene_set[h.split()[0]] = h.split()[2:]
    return gene_set['KEGG_REGULATION_OF_ACTIN_CYTOSKELETON']

def get_msigdb():
    gene_set = {}
    for gene_set_collection in os.listdir(msigdb):
        with open(f'{msigdb}/{gene_set_collection}') as gmt:
            h = gmt.read()
        gene_set[h.split()[0]] = h.split()[2:]

    return gene_set
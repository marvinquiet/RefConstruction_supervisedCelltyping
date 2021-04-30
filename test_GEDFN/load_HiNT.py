import os
import numpy as np
import pandas as pd

data_dir = "/home/wma36/gpu/data/GeneNetworks"
HiNTfile = data_dir+os.sep+"HomoSapiens_htb_hq.txt"

def generate_HiNT_adjacency():
    hint_df = pd.read_csv(HiNTfile, sep="\t")
    total_genes = set(hint_df['Gene_A'].tolist()).union(set(hint_df['Gene_B'].tolist())) # 10,898 genes
    adj_mat = np.zeros((len(total_genes), len(total_genes)))
    np.fill_diagonal(adj_mat, 1)
    adj_pd = pd.DataFrame(adj_mat, index=total_genes, columns=total_genes)

    for index, row in hint_df.iterrows():
        gene_A = row['Gene_A']
        gene_B = row['Gene_B']
        adj_pd.at[gene_A, gene_B] = 1

    adj_pd.to_csv(data_dir+os.sep+"HomoSapients_htb_hq_adj.txt", sep=" ")


def generate_pathway_adjacency(gmt_file):
    pathway_dict = {}
    gene_sets = set()
    with open(gmt_file, 'r') as f:
        for line in f:
            res = line.strip().split('\t')
            pathway_name = res[0]
            pathway_genes = res[2:]
            gene_sets = gene_sets.union(set(pathway_genes))
            pathway_dict[pathway_name] = pathway_genes
    gene_list = list(gene_sets)
    gene_list.sort()

    adj_mat = np.zeros((len(gene_list), len(pathway_dict)))
    adj_pd = pd.DataFrame(adj_mat, index=gene_list, columns=list(pathway_dict.keys()))

    for pathway, genes in pathway_dict.items():
        adj_pd.loc[genes, pathway] = 1
    adj_pd.to_csv(gmt_file.replace('.gmt', '')+'adj.txt', sep=" ")


def load_adjacency():
    adj_df = pd.read_csv(data_dir+os.sep+"HomoSapients_htb_hq_adj.txt", 
            sep=" ", index_col=0)
    return adj_df


if __name__ == "__main__":
    gmt_file = data_dir+os.sep+"c5.go.v7.2.symbols.gmt"
    generate_pathway_adjacency(gmt_file)

    gmt_file = data_dir+os.sep+"c2.cp.v7.2.symbols.gmt"
    generate_pathway_adjacency(gmt_file)

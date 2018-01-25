from __future__ import print_function

from pypuma import Puma

import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt

class AnalyzePuma(Puma):
    '''Network plot.'''
    def __init__(self, puma_data):
        '''Load variables from Puma.'''
        self.puma_results = puma_data.export_puma_results
        return None
    def top_network_plot(self, top = 100, file = 'puma_top_100.png'):
        '''Select top genes.'''
        subset_puma_results = self.puma_results.sort(['force'], ascending = [0])
        subset_puma_results = subset_puma_results[subset_puma_results.tf != subset_puma_results.gene]
        subset_puma_results = subset_puma_results[0:top]
        self.__shape_plot_network(subset_puma_results = subset_puma_results, file = file)
        return None
    def __shape_plot_network(self, subset_puma_results, file = 'puma.png'):
        '''Create plot.'''
        #reshape data for networkx
        unique_genes = list(set(list(subset_puma_results['tf'])+list(subset_puma_results['gene'])))
        unique_genes = pd.DataFrame(unique_genes)
        unique_genes.columns = ['name']
        unique_genes['index'] = unique_genes.index
        subset_puma_results = subset_puma_results.merge(unique_genes, how='inner', left_on='tf', right_on='name')
        subset_puma_results = subset_puma_results.rename(columns = {'index': 'tf_index'})
        subset_puma_results = subset_puma_results.drop(['name'], 1)
        subset_puma_results = subset_puma_results.merge(unique_genes, how='inner', left_on='gene', right_on='name')
        subset_puma_results = subset_puma_results.rename(columns = {'index': 'gene_index'})
        subset_puma_results = subset_puma_results.drop(['name'], 1)
        links = subset_puma_results[['tf_index', 'gene_index', 'force']]
        self.__create_plot(unique_genes = unique_genes, links = links, file = file)
        return None
    def __create_plot(self, unique_genes, links, file = 'puma.png'):
        '''Run plot.'''
        #plot
        g = nx.Graph()
        g.clear()
        plt.clf()
        g.add_nodes_from(unique_genes['index'])
        edges = []
        for i in range(0, len(links)):
            edges = edges + [(links.iloc[i]['tf_index'], links.iloc[i]['gene_index'], links.iloc[i]['force']/200)]
        g.add_weighted_edges_from(edges)
        labels = {}
        for i in range(0, len(unique_genes)):
            labels[unique_genes.iloc[i]['index']] = unique_genes.iloc[i]['name']
        pos = nx.spring_layout(g)
        nx.draw_networkx(g, pos, labels=labels,
                 node_size=20, font_size=3,
                 alpha=0.3, linewidth = 0.5, width =0.5)
        plt.axis('off')
        plt.savefig(file, dpi=300)
        return None

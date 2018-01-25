from __future__ import print_function

from pypuma import Puma
from pypuma import Lioness
from pypuma import AnalyzePuma

import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt

class AnalyzeLioness(Lioness):
    '''Network plot for Lioness data.'''
    def __init__(self, lioness_data):
        '''Load variables from lioness.'''
        self.export_puma_results = lioness_data.export_puma_results
        self.lioness_results = lioness_data.export_lioness_results
        return None
    def top_network_plot(self, column = 0, top = 100, file = 'lioness_top_100.png'):
        '''Select top genes.'''
        self.export_puma_results[['force']] = self.lioness_results.ix[:,column]
        plot = AnalyzePuma(self)
        plot.top_network_plot(top, file)
        return None

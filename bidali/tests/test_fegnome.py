from unittest import TestCase
from unittest.mock import Mock, MagicMock
import pandas as pd
from bidali.fegnome import enrichometer, fenrichmentscore
import matplotlib.pyplot as plt

# Test data for doctests
testRanks = pd.Series({'gene_a':1,'gene_b':2,'gene_c':3,'gene_d':4})
testGeneset = {'gene_b','gene_c'}

# Test cases
class test_fegnome(TestCase):
    def setUp(self):
        ranks = pd.Series(range(10),index=('gene{}'.format(i) for i in range(10)))
        self.ranks = MagicMock()
        self.ranks.__len__.side_effect = ranks.__len__
        self.ranks.min.side_effect = ranks.min
        self.ranks.max.side_effect = ranks.max
        self.ranks.__sub__.side_effect = ranks.__sub__
        self.ranks.index = ranks.index
        self.ranks.sort_values.return_value = self.ranks #return same mock when sorting
        
        genesUp = {'gene0','gene1','gene2','gene3','gene6','gene40'}
        self.genesUp = MagicMock()
        self.genesUp.__len__.side_effect = genesUp.__len__
        self.genesUp.__iter__.side_effect = genesUp.__iter__
        
    def tearDown(self):
        del self.ranks, self.genesUp

    def test_enrichometer(self):
        fig = enrichometer(ranks=self.ranks,genesUp=self.genesUp,genesDown=None)
        #plt.show()

    def test_fenrichmentscore(self):
        fenrichscores = fenrichmentscore(ranks=self.ranks,genesUp=self.genesUp,genesDown=None)
        print(fenrichmentscore)

#!/bin/env python3

from unittest.mock import Mock,MagicMock,patch,call
from unittest import TestCase
from bidali.visualizations import drawGeneEnvNetwork,drawCNAcircos
import matplotlib.pyplot as plt

class test_drawGeneEnvNetwork(TestCase):
    def setUp(self):
        self.gene = Mock()
        
    def tearDown(self):
        del self.gene

    def test_drawGeneEnvNetwork(self):
        with patch('lostdata.get_proteinNetworks') as nw:
            nw().__getattribute__ = nw().mockgetattribute
            nw().mockgetattribute.return_value = nw().stringnx
            fig = drawGeneEnvNetwork(self.gene)
            self.assertIsInstance(fig,plt.Figure)


class test_drawCNAcircos(TestCase):
    def setUp(self):
        self.cnaPositions = MagicMock(side_effect=[(36094885,83084062),(59577514,83084062)])
        
    def tearDown(self):
        del self.cnaPositions

    def test_drawCNAcircos(self):
        drawCNAcircos(self.cnaPositions,cnaTotal=10,chrRange=(26885980,83257441),genePositions={'BRIP1':61863521})


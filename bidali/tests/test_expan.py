from unittest import TestCase
from unittest.mock import Mock, MagicMock
import pandas as pd
from io import StringIO
from bidali.expan import Expan, GeneResult, GenesetResult

# Test cases
class test_expan(TestCase):
    def setUp(self):
        self.counts = StringIO(',s1,s2,s3,s4\ng1,1,2,5,6\ng2,4,5,0,0')
        self.metadata = StringIO(',ev1,ev2\ns1,a,a\ns2,a,b\ns3,b,b\ns4,b,a')
        self.annotations = pd.DataFrame({'Gene name':['g1gn','g2gn']}, index=['g1','g2'])

    def tearDown(self):
        del self.counts self.metadata

    def test_expanInit(self):
        expan = Expan(self.counts, self.metadata, annotations = self.annotations)
        self.assertIn('g1',expan.counts.index)
        self.assertIn('g1',expan.annotations.index)
        self.assertIn('s1',expan.counts.columns)
        self.assertIn('s1',expan.metadata.index)

from unittest import TestCase
from unittest.mock import Mock, MagicMock, patch
import pandas as pd
from io import StringIO
from bidali.expan import Expan, GeneResult, GenesetResult

# Test cases
class test_expan(TestCase):
    def setUp(self):
        self.counts = StringIO(',s1,s2,s3,s4\ng1,1,2,5,6\ng2,4,5,0,0')
        self.metadata = StringIO(',ev1,ev2\ns1,a,a\ns2,a,b\ns3,b,b\ns4,b,a')
        self.annotations = pd.DataFrame({'Gene name':['g1gn','g2gn']}, index=['g1','g2'])
        self.expan = Expan(self.counts, self.metadata, annotations = self.annotations)

    def tearDown(self):
        del self.counts, self.metadata, self.annotations

    def test_expanInit(self):
        self.assertIn('g1',self.expan.counts.index)
        self.assertIn('g1',self.expan.annotations.index)
        self.assertIn('s1',self.expan.counts.columns)
        self.assertIn('s1',self.expan.metadata.index)

    def test_designator(self):
        self.expan.designator('~ev1+ev2',reflevels={'ev1':'a','ev2':'a'})
        self.assertEqual(self.expan.design,'~ev1+ev2')
        self.assertListEqual(list(self.expan.designmatrix.columns),['(Intercept)','ev1b','ev2b'])

    def test_exdif(self):
        self.expan.designator('~ev1+ev2',reflevels={'ev1':'a','ev2':'a'})
        with patch('bidali.retro.DEA', return_value = (MagicMock(),self.expan.counts)) as dea:
            self.expan.exdif(contrasts = [2, 3], countfilter = 1)
            self.assertTrue(hasattr(self.expan,'results'))
            

#!/bin/env python3

from unittest.mock import Mock,MagicMock,patch,call
from unittest import TestCase
from bidali.LSD import storeDatasetLocally

class test_storeDatasetLocally(TestCase):
    def setUp(self):
        def simpleFunction(a,b=Mock()):
            return Mock()
        self.simpleFunction = simpleFunction
        test_storeDatasetLocally.simpleFunction = simpleFunction
        test_storeDatasetLocally.dependendObject1 = Mock()
        test_storeDatasetLocally.dependendObject2 = Mock()
        
    def tearDown(self):
        del self.simpleFunction

    def test_simpleFunction(self):
        wrappedFunction = storeDatasetLocally(self.simpleFunction)
        with patch('hashlib.md5') as mmd5, patch('inspect.getsource') as mgs, \
             patch('pickle.dump') as mpd, patch('pickle.load'), patch('builtins.open'):
            wrappedFunction(Mock())
        mgs.assert_called_once()

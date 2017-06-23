#!/bin/env python3

from unittest.mock import Mock,MagicMock,patch,call
from unittest import TestCase
from LSD import storeDatasetLocally

class test_storeDatasetLocally(TestCase):
    def setUp(self):
        def simpleFunction(a,b=Mock()):
            return Mock()
        self.simpleFunction = simpleFunction
        test_storeDatasetLocally.simpleFunction = simpleFunction
        test_storeDatasetLocally.dependendObject1 = Mock()
        test_storeDatasetLocally.dependendObject2 = Mock()
        
        def testFunction(a,b=Mock()):
            '''
            Dependency: test_storeDatasetLocally.simpleFunction
            Dependencies: test_storeDatasetLocally.dependendObject1 test_storeDatasetLocally.dependendObject2
            '''
            return Mock()

        self.testFunction = testFunction
        
    def tearDown(self):
        del self.simpleFunction

    def test_simpleFunction(self):
        wrappedFunction = storeDatasetLocally(self.simpleFunction)
        with patch('hashlib.md5') as mmd5, patch('inspect.getsource') as mgs, \
             patch('pickle.dump') as mpd:
            wrappedFunction(Mock())
        mgs.assert_called_once()
        mpd.assert_called()

    def test_functionWithDependency(self):
        wrappedFunction = storeDatasetLocally(self.testFunction)
        with patch('pickle.dump'): wrappedFunction(Mock())

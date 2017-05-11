#!/bin/env python3
"""
Module for Lab Speleman reporting
"""
from os.path import expanduser
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
import pandas as pd
from itertools import count

reportsDir = expanduser('~/Reports/')

class Report:
    """
    Contains as main attribute a list of sections.
    Defines methods of outputting the sections.
    """
    def __init__(self,title,intro='',conclusion='',outfile=None):
        self.sections = []
        self.title = title
        self.intro = intro
        self.conclusion = conclusion
        self.outfile = outfile

    def append(self, section):
        self.sections.append(section)

    def outputZip(self):
        """
        Outputs the report in a zip container.
        Figs and tabs as pngs and excells.
        """
        from zipfile import ZipFile
        import time
        if not self.outfile:
            self.outfile = reportsDir+time.strftime('%Y_%m_%d')+'.zip'
        with ZipFile(self.outfile, 'w') as zipcontainer:
            with zipcontainer.open('summary.txt',mode='w') as zipf:
                zipf.write('# {}\n\n{}\n{}'.format(
                    self.title,
                    self.intro,
                    ('\n## Conclusion\n' if self.conclusion else '')+self.conclusion
                ).encode())
            c = count(1)
            for section in self.sections:
                section.sectionOutZip(zipcontainer,'s{}/'.format(next(c)))

class Section:
    """
    Report sections.
    Defines methods dealing with the structure of
    sections, and section specific output.

    For adding subsections, a method is provided.
    """
    def __init__(self,title,text,
                 figures=[],tables=[],
                 subsections=[]):
        self.title = title
        self.p = text
        self.figs = figures
        self.tabs = tables
        self.subs = subsections
        self.checkValidity()

    def checkValidity(self):
        assert type(self.title) == str
        assert type(self.p) == str
        for f in self.figs:
            assert type(f) == Figure
        for t in self.tabs:
            assert type(t) == pd.DataFrame
        for s in self.subs:
            assert type(s) == Section
            s.checkValidity()
        
    def append(self,section):
        self.subs.append(section)

    def sectionOutZip(self,zipcontainer,zipdir=''):
        from io import StringIO
        with zipcontainer.open(zipdir+'section.txt',mode='w') as zipf:
            zipf.write('# {}\n{}'.format(self.title,self.p).encode())
        for f in self.figs:
            c = count(1)
            with zipcontainer.open(zipdir+'fig{}.png'.format(next(c)),mode='w') as zipf:
                f.savefig(zipf)
        for t in self.tabs:
            c = count(1)
            with zipcontainer.open(zipdir+'table{}.csv'.format(next(c)),mode='w') as zipf:
                b = StringIO()
                t.to_csv(b,sep=';',decimal=',')
                b.seek(0)
                zipf.write(b.read().encode())
        for s in self.subs:
            c = count(1)
            s.sectionOutZip(zipcontainer,'{}{}/'.format(zipdir,next(c)))

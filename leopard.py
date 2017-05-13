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
        self.title = title.strip()
        self.intro = intro.strip()
        self.conclusion = conclusion.strip()
        self.outfile = outfile

    def append(self,*args,toSection=None,**kwargs):
        """
        If toSection is None, section is appended to the main section list.
        Else if toSection is int or (int,int,...), it gets added to the subs (subsection)
        list of the specified section.

        *args and **kwargs are processed by Section class.
        See Section docs for supported arguments.
        """
        section = Section(*args,**kwargs)
        if toSection:
            if type(toSection) is int: toSection = (toSection,)
            self.sections[toSection[0]].append_subsection(section,toSection[1:])
        else:
            self.sections.append(section)

    def list(self):
        for i in range(len(self.sections)):
            self.sections[i].list(walkTrace=(i,))
        
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
                 figures=None,tables=None,
                 subsections=None):
        self.title = title.strip()
        self.p = text.strip()
        self.figs = figures if figures else []
        self.tabs = tables if tables else []
        self.subs = subsections if subsections else []
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
        
    def append_subsection(self,section,toSection=()):
        """
        toSection has to be a tuple specifying the section to
        which the subsection will be appended.
        """
        if toSection:
            self.subs[toSection[0]].append_subsection(section,toSection=toSection[1:])
        else:
            self.subs.append(section)

    @staticmethod
    def sectionWalker(section,callback,walkTrace,*args,**kwargs):
        """
        callback needs to be a function that handles different 
        Section elements appropriately
        walkTrace needs to be a tuple, indicate the route to the section, e.g. (1,2,0)
        """
        callback(section,walkTrace,*args,case='sectionmain',**kwargs)
        c = count(1)
        for f in section.figs:
            callback(section,walkTrace,*args,case='fig',element=f,**kwargs)
        c = count(1)
        for t in section.tabs:
            callback(section,walkTrace,*args,case='fig',element=t,**kwargs)
        c = count(1)
        for s in section.subs:
            Section.sectionWalker(s,callback,walkTrace+(next(c),),*args,**kwargs)

    def walkerWrapper(callback):
        def wrapper(*args,**kwargs):
            #args[0] => has to be the current walked section
            return Section.sectionWalker(args[0],callback,*args[1:],**kwargs)
        return wrapper

    @walkerWrapper
    def listContent(self,walkTrace,case=None):
        if case == 'sectionmain': print(walkTrace,self.title)   
        
    def list(self,walkTrace=()):
        def listContent(self,walkTrace,case=None):
            if case == 'sectionmain': print(walkTrace,self.title)
        Section.sectionWalker(self,listContent,walkTrace)

    def sectionOutZip(self,zipcontainer,zipdir=''):
        from io import StringIO
        with zipcontainer.open(zipdir+'section.txt',mode='w') as zipf:
            zipf.write('# {}\n{}'.format(self.title,self.p).encode())
        c = count(1)
        for f in self.figs:
            with zipcontainer.open(zipdir+'fig{}.png'.format(next(c)),mode='w') as zipf:
                f.savefig(zipf)
        c = count(1)
        for t in self.tabs:
            with zipcontainer.open(zipdir+'table{}.csv'.format(next(c)),mode='w') as zipf:
                b = StringIO()
                t.to_csv(b,sep=';',decimal=',')
                b.seek(0)
                zipf.write(b.read().encode())
        c = count(1)
        for s in self.subs:
            s.sectionOutZip(zipcontainer,'{}s{}/'.format(zipdir,next(c)))

# Helper functions
def makeFigFromFile(filename,*args,**kwargs):
    """
    Renders an image in a matplotlib figure, so it can be added to reports 
    args and kwargs are passed to plt.subplots
    """
    img = plt.imread(filename)
    fig,ax = plt.subplots(*args,**kwargs)
    ax.axis('off')
    ax.imshow(img)
    return fig

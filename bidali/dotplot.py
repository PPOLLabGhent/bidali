# -*- coding: utf-8 -*-
"""dotplot for biological syntony comparing

Module for comparitive genomics
"""
import numpy as np
import os, re
from scipy import sparse
from collections import OrderedDict
import bidali.seqanalysis as bsa
fasta_re = re.compile(r'^>(?P<seqid>.+)\n(?P<sequence>(?:[ACTG]+\n)+)', re.IGNORECASE|re.MULTILINE)

class DotPlot(object):
    """DotPlot class

    Currently using a non-sliding window, and trimming complete sequence to match this.
    Beginning and end of sequences can be biased.

    Args:
        fasta1 (str): Filename for fasta1
        fasta2 (str): Filename for fasta2
        window (int): Window size for comparing (not sliding/overlapping). Default 20.
        spacer (int): Space between two probing windows. Default 10000.
    """
    def __init__(self,fasta1,fasta2,window=20,spacer=10000):
        self.fasta1 = fasta1
        self.fasta2 = fasta2
        self.window = window
        self.spacer = spacer
        
        # seq1
        seq1 = open(self.fasta1).read()
        self.seq1 = ''
        self.seq1_contigs = OrderedDict()
        total_length = 0
        for contig in fasta_re.findall(seq1):
            self.seq1+=contig[1].replace('\n','')
            seqlen = len(self.seq1)
            self.seq1_contigs[contig[0]] = (total_length,seqlen)
            total_length = seqlen
        #Trim seq to match multiple of window size
        #self.seq1 = self.seq1[:len(self.seq1)-len(self.seq1)%self.window]
        self.seq1 = self.seq1.upper()
        
        # seq2
        seq2 = open(self.fasta2).read()
        self.seq2 = ''
        self.seq2_contigs = OrderedDict()
        total_length = 0
        for contig in fasta_re.findall(seq2):
            self.seq2+=contig[1].replace('\n','')
            seqlen = len(self.seq2)
            self.seq2_contigs[contig[0]] = (total_length,seqlen)
            total_length = seqlen
        #Trim seq to match multiple of window size
        #self.seq2 = self.seq2[:len(self.seq2)-len(self.seq2)%self.window]
        self.seq2 = self.seq2.upper()
        
        # Make dotmatrix
        self.dotmatrix = sparse.lil_matrix(#np.zeros(
            (int(len(self.seq1)/self.window),int(len(self.seq2)/self.window))
        )
        for s1i,kmer in enumerate(range(0,len(self.seq1),self.spacer)):
            s1i = int(kmer/self.window)
            kmer = self.seq1[kmer:kmer+self.window]
            kmer_c = bsa.recomplement(kmer)
            # kmer same strand
            startsearchpos = 0
            try:
                while startsearchpos < len(self.seq1):
                    pos = self.seq2.index(kmer,startsearchpos)
                    self.dotmatrix[s1i,int(pos/self.window)] = 1
                    startsearchpos=pos+self.window
            except ValueError:
                pass
            # kmer on opposing strand
            startsearchpos = 0
            try:
                while startsearchpos < len(self.seq1):
                    pos = self.seq2.index(kmer_c,startsearchpos)
                    self.dotmatrix[s1i,int(pos/self.window)] = -1
                    startsearchpos=pos+self.window
            except ValueError:
                pass

    def plot(self,markersize=5,colorReverseComplement='g'):
        """Plot dotmatrix
        
        Args:
            markersize (int): Marker size, default 5.
            colorReverseComplement (color or None): If None/false same color
              as reference strand. Default 'g'.
        """
        import matplotlib.pyplot as plt
        coo = self.dotmatrix.tocoo()
        fig = plt.figure(figsize=(10,10))
        ax = fig.add_subplot(111)
        if colorReverseComplement:
            selection = coo.data == 1
            ax.plot(coo.col[selection], coo.row[selection], 's', color='black', ms=markersize)
            ax.plot(
                coo.col[~selection], coo.row[~selection], 's', color=colorReverseComplement, ms=markersize
            )
        else:
            ax.plot(coo.col, coo.row, 's', color='black', ms=markersize)
        ax.set_xlim(0, coo.shape[1])
        ax.set_ylim(0, coo.shape[0])
        ax.xaxis.tick_top()
        ax.set_xlabel(os.path.basename(self.fasta2))
        ax.xaxis.set_label_position('top') 
        ax.invert_yaxis()
        ax.set_ylabel(os.path.basename(self.fasta1))
        plt.show()

    def sort_genome_according_to_reference(self,filename=None):
        """Sort fasta2 genome

        for every fasta2 sequence take median seq1 position for probes and median strand value
        to determine if complement is better fit or not
        output new fasta file
        """
        import pandas as pd
        coo = self.dotmatrix.tocoo()
        probes = pd.DataFrame({
            'seq1_pos': coo.row*self.window,
            'seq2_pos': coo.col*self.window,
            'strand':coo.data
        })
        def return_contig(pos):
            for contig in self.seq2_contigs:
                if self.seq2_contigs[contig][0] <= pos <= self.seq2_contigs[contig][1]:
                    return contig
            raise ValueError('No contig for position')
        probes['seq2_contig'] = probes.seq2_pos.apply(return_contig)
        contigs_sorted = probes.groupby('seq2_contig').median().sort_values('seq1_pos')
        if filename:
            seq = open(self.fasta2).read()
            seqs = {
                contig[0]:contig[1].strip()
                for contig in fasta_re.findall(seq)
            }
            with open(filename,'wt') as fout:
                for contigrow in contigs_sorted.iterrows():
                    complement = contigrow[1].strand == -1
                    fout.write('>{}{}\n'.format(
                        contigrow[0],
                        '' if not complement else ' [REV]'
                    ))
                    fout.write(
                        seqs[contigrow[0]] if not complement else bsa.recomplement(seqs[contigrow[0]])
                    )
                    fout.write('\n')
        else:
            return contigs_sorted

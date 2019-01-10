# Gabriel Besombes January 2019

# This file only includes the Wrapper class.

import pybedtools as pb
import pandas as pd
import pyfaidx as pf
import pysam as ps
import vcf
import os
from analysis_functions import *



class Wrapper:
    """
    Wrapper class for use in file opening in PrimerDesign.
    Wrapps different methods from file openers for different
    formats and different use cases.
    That alows for a unified syntax for different formats.
    """
        
    
    def __init__(self, file_name, use):
        """
        Create a Wrapper object for the specified file and use.
        """
        
        # This shows the supported format and opening method for each use
        self.supported = {
            "reference": {"fasta": self._fastaopen},
            "annotations": {"vcf": self._vcfopen,
                            "bed": self._bedtoolopen},
            "targets": {"bed": self._bedtoolopen,
                        "tmp": self._bedtoolopen}
        }[use]
        self.file_name = file_name
        self.use = use
        self.namesplitter()
        self.supported[self.type]()
    
    def namesplitter(self):
        """
        Split self.file_name to create self.path, self.gziped, self.name
        and self.type.
        """
        l = self.file_name.split("/")
        self.path = '/'.join(l[:-1])
        l = l[-1].split(".")
        if l[-1] == "gz":
            self.gziped = True
        else:
            self.gziped = False
        self.name = '.'.join(l[:-(1 + self.gziped)])
        self.type = l[-(1 + self.gziped)]
        self.analysed = None
        self.filtered = None
    
    def _vcfopen(self):
        """
        Create hidden readers self._reader and self._pyvcf_reader
        as well as the self.fetch() method.
        """
        self._reader = ps.VariantFile(self.file_name)
        self._pyvcf_reader = vcf.Reader(filename=self.file_name)
        self.fetch = self._reader.fetch
    
    def _bedtoolopen(self):
        """
        Create hidden reader self._reader, as well as the methods
        self.fetch(), self.analyse(), self.filter() and self.replace().
        """
        self._reader = pb.BedTool(self.file_name)
        self.fetch = self._bedfetch
        self.analyse = self._bedanalyse
        self.filter = self._bedfilter
        self.replace = self._bedreplace
        self._was_filtered = False
        
    def _fastaopen(self):
        """
        Create the hidden reader self._reader, as well as the method
        self.fecth().
        """
        self._reader = pf.Fasta(self.file_name)
        self.fetch = self._fastafetch
    
    def _bedfetch(self, chrom=None, start=None, stop=None):
        """
        Return all the intervals on chrom between start and stop.
        If chrom is left to None it will be done on all chromosomes.
        If start is left to None it will be done from the start of the chromosome.
        If stop is left to None it will be done until the end of the chromosome.
        """
        if not chrom:
            if not start and not stop:
                return(self._reader)
            elif not start:
                return(self._reader.filter(lambda b: b.stop<=stop))
            elif not stop:
                return(self._reader.filter(lambda b: b.start>=start))
            else:
                return(self._reader.filter(lambda b: b.start>=start and b.stop<=stop))
        else:
            if not start and not stop:
                return(self._reader.filter(lambda b: b.chrom==chrom))
            elif not start:
                return(self._reader.filter(lambda b: b.chrom==chrom and b.stop<=stop))
            elif not stop:
                return(self._reader.filter(lambda b: b.chrom==chrom and b.start>=start))
            else:
                return(self._reader.filter(lambda b: b.chrom==chrom and b.start>=start and b.stop<=stop))
    
    def _fastafetch(self, chrom=None, start=None, stop=None):
        """
        Return the sequence on chrom between start and stop.
        If chrom is left to None it will be done on all chromosomes.
        If start is left to None it will be done from the start of the chromosome.
        If stop is left to None it will be done until the end of the chromosome.
        """
        if not chrom:
            return([self._reader[X][start:stop].seq for X in self._reader.keys()])
        return(self._reader[chrom][start:stop].seq)
    
    def __getitem__(self, n):
        """
        Links to self._reader.__getitem__() if it exists.
        """
        return(self._reader.__getitem__(n))
    def __iter__(self):
        """
        Links to self._reader.__iter__() if it exists.
        """
        return(self._reader.__iter__())
    
    def _bedanalyse(self, pyvcf_reader, reference, SSRs=True, SSRs_variants=True,
                       PI=True, call_rate=True, write=False):
        """
        Analyse the content of each interval in self._reader.
        Store the resulting DataFrame in self.analysed
        If write="file", write the result to "file" as a csv.
        """
        columns = ["CHR", "START", "END"]
        if SSRs:
            columns.append("SSRs")
        if SSRs_variants:
            columns.append("SSRs_variants")
        if PI:
            columns.append("PI")
        if call_rate:
            columns.append("min_call_rate")
        df = pd.DataFrame(columns=columns)
        
        for targ in self._reader:
            d = {
                "CHR": targ.chrom,
                "START": targ.start, 
                "END": targ.stop
            }
            if SSRs_variants or PI or call_rate:
                annots = [X for X in pyvcf_reader.fetch(targ.chrom, targ.start, targ.stop)]
            
            if SSRs:
                d["SSRs"] = find_SSR(reference[targ.chrom][targ.start:targ.stop].seq)
            if SSRs_variants:
                d["SSRs_variants"] = SSR_count(annots)
            if PI:
                d["PI"] = pi(annots)
            if call_rate:
                d["min_call_rate"] = min_call_rate(annots)
                
            df = df.append(d, ignore_index=True)
            
        self.analysed = df
        
        if write:
            df.to_csv(write, index=False)
    
    def _bedfilter(self,  SSRs=True, SSRs_variants=True, PI=True, call_rate=True, write=False):
        """
        Filter self.analysed using the given arguments.
        self.analysed must be a pandas DataFrame, usually
        created by self.analyse(). Make sure to add a DataFrame in
        the right format or to run self.analyse().
        Store the resulting DataFrame in self.filtered.
        If write="file", write the result to "file" as a csv.
        """
        self.filtered = self.analysed.head(0)
        if SSRs == True:
            SSRs = {
                2: 10,
                3: 2,
                4: 2,
                5: 2,
                6: 2,
                7: 2,
                8: 2,
                9: 2,
                10: 2
            }
        if SSRs_variants == True:
            SSRs_variants = {
                2: 10,
                3: 2,
                4: 2,
                5: 2,
                6: 2,
                7: 2,
                8: 2,
                9: 2,
                10: 2
            }
        if PI == True:
            PI = 0.5
        if call_rate == True:
            call_rate = 0.5
        
        for i, row in self.analysed.iterrows():
            l = []
            
            if SSRs:
                l.append(False)
                for length in SSRs:
                    for pattern in row["SSRs"]:
                        if type(SSRs[length]) == tuple:
                            if (len(pattern) == length
                                and any(((SSRs[length][0] <= X <= SSRs[length][1]) for X in row["SSRs"][pattern]))):
                                l[-1] = True
                                break
                        else:
                            if (len(pattern) == length 
                                and any(((X >= SSRs[length]) for X in row["SSRs"][pattern]))):
                                l[-1] = True
                                break
                    else:
                        continue
                    break
                    
            if SSRs_variants:
                l.append(False)
                for length in SSRs_variants:
                    for reg in row["SSRs_variants"]:
                        for pattern in row["SSRs_variants"][reg]:
                            if type(SSRs_variants[length]) == tuple:
                                if (len(pattern) == length
                                    and any(((SSRs_variants[length][0] <= X <= SSRs_variants[length][1]) for X in row["SSRs_variants"][reg][pattern]))):
                                    l[-1] = True
                                    break
                            else:
                                if (len(pattern) == length
                                    and any(((X >= SSRs_variants[length]) for X in row["SSRs_variants"][reg][pattern]))):
                                    l[-1] = True
                                    break
                        else:
                            continue
                        break
                        
                    else:
                        continue
                    break
            
            if PI:
                if row["PI"] >= PI:
                    l.append(True)
                else:
                    l.append(False)
            
            if call_rate:
                if row["min_call_rate"] >= call_rate:
                    l.append(True)
                else:
                    l.append(False)
            
            if all(l):
                self.filtered = self.filtered.append(row, ignore_index=True)
             
        if write:
            self.filtered.to_csv(write, index=False)
        
        self._was_filtered = True
    
    def _bedreplace(self, name):
        """
        Write the content of self.filtered to a bed file.
        Replace self.file_name and associated attributes by
        this new file.
        """
        if self._was_filtered:
            self.filtered.to_csv(name, sep="\t", index=False, columns=["CHR", "START", "END"],
                                 header=False, doublequote=False)
            self.file_name = name
            self.namesplitter()
            self._reader = pb.BedTool(self.file_name)
        else:
            raise ValueError("self.filtered doesn't exist. You need to analyse and filter this file first.") 
          
            
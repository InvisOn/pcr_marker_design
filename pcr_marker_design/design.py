"""
Created on Sat Feb 20 10:29:48 2016

@author: cfljam

Bulk Design PCR Primers from NGS Variants
-----------------------------------------

Provided with

- a fasta reference
- primer design criteria in a dict
- annotations as bed intervals
- targets as bed intervals

return an iterable  of primer sets

"""


from pyfaidx import Fasta
from pybedtools import BedTool

class PrimerDesign:
    """A primer design object that is primed
    with genome reference and variant data
    """
    def __init__(self, reference):
        self.reference = Fasta(reference)


    def getseqslice(self,contig,start,end):
        myslice=self.reference[contig][start,end]
        return myslice.seq

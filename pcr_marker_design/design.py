"""
Created on Sat Feb 20 10:29:48 2016

@author: cfljam

Bulk Design PCR Primers from NGS Variants
-----------------------------------------

Provided with

- a fasta reference
- primer design criteria in a dict
- annotations as bed intervals
- variant annotations imported via pyvcf
- targets as bed intervals

return an iterable  of primer sets

"""


from pyfaidx import Fasta
from pybedtools import BedTool

class PrimerDesign:
    """A primer design object that is primed
    with genome reference and variant data
    """
    def __init__(self,reference,annot_file,desc):
        """
        Usage:  PrimerDesign(reference, annotation, description)
        Initialise a design object witha  reference assembly and
        annotation file(s)
        """
        self.reference = Fasta(reference)
        self.annotations=BedTool(annot_file)
        self.desc=desc
        self.genome=self.reference.filename.replace("fasta","fasta.fai")


    def getseqslicedict(self,target,max_size):
        """Pass a target to a designer and get a dictionary
        slice that we can pass to P3
        """
        target_int=target.slop(b=max_size,g=self.genome)
        offset=target_int[0].start
        sldic=dict(SEQUENCE_ID=self.desc)
        sldic['TARGET_ID']=str(target[0].chrom+ "_" + str(target_int[0].start) +"_" + str(target_int[0].end))
        sldic['SEQUENCE_TEMPLATE']=str(self.reference[target[0].chrom][target_int[0].start:target_int[0].end].seq)
        slice_annot=[(X.start -offset,X.length) for X in (self.annotations - target) if (X.chrom==target[0].chrom) & \
                     (X.start > target_int[0].start) & (X.end < target_int[0].end)]
        sldic['SEQUENCE_EXCLUDED_REGION']=slice_annot
        sldic['SEQUENCE_TARGET']= (target[0].start -offset,target[0].length)
        return sldic

class VcfPrimerDesign:
    """A primer design object that is primed
    with genome reference and vcf variant data
    """
    def __init__(self,reference,vcf_file,desc):
        """
        Usage:  PrimerDesign(reference, vcf.gz, description)
        Initialise a design object with a  reference assembly and
        variant file(s)
        """
        self.reference = Fasta(reference)
        self.vcf=vcf.Reader(filename=vcf_file)
        self.desc=desc
        self.genome=self.reference.filename.replace("fa","fa.fai")


    def getseqslicedict(self,target,max_size):
        """Pass a bed target to a designer and get a dictionary
        slice that we can pass to P3
        """
        target_int=target.slop(b=max_size,g=self.genome)
        target_chrom=target[0].chrom
        target_start=target_int[0].start
        target_end=target_int[0].end
        offset=target_int[0].start
        sldic=dict(SEQUENCE_ID=self.desc)
        sldic['TARGET_ID']=target_chrom + "_" + str(target_start) +"_" + str(target_end)
        sldic['SEQUENCE_TEMPLATE']=str(self.reference[target_chrom][target_start:target_end].seq)
        #slice_annot=[(X.start -offset,X.length) for X in (self.annotations - target) if (X.chrom==target[0].chrom) & \
         #           (X.start > target_int[0].start) & (X.end < target_int[0].end)]
        slice_vars=[target_chrom + " " + str(X.start)+ " " +str(X.end) for X in self.vcf.fetch(target_chrom,target_start,target_end)]
        slice_annot=BedTool("\n".join(slice_vars),from_string=True)
        slice_annot=slice_annot-target
        sldic['SEQUENCE_EXCLUDED_REGION']=[(X.start,X.length) for X in slice_annot]
        sldic['SEQUENCE_TARGET']= (target[0].start -offset,target[0].length)
        return sldic

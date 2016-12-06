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
from pcr_marker_design import run_p3 as P3
import vcf
import re

class GffPrimerDesign:
    """A primer design object that is primed
    with genome reference and vgff/gvf variant data
    """
    def __init__(self,reference,annot_file,desc):
        """
        Usage:  PrimerDesign(reference, gff annotation, variants, description)
        Initialise a design object with a  reference assembly ,
        gff annotation file(s)
        """
        self.reference = Fasta(reference)
        self.annotations=BedTool(annot_file)
        self.desc=desc
        self.genome=re.sub("fasta$","fasta.fai",re.sub("fa$","fa.fai",self.reference.filename))


    def getseqslicedict(self,target,max_size):
        """Pass a target to a designer and get a dictionary
        slice that we can pass to P3
        """
        target_int=target.slop(b=max_size,g=self.genome)
        offset=target_int[0].start
        sldic=dict(SEQUENCE_ID=self.desc)
        sldic['TARGET_ID']=str(target[0].chrom+ "_" + str(target_int[0].start) +"_" + str(target_int[0].end))
        sldic['SEQUENCE_TEMPLATE']=str(self.reference[target[0].chrom][target_int[0].start:target_int[0].end].seq)
<<<<<<< HEAD
        ### Do this with BedTool.allhits
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
        Usage:  PrimerDesign(reference, vcf, description)
        Initialise a design object with a  reference assembly ,
        and Tabixed vcf file(s)
        """
        self.reference = Fasta(reference)
        self.annotations=BedTool(vcf_file).tabix()
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
        ### Do this with BedTool.allhits
        slice_annot=[(X.start -offset,X.length) for X in (self.annotations - target) if (X.chrom==target[0].chrom) & \
=======
        slice_annot=[(X.start-offset,X.length) for X in (self.annotations - target) if (X.chrom==target[0].chrom) & \
>>>>>>> f35dc4fb5e24da002a1c2815e7a986162c93c102
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
        self.annot=vcf.Reader(filename=vcf_file)
        self.desc=desc
        self.genome=re.sub("fasta$","fasta.fai",re.sub("fa$","fa.fai",self.reference.filename))



    def getseqslicedict(self,target,max_size):
        """Pass a bed target to a designer and get a dictionary
        slice that we can pass to P3
        """
        target_int=target.slop(b=max_size,g=self.genome)
        target_chrom=target[0].chrom
        target_start=target_int[0].start
        target_end=target_int[0].end
        sldic=dict(SEQUENCE_ID=self.desc)
        sldic['TARGET_ID']=target_chrom + "_" + str(target_start) +"_" + str(target_end)
        sldic['SEQUENCE_TEMPLATE']=str(self.reference[target_chrom][target_start:target_end].seq)
        slice_vars=[target_chrom + " " + str(X.start)+ " " +str(X.end) for X in self.annot.fetch(target_chrom,target_start,target_end)]
        slice_annot=BedTool("\n".join(slice_vars),from_string=True)
        slice_annot=slice_annot-target
        sldic['SEQUENCE_EXCLUDED_REGION']=[(X.start - target_start,X.length) for X in slice_annot]
        sldic['SEQUENCE_TARGET']= (target[0].start - target_start,target[0].length)
        return sldic

def designfromvcf(bedtargets,VCFdesigner,max_size,min_size):
    """
    usage: bedfile of targets,designer obj, max , min
    pass targets as bedtool to a designer
    return a list of dicts
    """
    P3.p3_globals['PRIMER_PRODUCT_SIZE_RANGE']=[[min_size,max_size]]
    designdict=[VCFdesigner.getseqslicedict(BedTool([b]),max_size) for b in bedtargets]
    PCR_result=[P3.run_P3(X,P3.p3_globals) for X in designdict]
    return PCR_result

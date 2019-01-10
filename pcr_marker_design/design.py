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


from pyfaidx import Fasta, FastaVariant
from pybedtools import BedTool
from pcr_marker_design import run_p3 as P3
from pcr_marker_design import  umelt_service as um

import vcf
import re
import warnings


class PrimerDesign:
    """A primer design object that is primed
    with genome reference and variant data
    """
    def __init__(self, reference, annot_file, desc):
        """
        Usage:  PrimerDesign(reference, annotation, description)
        Initialise a design object witha  reference assembly and
        annotation file(s)
        """
        self.reference = Fasta(reference)
        self.annotations = BedTool(annot_file)
        self.desc = desc
        self.genome = re.sub("fasta$", "fasta.fai", re.sub("fa$", "fa.fai", self.reference.filename))

    def getseqslicedict(self, target, max_size):
        """Pass a bedtool target to a designer and get a dictionary
        slice that we can pass to P3
        """
        target_int = target.slop(b=max_size, g=self.genome)
        offset = target_int[0].start
        target_len = target_int[0].length
        sldic = dict(SEQUENCE_ID=self.desc)
        ### Target in  region format
        #### ie CHR:start-{start + length -1}
        sldic['TARGET_ID'] = str(target[0].chrom + ":" +
                                 str(target[0].start +1) + "-" +
                                 str(target[0].end))
        ### Pass the offset to allow correction
        sldic['REF_OFFSET']=offset
        sldic['SEQUENCE_TEMPLATE'] = str(self.reference[target[0].chrom][target_int[0].start:target_int[0].end].seq)
        slice_annot = [(X.start-offset, X.length) for
                       X in (self.annotations - target) if
                       (X.chrom == target[0].chrom) & (X.start > target_int[0].start) & (X.end < target_int[0].end)]
        sldic['SEQUENCE_EXCLUDED_REGION'] = slice_annot
        sldic['SEQUENCE_TARGET'] = (target[0].start - offset, target[0].length)
        return sldic


class VcfPrimerDesign:
    """A primer design object that is primed
    with genome reference and vcf variant data
    """
    def __init__(self, reference, vcf_file, desc):
        """
        Usage:  PrimerDesign(reference, vcf.gz, description)
        Initialise a design object with a  reference assembly and
        variant file(s)
        """
        self.reference = Fasta(reference)
        ## Following to be gagged https://docs.python.org/2/library/warnings.html#temporarily-suppressing-warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            self.alt=FastaVariant(reference,vcf_file,het=True, hom=True,sample=None, as_raw=True)
        self.annot = vcf.Reader(filename=vcf_file) ## Do we need bot of these? FastaVariant may suffice for snps
        self.desc = desc
        self.genome = re.sub("fasta$", "fasta.fai", re.sub("fa$", "fa.fai", self.reference.filename))

    def getseqslicedict(self, interval , max_size, flanking=True):
        """Pass an interval target to a designer and get a dictionary
        slice that we can pass to P3. Default is for design flanking a target.
        """
        ## Grab a slice for design
        if flanking:
            lst=[interval.chrom,interval.start,interval.stop]
            interval_bed=BedTool("\t".join([str(X) for X in lst]),from_string=True)
            target_int = interval_bed.slop(b=max_size, g=self.genome)[0]
        else:
            target_int = interval
        target_chrom = target_int.chrom
        target_start = target_int.start
        target_end = target_int.end
        offset = target_int.start  ## so we can adjust against the reference
        ### this ID is for the slice passed for design
        sldic = dict(SEQUENCE_ID=target_chrom + ":" + str(target_start) + "-" + str(target_end))
        sldic['REF_OFFSET'] = offset
        ### This id is for the original target
        sldic['TARGET_ID'] = target_chrom + ":" + str(interval.start ) + "-" + str(interval.end-1)
        ### Cut out the sequence from the index
        sldic['SEQUENCE_TEMPLATE'] = str(self.reference[target_chrom][target_start:target_end].seq)
        ### Build a list of the variant features for masking
        slice_vars = [target_chrom + " " + str(max(X.start,target_start)) + " " + str(min(X.end,target_int.end)) \
                      for X in self.annot.fetch(target_chrom, target_start - 200, target_end) \
                      if X.end > target_start]  # <----  Quick fix here. Need more work
        ### and turn these into a bedtool
        slice_annot = BedTool("\n".join(slice_vars), from_string=True)
        ### Remove any annotations overlapping with target
        slice_annot = slice_annot.subtract(interval_bed,A=True) ### Check this!!!
        sldic['SEQUENCE_EXCLUDED_REGION'] = [(X.start - offset, X.length) for X in slice_annot]
        if flanking:
            sldic['SEQUENCE_TARGET'] = (max_size, interval.length)
        return sldic

    def meltSlice(self, region):
        """Apply variants to an amplicon region and pass
        ref and alt consensus to uMelt web service, returning a tuple of (ref_Tm, alt_Tm)
        """
        target=region.split(':')
        coord=[int(X)  for X in target[1].split('-')]
        target_chrom=target[0]
        target_start=coord[0] -1
        target_end=coord[1]
        ref_seq=str(self.reference[target_chrom][target_start:target_end].seq)
        alt_seq=self.alt[target_chrom][target_start:target_end]
        ## Melt both
        umelt = um.UmeltService()
        refmelt = um.MeltSeq(ref_seq)
        altmelt=um.MeltSeq(alt_seq)
        ref_melt_Tm = umelt.get_helicity_info(umelt.get_response(refmelt)).get_melting_temp()
        alt_melt_Tm = umelt.get_helicity_info(umelt.get_response(altmelt)).get_melting_temp()
        return (ref_melt_Tm,alt_melt_Tm)



def designfromvcf(bedtargets, VCFdesigner, max_size, min_size):
    """
    usage: bedTool of targets,designer obj, max , min
    pass targets as bedtool to a designer
    return a list of dicts
    """
    P3.p3_globals['PRIMER_PRODUCT_SIZE_RANGE'] = [[min_size, max_size]]
    designdict = [VCFdesigner.getseqslicedict(BedTool([b]), max_size) for b in bedtargets]
    PCR_result = [P3.run_P3(X, P3.p3_globals) for X in designdict]
    return PCR_result

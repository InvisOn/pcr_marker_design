


import pybedtools as pb
import pandas as pd
import pyfaidx as pf
import pysam as ps
import vcf
import re
import primer3 as p3
import os
from wraper import *



class PrimerDesign:
    """
    A class that contains different usefull tools for
    designing, analysing and filtering primers.
    """
    
    def __init__(self, 
                 reference_file, 
                 annotation_file, 
                 description, 
                 targets_file=None, 
                 output_dir=None,
                 amplicon_size_range=None,
                 primer_size_range=None,
                 p3_globals=None):
        """
        reference_file : should be a fasta file.
        annotation_file : should be a bed or vcf file.
        targets_file : should be a bed file.
        output_file : this file can exist or will be created.
        NOTE : depending on the format used some features may
        not be available.
        
        amplicon_size_range : should be a tuple. If None given
        the default value used is (200, 300).
        primer_size_range : should be a tuple. If None given
        the default value used is (20, 25).
        NOTE : if a tagets_file is given it should use the same
        parameters for amplicon and primer size.
        
        p3_globals : dictionary of global parameters for primer3.
        If None given the default value used is : {
            'PRIMER_OPT_SIZE': 20,
            'PRIMER_PICK_INTERNAL_OLIGO': 1,
            'PRIMER_INTERNAL_MAX_SELF_END': 8,
            'PRIMER_MIN_SIZE': self.primer_size_range[0],
            'PRIMER_MAX_SIZE': self.primer_size_range[1],
            'PRIMER_OPT_TM': 60.0,
            'PRIMER_MIN_TM': 57.0,
            'PRIMER_MAX_TM': 63.0,
            'PRIMER_MIN_GC': 20.0,
            'PRIMER_MAX_GC': 80.0,
            'PRIMER_MAX_POLY_X': 100,
            'PRIMER_NUM_RETURN': 5,
            'PRIMER_INTERNAL_MAX_POLY_X': 100,
            'PRIMER_SALT_MONOVALENT': 50.0,
            'PRIMER_DNA_CONC': 50.0,
            'PRIMER_MAX_NS_ACCEPTED': 0,
            'PRIMER_MAX_SELF_ANY': 12,
            'PRIMER_MAX_SELF_END': 8,
            'PRIMER_PAIR_MAX_COMPL_ANY': 12,
            'PRIMER_PAIR_MAX_COMPL_END': 8,
            'PRIMER_PRODUCT_SIZE_RANGE': self.amplicon_size_range,
            'PRIMER_NUM_RETURN' : 2
        }
        
        """
        # Set the reference, annotations and description using the setter methods
        self.reference = reference_file
        self.annotations = annotation_file
        self.description = description
        
        # Set the amplicon and primer size range to the default or given values
        if not amplicon_size_range:
            self.amplicon_size_range = (200, 300)
        else:
            try:
                self.amplicon_size_range = tuple(amplicon_size_range)
                if not len(self.amplicon_size_range) == 2:
                    raise ValueError("""'amplicon_size_range' must be of length 2 but was set to {} of length {}"""\
                                    .format(amplicon_size_range, len(amplicon_size_range)))
            except Exception as e:
                raise type(e)("""The value {} isn't valid for 'amplicon_size_range'.
            Try something similar to the default value of: {}""".format(amplicon_size_range, (200, 300)))
                
        if not primer_size_range:
            self.primer_size_range = (20, 25)
        else:
            try:
                self.primer_size_range = tuple(primer_size_range)
                if not len(self.primer_size_range) == 2:
                    raise ValueError("""'primer_size_range' must be of length 2 but was set to {} of length {}"""\
                                    .format(primer_size_range, len(primer_size_range)))
            except Exception as e:
                raise type(e)("""The value {} isn't valid for 'primer_size_range'.
            Try something similar to the default value of: {}""".format(primer_size_range, (20, 25)))
                
        # Set a few default values and update the other values for the primer3 parameters
        # NEED TO CHECK IF THESE ARE DEFAULT PARAM FOR PRIMER3 <-----------------------------------
        self.p3_globals = {
            'PRIMER_OPT_SIZE': 20,
            'PRIMER_PICK_INTERNAL_OLIGO': 1,
            'PRIMER_INTERNAL_MAX_SELF_END': 8,
            'PRIMER_MIN_SIZE': self.primer_size_range[0],
            'PRIMER_MAX_SIZE': self.primer_size_range[1],
            'PRIMER_OPT_TM': 60.0,
            'PRIMER_MIN_TM': 57.0,
            'PRIMER_MAX_TM': 63.0,
            'PRIMER_MIN_GC': 20.0,
            'PRIMER_MAX_GC': 80.0,
            'PRIMER_MAX_POLY_X': 100,
            'PRIMER_NUM_RETURN': 5,
            'PRIMER_INTERNAL_MAX_POLY_X': 100,
            'PRIMER_SALT_MONOVALENT': 50.0,
            'PRIMER_DNA_CONC': 50.0,
            'PRIMER_MAX_NS_ACCEPTED': 0,
            'PRIMER_MAX_SELF_ANY': 12,
            'PRIMER_MAX_SELF_END': 8,
            'PRIMER_PAIR_MAX_COMPL_ANY': 12,
            'PRIMER_PAIR_MAX_COMPL_END': 8,
            'PRIMER_PRODUCT_SIZE_RANGE': self.amplicon_size_range,
            'PRIMER_NUM_RETURN' : 2
        }
        if p3_globals:
            self.p3_globals.update(p3_globals)
        
        # If no target file is given, create one using self.gettargets()
        if not targets_file:
            # Need to modify the region here in final code
#             self.gettargets(write="targets.bed")
#             self.targets = "targets.bed"
            self._have_targets = False
        else:
            self.targets = targets_file
        
        # Set the output directory to default or to the given value
        if not output_dir:
            self.output_dir = "./"
        else:
            self.output_dir = output_dir
    
    # Getter end setter methods
    # Use hidden attributes to hold the Wrapper onjects for each file
    @property
    def reference(self):
        return(self._ref)
    @reference.setter
    def reference(self, reference_file):
        self._ref = Wraper(reference_file, "reference")
        
    @property
    def annotations(self):
        return(self._ann)
    @annotations.setter
    def annotations(self, annotation_file):
        self._ann = Wraper(annotation_file, "annotations")
        
    @property
    def targets(self):
        return(self._tar)
    @targets.setter
    def targets(self, targets_file):
        self._tar = Wraper(targets_file, "targets")
        
    def gettargetdict(self, target):
        """
        Given a target as target=["chromosome", start, end]
        output a dictionary of parameters for primer3.
        """
        chrom, start, end = target
        
        # Check if target is small enough to design around
        if not (end-start) < self.amplicon_size_range[1] - 2*self.primer_size_range[0]:
            raise ValueError("""The target '{}' is in the wrong size range.
Make sure to use targets that were designed with the same amplicon and primer size range:
            - amplicon size range = {}
            - primer size range = {}""".format(target, self.amplicon_size_range, self.primer_size_range))
        
        # Check the annotation type
        elif self.annotations.type == "vcf":
            # Slop around the target and cut the parts that are outside of the chromosome
            # to make the working interval for the design
            interval = (max(0,
                            end+self.primer_size_range[0]-self.amplicon_size_range[1]),
                        min(len(self.reference[chrom]),
                            start-self.primer_size_range[0]+self.amplicon_size_range[1]))
            
            # Create a list of regions to exclude for the design
            l=[]
            for rec in self.annotations.fetch(chrom, interval[0], interval[1]):
                # If the start of the feature is before the start of the target
                if rec.start < start:
                    rec_start = max(interval[0], rec.start)
                    rec_end = min(start, rec.stop)
                    # Add the part of the feature that is inbetween the start of the interval and the target
                    l.append((rec_start-interval[0], rec_end-rec_start))
                # If the end of the feature is after the end of the target
                if rec.stop > end:
                    rec_start = max(end, rec.start)
                    rec_end = min(interval[1], rec.stop)
                    # Add the part of the feature that is inbetween the target and the end of the interval
                    l.append((rec_start-interval[0], rec_end-rec_start))
            
            # Return the corresponding dictionnary for primer3
            return({
                "SEQUENCE_ID": self.description,
                "TARGET_ID": "{}:{}-{}".format(chrom, start, end),
                "REF_OFFSET": interval[0],
                "SEQUENCE_TEMPLATE": str(self.reference[chrom][interval[0]:interval[1]].seq),
                "SEQUENCE_TARGET": (start - interval[0], end - start),
                "SEQUENCE_EXCLUDED_REGION": l
            })
        
    
    def run_p3(self, region=None):
        """
        Run primer3 using self.p3_globals and the result from
        self.gettargetdict() for each target in `region`.
        Save the result as a pandas dataframe containing the dictionary output from
        primer3 for each target in self.p3_out.
        If no region is specified run all the targets.
        """
        if self.output_dir != "stdout":
            self.p3_out = pd.DataFrame()
            if not region:
                l = []
                for targ in self.targets.fetch():
                    target_dict = self.gettargetdict([targ.chrom, targ.start, targ.stop])
                    l.append(p3.bindings.designPrimers(target_dict, self.p3_globals))
                    l[-1]["REF_OFFSET"] = target_dict["REF_OFFSET"]
                    l[-1]["CHROMOSOME"] = targ.chrom

                self.p3_out.append(pd.DataFrame(l), ignore_index=True)

            else:
                l = []
                for chrom in region:
                    for interval in region[chrom]:
                        start, end = interval
                        for targ in self.targets.fetch(chrom, start, end):
                            target_dict = self.gettargetdict([targ.chrom, targ.start, targ.stop])
                            l.append(p3.bindings.designPrimers(target_dict, self.p3_globals))
                            l[-1]["REF_OFFSET"] = target_dict["REF_OFFSET"]
                            l[-1]["CHROMOSOME"] = targ.chrom

                self.p3_out.append(pd.DataFrame(l), ignore_index=True)
                
        else:
            if not region:
                for targ in self.targets.fetch():
                    target_dict = self.gettargetdict([targ.chrom, targ.start, targ.stop])
                    d = {
                        "REF_OFFSET": target_dict["REF_OFFSET"],
                        "CHROMOSOME": targ.chrom
                    }
                    d.update(p3.bindings.designPrimers(target_dict, self.p3_globals))
                    print(d)

            else:
                for chrom in region:
                    for interval in region[chrom]:
                        start, end = interval
                        for targ in self.targets.fetch(chrom, start, end):
                            target_dict = self.gettargetdict([targ.chrom, targ.start, targ.stop])
                            l.append(p3.bindings.designPrimers(target_dict, self.p3_globals))
                            l[-1]["REF_OFFSET"] = target_dict["REF_OFFSET"]
                            l[-1]["CHROMOSOME"] = targ.chrom
                            print("\"{}\"".format("\",\"".join([str(X) for X in l])))

        
    def cleanoutput(self, columns=None):
        """
        Get a list of dictionaries to a list of list with
        the first one as the header.
        Intended use : clean up the output from self.run_p3()
        """
        self.amplicons_df = pd.DataFrame()
        self.primers_df = pd.DataFrame()
        
        for i,row in self.p3_out.iterrows():
            for j in range(0, row['PRIMER_PAIR_NUM_RETURNED']):
                amplicon = row.loc[
                    ((not any((str(Y) in X for Y in range(0, 10))) or str(j) in X)
                     for X in row.index)
                ]
                amplicon["START"] = row["PRIMER_LEFT_" + str(j)][0] + row["REF_OFFSET"]
                amplicon["END"] = row["PRIMER_RIGHT_" + str(j)][0] + row["REF_OFFSET"] + 1
                amplicon["PRIMER_LEFT_REGION"] = "{}:{}-{}".format(
                    row["CHROMOSOME"],
                    row["PRIMER_LEFT_" + str(j)][0] + row["REF_OFFSET"],
                    row["PRIMER_LEFT_" + str(j)][0] + row["REF_OFFSET"] + row["PRIMER_LEFT_" + str(j)][1]
                )
                amplicon["PRIMER_INTERNAL_REGION"] = "{}:{}-{}".format(
                    row["CHROMOSOME"],
                    row["PRIMER_INTERNAL_" + str(j)][0] + row["REF_OFFSET"],
                    row["PRIMER_INTERNAL_" + str(j)][0] + row["REF_OFFSET"] + row["PRIMER_INTERNAL_" + str(j)][1]
                )
                amplicon["PRIMER_RIGHT_REGION"] = "{}:{}-{}".format(
                    row["CHROMOSOME"],
                    row["PRIMER_RIGHT_" + str(j)][0] + row["REF_OFFSET"] - row["PRIMER_RIGHT_" + str(j)][1] + 1,
                    row["PRIMER_RIGHT_" + str(j)][0] + row["REF_OFFSET"] + 1
                )
                amplicon.index = [X.replace("_" + str(j), "") for X in amplicon.index]
                self.amplicons_df = self.amplicons_df.append(amplicon, ignore_index=True)
                
                primer_left = row.loc[
                    ((not any((str(Y) in X for Y in range(0, 10))) or "LEFT_{}".format(str(j)) in X) 
                     for X in row.index)
                ]
                primer_left["START"] = row["PRIMER_LEFT_" + str(j)][0] + row["REF_OFFSET"]
                primer_left["END"] = row["PRIMER_LEFT_" + str(j)][0] + row["REF_OFFSET"] + row["PRIMER_LEFT_" + str(j)][1]
                primer_right = row.loc[
                    ((not any((str(Y) in X for Y in range(0, 10))) or "RIGHT_{}".format(str(j)) in X) 
                     for X in row.index)
                ]
                primer_right["START"] = row["PRIMER_RIGHT_" + str(j)][0] + row["REF_OFFSET"] - row["PRIMER_RIGHT_" + str(j)][1] + 1
                primer_right["END"] = row["PRIMER_RIGHT_" + str(j)][0] + row["REF_OFFSET"] + 1
                primer_right.index = [X.replace("_RIGHT_" + str(j), "") for X in primer_right.index]
                primer_left.index = [X.replace("_LEFT_" + str(j), "") for X in primer_left.index]
                self.primers_df = self.primers_df.append(primer_left, ignore_index=True)
                self.primers_df = self.primers_df.append(primer_right, ignore_index=True)
        
        self.amplicons_df[["START", "END"]] = self.amplicons_df[["START", "END"]].astype(int)
        self.amplicons_df.drop_duplicates(inplace=True)
        self.amplicons_df.sort_values("REF_OFFSET", inplace=True)
        if columns:
            self.amplicons_df.to_csv("{}/amplicons_{}.csv".format(self.output_dir,
                                                                  self.targets.name),
                                     columns=columns,
                                     index=False)
        else:
            self.amplicons_df.to_csv("{}/amplicons_{}.csv".format(self.output_dir,
                                                                  self.targets.name),
                                     index=False)
        self.amplicons_df.to_csv("{}/amplicons_{}.bed".format(self.output_dir,
                                                              self.targets.name),
                                 columns=["CHROMOSOME", "START", "END"],
                                 header=False, index=False, sep="\t")
        self.amplicons = Wraper("{}/amplicons_{}.bed".format(self.output_dir,
                                                              self.targets.name),
                                "targets")
        
        self.primers_df[["START", "END"]] = self.primers_df[["START", "END"]].astype(int)
        self.primers_df.drop_duplicates(inplace=True)
        self.primers_df.sort_values(["CHROMOSOME", "START"], inplace=True)
        if columns:
            self.primers_df.to_csv("{}/primers_{}.csv".format(self.output_dir,
                                                              self.targets.name),
                                   columns=columns,
                                   index=False)
        else:
            self.primers_df.to_csv("{}/primers_{}.csv".format(self.output_dir,
                                                              self.targets.name),
                                   index=False)
        self.primers_df.to_csv("{}/primers_{}.bed".format(self.output_dir,
                                                          self.targets.name),
                               columns=["CHROMOSOME", "START", "END"],
                               header=False, index=False, sep="\t")
        self.primers = Wraper("{}/primers_{}.bed".format(self.output_dir,
                                                        self.targets.name),
                              "targets")
    
    def auto(self, region=None):
        """
        Run the whole primer design pipeline on the targets in `region`:
            -get each target's parameters for primer3 with self.gettargetdict()
            -run primer3 with self.run_p3()
            -clean up the output with self.cleanoutput()
            -output as csv to self.output_file with self.output()
        """
        if not self._have_targets:
            self.gettargets(region=region, write="targets.bed")
            self.targets = "{}/targets.bed".format(self.output_dir)
            self._have_targets = True
        self.run_p3(region=region)
        self.cleanoutput()
        
    def gettargets(self, region=None, write=False):
        """
        pass
        """
        l = [[False,0,0]]
        l2 = [['CHR','START','END']]
        target_size_range = [self.amplicon_size_range[0] - 2*self.primer_size_range[1],
                             self.amplicon_size_range[1] - 2*self.primer_size_range[0]]
        if write:
            f = open(write, "w+")
            
        if not region:
            for rec in self.annotations.fetch():
                if rec.start < l[-1][2] + self.primer_size_range[0] and rec.chrom == l[-1][0]:
                    l[-1][2] = rec.stop
                else:
                    i=1
                    while (i+1<len(l) and
                           l[-1][2] - target_size_range[1] < l[-i][1] and
                           l[-i][0] == l[-1][0]):
                        if l[-i][1] < l[-1][2] - target_size_range[0]:
                            target = [l[-i][1],l[-1][2]]
                            l2.append([l[-1][0], target[0], target[1]])
                            if write:
                                f.writelines("{}\t{}\t{}\n".format(l[-1][0], target[0], target[1]))
                        i+=1
                    if rec.chrom != l[-1][0]:
                        l = [[False,0,0]]
                    l.append([rec.chrom, rec.start, rec.stop])
        else:
            for chromosome in region:
                for reg in region[chromosome]:
                    if reg[0]:
                        l = [['CHR','START','END'],[False, 0, reg[0]]]
                    else:
                        l = [['CHR','START','END'],[False, 0, 0]]
                    for rec in self.annotations.fetch(contig=chromosome, start=reg[0], stop=reg[1]):
                        if rec.start < l[-1][2] + self.primer_size_range[0]:
                            l[-1][2] = rec.stop
                        else:
                            i=1
                            while (i+1<len(l) and
                                   l[-1][2] - target_size_range[1] < l[-i][1]):
                                if l[-i][1] < l[-1][2] - target_size_range[0]:
                                    target = [l[-i][1],l[-1][2]]
                                    l2.append([l[-1][0], target[0], target[1]])
                                    if write:
                                        f.writelines("{}\t{}\t{}\n".format(l[-1][0], target[0], target[1]))
                                i+=1
                            l.append([rec.chrom, rec.start, rec.stop])
        
        if write:
            f.close()
        else:
            return(l2)
    
    def targetsanalyse(self, SSRs=True, SSRs_variants=True, PI=True, call_rate=True, write=False):
        """
        Analyse the targets bed file according to parameters and using
        the data in the vcf self.annotations file.
        SSRs : If True, analyse the SSRs content of each target.
            Result for each target will look like : 
            {"PATERN": [int(repeats), ...], ...}
        SSRs_variants : If True, analyse the SSRs content of each feature
            in the vcf for each target. 
            Result for each target will look like :
            {"CHR:start-end": {"PATERN": [int(repeats), ...], ...}, ...}
        PI : If True, calculate the nucleotide diversity as a sum of the
            nucleotide diversity for each feature in the vcf.
            Result will be a float.
        call_rate : If True, analyse the minimum call rate in the target region
            in the vcf.
        write : If not false, use this given string for the name of the output
            csv file.
        """
        self.targets.analyse(self.annotations._pyvcf_reader, self.reference, SSRs=SSRs,
                             SSRs_variants=SSRs_variants, PI=PI, call_rate=call_rate, write=write)
    
    def targetsfilter(self,  SSRs=True, SSRs_variants=True, PI=True, call_rate=True, write=False):
        """
        pass
        """
        self.targets.filter(SSRs=SSRs, SSRs_variants=SSRs_variants,
                            PI=PI, call_rate=call_rate, write=write)
    
    def targetsreplace(self, name):
        """
        pass
        """
        self.targets.replace(name)
            
    def targetsdirectreplace(self, name, SSRs=True, SSRs_variants=True, PI=True, call_rate=True,
                             write_analysed=False, write_filtered=False):
        """
        pass
        """
        l=[True if X else False for X in [SSRs, SSRs_variants, PI, call_rate]]
        self.targets.analyse(self.annotations._pyvcf_reader, self.reference, SSRs=l[0],
                             SSRs_variants=l[1], PI=l[2], call_rate=l[3], write=write_analysed)
        self.targets.filter(SSRs=SSRs, SSRs_variants=SSRs_variants,
                            PI=PI, call_rate=call_rate, write=write_filtered)
        self.targets.replace(name)
    
    def ampliconsanalyse(self, SSRs=True, SSRs_variants=True, PI=True, call_rate=True, write=False):
        """
        pass
        """
        self.amplicons.analyse(self.annotations._pyvcf_reader, self.reference, SSRs=SSRs,
                             SSRs_variants=SSRs_variants, PI=PI, call_rate=call_rate, write=write)
    
    def ampliconsfilter(self,  SSRs=True, SSRs_variants=True, PI=True, call_rate=True, write=False):
        """
        pass
        """
        self.amplicons.filter(SSRs=SSRs, SSRs_variants=SSRs_variants,
                            PI=PI, call_rate=call_rate, write=write)
    
    def ampliconsreplace(self, name):
        """
        pass
        """
        self.amplicons.replace(name)
            
    def ampliconsdirectreplace(self, name, SSRs=True, SSRs_variants=True, PI=True, call_rate=True,
                             write_analysed=False, write_filtered=False):
        """
        pass
        """
        l=[True if X else False for X in [SSRs, SSRs_variants, PI, call_rate]]
        self.amplicons.analyse(self.annotations._pyvcf_reader, self.reference, SSRs=l[0],
                             SSRs_variants=l[1], PI=l[2], call_rate=l[3], write=write_analysed)
        self.amplicons.filter(SSRs=SSRs, SSRs_variants=SSRs_variants,
                            PI=PI, call_rate=call_rate, write=write_filtered)
        self.amplicons.replace(name)
    
            
    
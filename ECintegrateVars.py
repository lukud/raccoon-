#!/usr/bin/env python3

import sys
import os
import re
import logging
import vcf
from ECUtils import *
import argparse
import pysam
from collections import defaultdict
from collections import Counter

class SNVFilter(vcf.filters.Base):
    '''
    controlls pass filtering for snv
    filtersettings are hardcoded class attribtues..
    filters simply return true for pass and false for fail

    '''
    MIN_QUAL=30     #minimum variant qual (phred)
    MIN_QD=2        #minimum quality by depth
    def __init__(self):
        pass
    def __call__(self, Record):
        if Record.QUAL < self.MIN_QUAL:  return False
        #some records have missing QD values for obscure reasons, skip those
        try:
            if Record.INFO['QD'] < self.MIN_QD:  return False
        except:
            return False
        
        return True
    
class IndelFilter(vcf.filters.Base):
    '''
    controlls pass filtering for indels
    filtersettings are hardcoded class attribtues..
    filters simply return true for pass and false for fail

    '''
    MIN_QUAL=30     #minimum variant qual (phred)
    MIN_QD=2        #minimum quality by depth
    MAXLEN=6
    
    def __init__(self):
        pass
    
    def __call__(self, Record):
        if Record.QUAL < self.MIN_QUAL: return False
        #some records have missing QD values for obscure reasons, skip those
        try:
            if Record.INFO['QD'] < self.MIN_QD:  return False
        except:
            return False
        # add 1 to maxlen as the representaiton includes
        # 1 base upstream by vcf specification
        if len(Record.REF) > self.MAXLEN+1:    return False #should this really
        for i in Record.ALT:                                #be that absolute?
            if len(i.sequence)>self.MAXLEN+1:  return False        
        
        return True
    
    
class IntegrateTrackVariants(object):
    
    def __init__(self, reference, variantFiles, mappings, outPrefix):
        self.reference=     reference
        self.variantFiles=  variantFiles
        self.mappings=      pysam.AlignmentFile(mappings, 'rb')
        self.outPrefix=     outPrefix

        self.sequenceDict={}
        self.outTracker=set()
        self.statTracker=defaultdict(int)
        self.varTrack=open(os.path.join(outPrefix,\
                                        'varTrack.tsv'),'w')
        self.statCounter=Counter()

        #load reference
        logging.info("Loading reference")        

        for sequence in load_fasta(reference):
            self.sequenceDict[sequence.identifier]=sequence
        
        #open file to track inserted variants
        self.varTracker=open(os.path.join(outPrefix,\
                                        'varTrack.tsv'),'w')
        #open file to write modified fastas
        self.varIntegration=open(os.path.join(outPrefix, \
                                              'reference.varcall.integrated.fa'),'w')
        
        for vcfFile in self.variantFiles.split(","):
            modScaffs, stats=self.integrate_variants(vcfFile)
            self.outTracker=self.outTracker.union(modScaffs)
            self.statCounter+=stats

        for sequence in self.sequenceDict:
            if sequence not in self.outTracker:
                self.varIntegration.write(self.sequenceDict[sequence]\
                                          .get_fasta_string())
        for i in self.statCounter:
            print("{}\t{}".format(i,self.statCounter[i]), file=sys.stderr)
    
    def integrate_variants(self, vcfFile):
        '''
        main method to integrate variants into the sequence
        '''
        
        sample=         None
        identifier=     None
        modSequence=    ''
        refSliceStart=  0
        refSliceStop=   0
        hetRef=         re.compile('0[/\|]\d+')
        integratedSNV=  0
        integratedIndel=0
        nonCovered=     0
        filtered=       0
        outputScaff=    set()
        
        #get filter intstances:
        SNVFilterPass=SNVFilter()
        IndelFilterPass=IndelFilter()
        
        #open vcf iterator
        vcf_reader = vcf.Reader(open(vcfFile, 'r'))
        for Record in vcf_reader:
            #get sample name in first iteration, there's only one
            if sample is None:
                call=Record.samples
                for entry in call:
                    sample=entry.sample
                    
            #check if we're on a new scaffold
            if Record.CHROM != identifier:
                #make sure it's not the first one
                if modSequence:
                    modSequence+=self.sequenceDict[identifier].\
                                sequence[refSliceStart:]
                    Seq=DNASequence(identifier, modSequence)
                    self.varIntegration.write(Seq.get_fasta_string())
                    modSequence=''
                    refSliceStart=  0
                    refSliceStop=   0
                    outputScaff.add(identifier)
                identifier=Record.CHROM
                    
            #get the current genotype
            genotype=Record.genotype(sample)['GT']

            ###filtering
            
            #filter out heterozygous calls with reference variants
            if  re.match(hetRef, genotype):
                continue
            #FILTER EXPRESSIONS FOR SNV
            if Record.is_snp and SNVFilterPass(Record):
                    integratedSNV+=1
            #FILTER EXPRESSION FOR INDELS
            elif Record.is_indel and IndelFilterPass(Record):
                    integratedIndel+=1
            #skip variants that don't pass the filter
            else:
                filtered+=1
                continue
            
            ###integrate variant to seq
            
            # check if were dealing with a haploid sample
            # (hack-ishly) and just take the alt if so
            if len(genotype)==1:
                variant=Record.ALT[0].sequence
            #if it is homozygous alt, just integrate alt:
            elif len(Record.ALT) == 1:
                variant=Record.ALT[0].sequence
            else:
                #mayority rule pick from AD tag (allelic depth
                allelicDepths=Record.genotype(sample)['AD']
                maxDepth=allelicDepths.index(max(allelicDepths))
                maxDepth=maxDepth-1 #to index the Record ALT list
                #print(maxDepth,Record.ALT)
                variant=Record.ALT[maxDepth].sequence
                
            #get coverage for the intervall
            coverage=get_perfect_coverage(self.mappings, Record.CHROM,
                             Record.start, Record.end)
                
            if coverage==0:
                nonCovered+=1
                
            #get sequence from last point of integration up to current one
            refSliceStop=Record.start
            modSequence+=self.sequenceDict[Record.CHROM].\
                        sequence[refSliceStart:refSliceStop]
            #keep track of the new coordinates of variant insertion
            varInsStart=len(modSequence)-1 #0based index
            #get rid of splip on first base of chromosome
            if varInsStart<0:
                varInsStart=0
            varInsEnd=varInsStart+len(variant)
            #plug in variant
            modSequence+=variant
            #set new start
            refSliceStart=Record.end
            self.varTrack.write('{}\t{}\t{}\t{}\t{}\t{}\n'\
                                        .format(Record.CHROM, varInsStart,\
                                                varInsEnd, variant,\
                                                Record.REF, coverage))
            
        #get sequence from last var untill end of sequence
        #try-except to catch exception that occurs if no variants were called
        try:
            modSequence+=self.sequenceDict[Record.CHROM].\
                        sequence[refSliceStart:]
        except UnboundLocalError:
            logging.warning('Looks like no variants have been called!')

        #output last scaffold:
        if modSequence:
            Seq=DNASequence(identifier, modSequence)
            self.varIntegration.write(Seq.get_fasta_string())
            outputScaff.add(identifier)

        stats=Counter()
        
        stats["INTEGRATED SNV"]=integratedSNV
        stats["INTEGRATED INDELS"]=integratedIndel
        stats["UNCOVERED VARS"]=nonCovered
        stats["FILTERED VARS"]=filtered
        stats["ASSEMBLY LENGTH"]=len(modSequence)
        
        return outputScaff, stats

    
if __name__ == '__main__':
    #set up logger
    logFormat = "%(asctime)s [%(levelname)s] %(message)s"
    logging.basicConfig( stream=sys.stderr, format=logFormat)
    
    #parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('reference', type=str)
    parser.add_argument('variants', type=str)
    parser.add_argument('alignments', type=str)
    parser.add_argument('outPrefix', type=str)
    parser.add_argument('-v', '--verbosity', type=str, \
                        choices=['info','WARNING','ERROR', 'DEBUG'], default='info')
    args=parser.parse_args()
    logging.basicConfig( stream=sys.stderr, format=logFormat, level=args.verbosity)
    i=IntegrateTrackVariants(args.reference, args.variants, args.alignments,\
                             args.outPrefix)
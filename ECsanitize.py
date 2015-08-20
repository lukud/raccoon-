#!/usr/bin/env python3

import sys
import os
import re
import logging
import argparse
import pysam
from ECUtils import *

class SanitizeVariants(object):
    '''
    This class contains the main methods for sanitizing potential error
    corrections based on read depth counts of near perfect read matches
    '''
    
    def __init__(self, bam, reference, varTrack, outfilePrefix):
        self.outputScaffolds=set()
        
        input_length=0
        self.outfilePrefix=outfilePrefix
        #parse fasta into dict of seqrecords:
        self.sequences=dict()
        logging.info('loading in reference')
        for Sequence in load_fasta(reference):
            self.sequences[Sequence.identifier]=Sequence
            input_length+=len(Sequence.sequence)
        self.Alignments=pysam.AlignmentFile(bam, 'rb')
        self.varTrack=open(varTrack,'r')
        self.varSanitiation=open(os.path.join(outfilePrefix, 'reference.sanitizedVariants.fa'),'w')
        #run the whole thing
        self.sanitize_variants()
        self.Alignments.close()
        #output scaffolds on which no corrections have been made
        for sequence in self.sequences:
            if sequence not in self.outputScaffolds:
                self.varSanitiation.write(self.sequences[sequence]\
                                          .get_fasta_string())
        self.varSanitiation.close()
        
    def sanitize_variants(self):
        '''
        main method to sanitize variants 
        '''
        
        modSequence         =''
        currentIdentifier   =''
        refSliceStart       =0
        refSliceStop        =0
        nVars               =0
        rejectedVars        =0
        input_length=0
        
        for entry in self.varTrack:
            nVars+=1
            identifier, start, end, variant, \
            subVariant, coverage=entry.split('\t')
            start=int(start)
            end=int(end)
            coverage=int(coverage)

            #check if we're on a new scaffold
            if identifier != currentIdentifier:
                #make sure it's not the first one
                if modSequence:
                    #append remaining sequence
                    modSequence+=self.sequences[currentIdentifier].\
                                    sequence[refSliceStart:]
                    Seq=DNASequence(currentIdentifier, modSequence)                    
                    self.varSanitiation.write(Seq.get_fasta_string())
                    self.outputScaffolds.add(Seq.identifier)
                    refSliceStart=0
                    refSliceStop=0
                    modSequence=''
                currentIdentifier=identifier
            
            newCoverage=get_perfect_coverage(self.Alignments, identifier,
                start, end)
            if newCoverage < coverage or (coverage==0 and newCoverage==0):
                refSliceStop=start
                modSequence+=self.sequences[identifier].\
                        sequence[refSliceStart:refSliceStop]
                refSliceStart=end
                #reinsert previous allel
                modSequence+=subVariant
                
                rejectedVars+=1
        
        #get last entry:
        modSequence+=self.sequences[identifier].\
                        sequence[refSliceStart:]
        Seq=DNASequence(currentIdentifier, modSequence)
        self.varSanitiation.write(Seq.get_fasta_string())
        self.outputScaffolds.add(Seq.identifier)
        print('REJECTED {}'.format(rejectedVars))
        print(nVars)
        
if __name__ == '__main__':
    
    logFormat = "%(asctime)s [%(levelname)s] %(message)s"
    logging.basicConfig( stream=sys.stderr, format=logFormat)
    
    #parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('reference', type=str)
    parser.add_argument('alignments', type=str)
    parser.add_argument('vartrack', type=str)
    parser.add_argument('outfilePrefix', type=str)
    parser.add_argument('-v', '--verbosity', type=str, \
                        choices=['info','WARNING','ERROR', 'DEBUG'], default='info')
    args=parser.parse_args()
    logging.basicConfig( stream=sys.stderr, format=logFormat, level=args.verbosity)
    import time
    start_time = time.time()

    i=SanitizeVariants(args.alignments, args.reference, args.vartrack, args.outfilePrefix)
    print("--- %s seconds ---" % (time.time() - start_time))
#!/usr/bin/env python
import sys
import os
import re
import logging
#import pysam

class DNASequence(object):
    """
    An object to hold a DNA sequence with some methods to manipulate it
    """
    
    alphabet = set('ACTGNactgn')
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N',\
                  'a': 't','c' : 'g', 'g': 'c','t': 'a','n':'n'}
    
    def __init__(self, identifier, sequence):
        self.identifier = identifier
        self.sequence = sequence
        
        for letter in self.sequence:
            if letter not in self.alphabet:
                raise IncorrectSequenceLetter(letter,\
                                              self.__class__.__name__)
    def __len__(self):
        return len(self.sequence)
    
    def __cmp__(self, other):
        """
        return longer sequence of two
        """
        if len(self)==len(other):   return 0
        elif len(self)<len(other):  return -1
        elif len(self)>len(other):  return 1

    def get_fasta_string(self):
        """
        return a nicely formated fasta string
        """
        splitted_seq = []
        for x in range(0,len(self.sequence),80):
            splitted_seq.append(self.sequence[x:x+80])
        return ">%s\n%s\n" %(self.identifier,"\n".join(splitted_seq))
    
    def get_revcomplement(self):
        """
        return a DNASequence instance of the complement of the current sequence
        """
        complementSeq = "".join([ self.complement[letter] for letter in self.sequence ])
        return DNASequence(identifier = self.identifier+"_complement",\
                           sequence = complementSeq[::-1] )

class IncorrectSequenceLetter(ValueError):

        def __init__(self, letter, classname):
                self.message = "The sequence item %s is not found in the alphabet of class %s\n" %(letter, classname)

def load_fasta(fasta_filename):

        fd = open(fasta_filename,"r")
        sequence = ""
        for line in fd:
                if line[0]==">":
                        if len(sequence)>0:
                                try:
                                        yield DNASequence(identifier, sequence)
                                except IncorrectSequenceLetter as e:
                                        logging.warning(e.message)
                        identifier = line[1:].strip()
                        sequence = ""
                else:
                        sequence+=line.strip()
        fd.close()

        if len(sequence)>0:
                try:
                        yield DNASequence(identifier, sequence)
                except IncorrectSequenceLetter as e:
                        logging.warning(e.message)

def get_perfect_coverage(Alignments, identifier, start, end):
    """
    Returns number of reads with a full alignment (CIGAR: readLengthM )
    in the intervall chromsome:start-stop
    Alignments is an instance of a pysam.AlignmentFile
    """
    
    perfectCoverage=0
    for read in Alignments.fetch(reference=identifier, start=start, end=end):
        #check for reads that match over the whole length with no mismatches
        if read.cigarstring=='{}M'.format(len(read.query_sequence)) and\
        read.get_tag('NM') <2:
            perfectCoverage+=1
        
    return perfectCoverage
    
    
def partition_to_re_chunks(lenghtIndex, nChunks):
    """
    A simple greedy function to partition a fasta length index (given as a
    list of name-length tuples) into nChunks of roughly equal sums of lengths
    adapted from:
    http://stackoverflow.com/questions/6855394/splitting-list-in-chunks-of-balanced-weight
    """
    
    lenghtIndex = sorted(lenghtIndex, key=lambda x: x[1])
    partitions=[[] for i in range(nChunks)]
    sums= {i:0 for i in range(nChunks)}
    current=0
    for element in lenghtIndex:
        for i in sums:
            if current == sums[i]:
                partitions[i].append(element)
                break
        sums[i]+= element[1]
            #get smallest sum
        current =min(sums.values())

    return partitions
    

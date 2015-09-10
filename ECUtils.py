#!/usr/bin/env python
import sys
import os
import re
import logging
from itertools import islice
import pysam

#################
#####CLASSES#####
#################

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
                
class ParseFastQ(object):
    """Returns a read-by-read fastQ parser analogous to file.readline()"""
    """taken from: https://scipher.wordpress.com/2010/05/06/simple-python-fastq-parser/"""
    def __init__(self,filePath,headerSymbols=['@','+']):
        """Returns a read-by-read fastQ parser analogous to file.readline().
        Exmpl: parser.next()
        -OR-
        Its an iterator so you can do:
        for rec in parser:
            ... do something with rec ...
 
        rec is tuple: (seqHeader,seqStr,qualHeader,qualStr)
        """
        if filePath.endswith('.gz'):
            self._file = gzip.open(filePath)
        else:
            self._file = open(filePath, 'rU')
        self._currentLineNumber = 0
        self._hdSyms = headerSymbols
         
    def __iter__(self):
        yield self
     
    def next(self):
        """Reads in next element, parses, and does minimal verification.
        Returns: tuple: (seqHeader,seqStr,qualHeader,qualStr)"""
        # ++++ Get Next Four Lines ++++
        elemList = []
        for i in range(4):
            line = self._file.readline()
            self._currentLineNumber += 1 ## increment file position
            if line:
                elemList.append(line.strip('\n'))
            else: 
                elemList.append(None)
         
        # ++++ Check Lines For Expected Form ++++
        trues = [bool(x) for x in elemList].count(True)
        nones = elemList.count(None)
        # -- Check for acceptable end of file --
        if nones == 4:
            raise StopIteration
        # -- Make sure we got 4 full lines of data --
        assert trues == 4,\
               "** ERROR: It looks like I encountered a premature EOF or empty line.\n\
               Please check FastQ file near line number %s (plus or minus ~4 lines) and try again**" % (self._currentLineNumber)
        # -- Make sure we are in the correct "register" --
        assert elemList[0].startswith(self._hdSyms[0]),\
               "** ERROR: The 1st line in fastq element does not start with '%s'.\n\
               Please check FastQ file near line number %s (plus or minus ~4 lines) and try again**" % (self._hdSyms[0],self._currentLineNumber) 
        assert elemList[2].startswith(self._hdSyms[1]),\
               "** ERROR: The 3rd line in fastq element does not start with '%s'.\n\
               Please check FastQ file near line number %s (plus or minus ~4 lines) and try again**" % (self._hdSyms[1],self._currentLineNumber) 
        # -- Make sure the seq line and qual line have equal lengths --
        assert len(elemList[1]) == len(elemList[3]), "** ERROR: The length of Sequence data and Quality data of the last record aren't equal.\n\
               Please check FastQ file near line number %s (plus or minus ~4 lines) and try again**" % (self._currentLineNumber) 
         
        # ++++ Return fatsQ data as tuple ++++
        return tuple(elemList)
    

###################
#####FUNCTIONS#####
###################


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

def load_fastq(fastqpath):
    if fastqpath.endswith('gz') or fastqpath.endswith('gzip'):
        fastqfile=gzip.open(fastqpath,'rb')
    else:
        fastqfile=open(fastqpath,'r')
        
    fastqiter = filter(lambda l: l, fastqfile)  # skip blank lines
    fastqiter = (l.strip('\n') for l in fastqiter)  # strip trailing newlines
    while True:
        values = list(islice(fastqiter, 4))
        if len(values) == 4:
            header1,seq,header2,qual = values
        elif len(values) == 0:
            raise StopIteration
        else:
            raise EOFError("Failed to parse four lines from fastq file!")

        if header1.startswith('@') and header2.startswith('+'):
            yield header1, seq, header2, qual
        else:
            raise ValueError("Invalid header lines: %s and %s" % (header1, header2))
        

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

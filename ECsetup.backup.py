#!/usr/bin/env python3

import argparse
import sys
import os
import re
import logging
import gzip

class Setup(object):
    '''
    Setup class to chop input reads and create a reference dictionary
    '''
    
    def __init__(self, reference, readFile1, readFile2, \
                 basedir, nPairs, nChunks=None):
        self.reference=reference
        self.readFile1=readFile1
        self.readFile2=readFile2
        self.nPairs=nPairs
        self.basedir=basedir
        self.nChunks=nChunks
        
        if self.readFile2=="None":
            self.readFile2=None
        
    def run(self):
        self.check_fasta_header()
        if self.nChunks is None or self.nChunks==0:
            suffix=''
            if self.readFile1.endswith('gz')\
            or self.readFile1.endswith('gzip'):
                suffix='.gz'
            simLinkP1=os.path.join(self.basedir, 'reads.p1.fq'+suffix)
            simLinkP2=os.path.join(self.basedir, 'reads.p2.fq'+suffix)
            os.symlink(self.readFile1, simLinkP1)
            os.symlink(self.readFile2, simLinkP2)
        else:
            self.chunkReads()
    
    def check_fasta_header(self):
        
        p=re.compile('[^a-zA-Z0-9_.\-!?=+():#]')
        with open(self.reference) as r:
            for line in r:
                if line.startswith('>'):
                    line=line.rstrip()
                    if re.search(p, line.lstrip('>')) is not None:
                        logging.critical("The headers of the reference may "\
                                         "only contain the follwing characters:"\
                                         " a-z, A-Z, 0-9, _.!?=+()- . "\
                                         "Make sure there is no whitespace!\n"\
                                         "Conflicting header: {} ".format(line))
                        sys.exit(1)
        
        
        
    def chunkReads(self):
        """
        super crude function to chunk reads. No fastq validation whatsoever
        is done, fastq is assumed to be 4 lined
        """
        def splitFile(readsPerChunk, fileHandle, basedir, prefix, gz):
            lineCount=0
            dest=None
            for line in fileHandle:
                if lineCount % (readsPerChunk*4)==0:
                    if dest:
                        dest.close()
                    outFileName='{}.{}.fq'.format(lineCount, prefix)
                    wmode="wb" if gz else "w"
                    dest=open(os.path.join(self.basedir, outFileName), wmode)
                dest.write(line)
                lineCount+=1
        gz=False
        pairsPerChunk=(self.nPairs//self.nChunks)+1
        logging.info('splitting readlines in chunks of {} reads'.format(pairsPerChunk))
        
        if self.readFile1.endswith('gz')\
        or self.readFile1.endswith('gzip'):
            logging.info('Readfiles seem to be gziped')
            p1=gzip.open(self.readFile1, 'rb')
            if self.readFile2 is not None:
                p2=gzip.open(self.readFile2, 'rb')
            gz=True
        else:
            p1=open(self.readFile1)
            if self.readFile2 is not None:
                p2=open(self.readFile2)
        
        splitFile(pairsPerChunk, p1, self.basedir, 'p1', gz)
        if self.readFile2 is not None:
            splitFile(pairsPerChunk, p2, self.basedir, 'p2', gz)
        
        p1.close()
        if self.readFile2 is not None:
            p2.close()
        
if __name__ == '__main__':
    
    logFormat = "%(asctime)s [%(levelname)s] %(message)s"
    logging.basicConfig( stream=sys.stderr, format=logFormat)
    
    #parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('reference', type=str)
    parser.add_argument('pair1', type=str)
    parser.add_argument('pair2', type=str)
    parser.add_argument('outDir', type=str)
    parser.add_argument('nPairs', type=int)
    parser.add_argument('-nChunks', type=int, default=None)
    parser.add_argument('-v', '--verbosity', type=str, \
                        choices=['INFO','WARNING','ERROR', 'DEBUG'])
    args=parser.parse_args()
    i=Setup(args.reference, args.pair1, args.pair2, \
            args.outDir, args.nPairs, args.nChunks)
    i.run()

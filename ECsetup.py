#!/usr/bin/env python3

import argparse
import sys
import os
import re
import logging
import gzip
from ECUtils import load_fastq

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
            self.chunkReads(self.readFile1, 1)
            if self.readFile2 is not None:
                self.chunkReads(self.readFile2, 2)
    
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
        
    def chunkReads(self, readfile, pair):
        '''
        Method to divide the reads into equaly sized chunks
        for scattered mapping.
        '''
        #I should be counting readpairs somwhere here
        readsPerChunk=(self.nPairs//self.nChunks)+1
        print('{} divided by {} is {}'.format(self.nPairs,self.nChunks, readsPerChunk))
        logging.info('splitting read files into chunks of {} reads'.format(readsPerChunk))
        
        chunk=0
        nReads=0
        readTracker=0
        chunkFile=None
        prefix='p1' if pair==1 else 'p2'
        for read in load_fastq(readfile):
            if readTracker%readsPerChunk==0:
                readTracker=0
                if chunkFile:
                    chunkFile.close()
                chunkFileName='{}.{}.fq'.format(nReads, prefix)
                chunkFile=open(os.path.join(self.basedir,chunkFileName), 'w')
            chunkFile.write('\n'.join(read)+"\n")
            nReads+=1
            readTracker+=1
        chunkFile.close()
                               
                    
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

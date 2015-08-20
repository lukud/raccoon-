#!/usr/bin/env python


from __future__ import division, print_function
from string import Template
from CommandRunner import *
from ECUtils import *
import tempfile
import xml.etree.ElementTree as ET
import logging
import sys,os,re
import subprocess
import signal
import argparse
import string
import glob
import traceback


class Protocol(object):
    
    def __init__(self,protocol):
        """
        Take XML input and parse it. Set Instance atributes to what 
        you want them to be
        """
        self.protocol=protocol
        self.parseProtocol()
        
    def parseProtocol(self):
        
        try:
            p=ET.parse(self.protocol)
        except:
            logging.critical("Input protocol not formatted correctly at {}".format(self.protocol))
            sys.exit(1)
        
        #parse the input and check sanity
        Reference=p.find('reference')
        if Reference is None:
            logging.error('Reference not correctly provided in protocol at {}'.format(self.protocol))
        self.reference=Reference.text
        if not os.path.exists(self.reference):
            logging.critical('The provided reference at {} does not exist'.format(self.reference))
            sys.exit(1)
            
        InBaseDir=p.find('input')
        if InBaseDir is None:
            logging.critical('Input directory not correctly provided in protocol at {}'.format(self.protocol))
            sys.exit(1)
        if "baseDir" in InBaseDir.attrib:
            self.inBaseDir=InBaseDir.attrib["baseDir"]
        else:
            self.inBaseDir=''

        Pair1=InBaseDir.find('p1')
        Pair2=InBaseDir.find('p2')
        NPairs=InBaseDir.find('nPairs')
        
        if Pair1 is None:
            logging.critical('Read file not correctly provided in protocol at {}'.format(self.protocol))
            sys.exit(1)
        self.pair1=os.path.join(self.inBaseDir,Pair1.text)
        if Pair2 is None:
            self.pair2=None
        else:
            self.pair2=os.path.join(self.inBaseDir,Pair2.text)
        if NPairs is None:
            logging.critical('Must provide number of read pairs in protocol at {}\n'\
                         '\tEg: run cat yourReadFile.fq | grep \'^@\' | wc -l'.format(self.protocol))
            sys.exit(1)
            
        #self.pair1=os.path.join(self.inBaseDir,Pair1.text)
        self.nPairs=int(NPairs.text)
        absentPair=None
        if not os.path.exists(self.pair1):
            absentPair=self.pair1
        elif self.pair2 is not None and not os.path.exists(self.pair2):
            absentPair=self.pair2
        if absentPair:
            logging.critical('The provided readfile at {} does not exist'.format(absentPair))
            sys.exit(1)
        ploidy=p.find('ploidy')
        if ploidy is None:
            logging.warn('No ploidy was provided in the protocol. Assuming a diploid sample!')
            self.ploidy=2
        else:
            self.ploidy=ploidy.text
        Outdir=p.find('outputDir')
        if Outdir is None:
            self.outDir=os.getcwd()
        else:
            self.outDir=Outdir.text
        NThreads=p.find('threads')
        if NThreads is not None:
            self.nThreads=NThreads.text
        else:
            self.nThreads='1'
        Cluster=p.find('cluster')
        if Cluster==None:
            self.cmdTemplate=None
            self.nJobs=0
        else:
            if Cluster.find('nJobs') is not None:
                self.nJobs=int(Cluster.find('nJobs').text)
            else:
                self.nJobs=0
            if Cluster.find('template') is None:
                logging.critical('Cluster tag provided but template absent at {}'.format(self.protocol))
                sys.exit(1)
            else:
                self.cmdTemplate=Cluster.find('template').text
                
        paths=p.find('paths')
        #where is the base directory of our scripts?
        if paths.find('scripts') is not None:
            self.scriptBase=paths.find('scripts').text
        else:
            logging.critical("You need to provide the 'scripts' tag within the"\
                             " 'paths' tag in the protocol!")
        #set default paths
        self.bwa='bwa'
        self.samtools='samtools'
        self.picardtools='picard.jar'
        self.gatk='gatk'
        self.pyhton3='python3'
        self.java='java'
        #reset defaults if paths were defined
        if paths is not None:
            if paths.find('bwa') is not None:
                self.bwa=paths.find('bwa').text
            if paths.find('samtools') is not None:
                self.samtools=paths.find('samtools').text
            if paths.find('picardtools') is not None:
                self.picardtools=paths.find('picardtools').text
            if paths.find('gatk') is not None:
                self.gatk=paths.find('gatk').text
            if paths.find('java') is not None:
                self.java=paths.find('java').text
            if paths.find('python3') is not None:
                self.python3=paths.find('python3').text

class StageDriver(object):
    
    '''
    This object holds the command construction methods for all stages of the
    pipeline. These methods don't do any computation by themselves, but build
    the commands to call the scripts that do.
    '''
    
    def __init__(self, protocolFile, stage, piped):

        #how are we and our protocol called?
        self.name=sys.argv[0]
        self.protocolName=protocolFile
        #parse the protocol
        self.stage=stage
        self.piped=piped
        self.MyProtocol=Protocol(protocolFile)
        self.runCmd = CommandRunner(self.MyProtocol.cmdTemplate, \
                                    self.MyProtocol.nJobs)
        
        
        #detect in what iteration we are currently at
        iterationFolders=[i.replace('ITERATION_', '') \
                          for i in glob.glob(os.path.join(self.MyProtocol.outDir,\
                                                          'ITERATION_*'))]
        #catch exception if theres no iterationfolder yet
        try:
            self.iteration=max(sorted(map(lambda f: int(f), \
                                          map(lambda f: os.path.basename(f),\
                                              iterationFolders))))
        except ValueError:
            self.iteration=1
            
        print('Working on iteration {}'.format(self.iteration))
        self.baseDir=os.path.join(self.MyProtocol.outDir,\
                                  "ITERATION_{}".format(self.iteration))
        
        stages=['setup', 'index', 'map', 'merge', 'prepvarcall', 'varcall', \
                'varintegration', 'reindex', 'remap', 'remerge', 'correction']
        
        #run stage
        if stage=="setup":
            commands, wDir=self.setup()
        elif stage=="index":
            commands, wDir=self.index(self.iteration, reindex=False, piped=self.piped)
        elif stage=="reindex":
            commands, wDir=self.index(self.iteration, reindex=True, piped=self.piped)
        elif stage=="map":
            commands, wDir=self.mapping(self.iteration, remap=False)
        elif stage=="remap":
            commands, wDir=self.mapping(self.iteration, remap=True)
        elif stage=="merge":
            commands, wDir=self.merge_bam(self.iteration, remerge=False, piped=self.piped)
        elif stage=="remerge":
            commands, wDir=self.merge_bam(self.iteration, remerge=True, piped=self.piped)
        elif stage=="prepvarcall":
            commands, wDir=self.prep_varcall(self.iteration, piped=self.piped)
        elif stage=="varcall":
            commands, wDir=self.varcall(self.iteration)
        elif stage=="varintegration":
            commands, wDir=self.varintegration(self.iteration, piped=self.piped)
        elif stage=="correction":
            commands, wDir=self.correction(self.iteration, piped=self.piped)
        elif stage=="mapprep":
            commands, wDir=self.mapprep(self.iteration, piped=self.piped)
        else:
            logging.warning('Man..this error should really really not occur...')
            #because calls to non existing stages should be filtered by argparser
            sys.exit(1)
        
        #autocall next stage for unscattered stages
        if self.piped and wDir==None:
            try:
                nextStage=stages[stages.index(stage)+1]
            except IndexError:
                #if the stage is correction, the next one is index in next iter 
                nextStage='index'
            nextStageCmd=self.selfCall(nextStage)
            commands=self.constructCommand(stage, commands+nextStageCmd)
            wDir=self.baseDir
            
        
        logging.debug("CommandRunner Returned: " + 
            str(self.runCmd(commands, wDir, self.stage )) )
        
        logging.info("Finished %s Stage: %s" % (self.runCmd.runType, self.stage))
        
        
    ##############################
    ####STAGE DRIVER FUNCTIONS####
    ##############################
            
        
    def setup(self):
        '''
        create job to check headers of reference and
        chunk up readfiles if necesary
        '''
        readDir=os.path.join(self.MyProtocol.outDir,'READS')
        retCmd = list()
        
        setupScript=os.path.join(self.MyProtocol.scriptBase, \
                                     'ECsetup.py')
        if not os.path.exists(readDir):
            os.makedirs(readDir)
        else:
            logging.critical("Directory {} already exists. Please, don't touch "\
                             "any of the folderstructures created and refrain "\
                             "from storing any files within them!".format(readDir))
            sys.exit(1)
        cmd=string.Template('${ECSetupPath} ${reference} '\
                            '${pair1} ${pair2} ${outDir} ${nPairs}').\
        substitute({'ECSetupPath':  setupScript,\
                    'reference':    self.MyProtocol.reference,\
                    'pair1':        self.MyProtocol.pair1,\
                    'pair2':        self.MyProtocol.pair2,\
                    'outDir':       readDir,\
                    'nPairs':       self.MyProtocol.nPairs})
        if self.MyProtocol.nJobs is not None:
            cmd+=' -nChunks {}'.format(self.MyProtocol.nJobs)
        cmd+=';\n'
        #why does this not work??:
        #if piped:
        #    return cmd, None
        logging.debug("Running command: {}".format(cmd))
        retCmd.append(Command(cmd, 'setup', os.path.join(readDir,"setup.out"), \
                        os.path.join(readDir,"setup.err")))
        return retCmd, readDir
    
    def index(self, iteration, reindex=False, piped=False):
        '''
        create job to index a reference
        reference nomenclature:
        first pass:     reference.fa
        after varcall:  reference.varcall.integrated.fa
        after qc:       reference.varcall.PASS.fa
        '''

        retCmd=list()
        #!!create sequence dict and fai !!
        if not reindex:
            reference=os.path.join(self.baseDir,'reference.fa')
            if iteration==1:
                if not os.path.exists(self.baseDir):
                    os.mkdir(self.baseDir)
                os.symlink(self.MyProtocol.reference, reference)
                #build command
            elif iteration!=1:
                lastPass=os.path.join(self.MyProtocol.outDir, \
                                         'ITERATION_{}'.format(iteration-1),\
                                         'reference.sanitizedVariants.fa')
                os.symlink(lastPass, reference)
        else:
            reference=os.path.join(self.baseDir, 'reference.varcall.integrated.fa')
            
        if reference.endswith('fasta'):
            seqdict=reference.replace(".fasta", ".dict")
        elif reference.endswith('.fa'):
            seqdict=reference.replace(".fa", ".dict")
        
        if reindex:
            cmd=string.Template('${bwa} index ${reference};\n')\
            .substitute({'bwa':         self.MyProtocol.bwa,
                         'reference':   reference})
            
        else:
            cmd=string.Template('${bwa} index ${reference}\n'\
                                '${samtools} faidx ${reference};\n'\
                                '${java} -jar ${picard} CreateSequenceDictionary '\
                                'R=${reference} O=${seqdict};\n')\
            .substitute({'bwa':         self.MyProtocol.bwa,
                         'reference':   reference,
                         'samtools':    self.MyProtocol.samtools,
                         'java':        self.MyProtocol.java,
                         'picard':      self.MyProtocol.picardtools,
                         'seqdict':     seqdict})
            
        #if were piping stages, just retrun the command string    
        if piped:
            return cmd, None
        #else construct the command actual command object
        j='reindex' if reindex else 'index'
        stdout=os.path.join(self.baseDir, j+".out")
        stderr=os.path.join(self.baseDir, j+".err")
        retCmd.append(Command(cmd,j,stdout, stderr))            
        
        return retCmd, self.baseDir
            
        
    def mapping(self, iteration, remap=False):
        '''create mapping jobs for all pairs in reads folder'''
        retCmd=list()
        if remap:
            outfileSuffix='_remap.bam'
            reference=os.path.join(self.baseDir, 'reference.varcall.integrated.fa')
        else:
            outfileSuffix='_map.bam'
            reference=os.path.join(self.baseDir,'reference.fa')
            
        if not os.path.exists(self.baseDir):
            os.mkdir(self.baseDir) 
        #build commands
        cmdTemplate=string.Template("${bwa} mem -t ${threads} -R"\
                                    "'@RG\\tID:dummy\\tSM:dummy\\tLB:dummy' "\
                                    "${reference} ${pair1} ${pair2} |"\
                                    "${samtools} view -@ ${threads} -Shb - |"\
                                    "${samtools} sort -@ ${threads} -O bam "\
                                    "-T tmpBam -o ${outfile} -;\n")
        
        readDir=os.path.join(self   .MyProtocol.outDir,'READS')
        pairNames=set()
        if not os.path.exists(readDir):
            logging.critical("Read directory does not exist. Did you run the setup stage?")
            sys.exit(1)
        pairs=glob.glob(os.path.join(readDir, '*.fq'))
        if not pairs:
            logging.critical("There are no read pairs in {}".format(readDir))
            sys.exit(1)
        #get basenames of pairs and remove .p[12].fq
        pairs=set(os.path.basename(re.sub('.p[12].fq', '', pair)) for pair in pairs)
        for pairPrefix in pairs:
            #check if we're dealing with single end reads
            if self.MyProtocol.pair2==None:
                pair2=''
            else:
                pair2=os.path.join(readDir,pairPrefix+'.p2.fq')
            cmd=cmdTemplate.substitute({'bwa':      self.MyProtocol.bwa,\
                                        'samtools': self.MyProtocol.samtools,\
                                        'threads':  self.MyProtocol.nThreads,\
                                        'reference':reference,\
                                        'pair1':    os.path.join(readDir,\
                                                                 pairPrefix\
                                                                 +'.p1.fq'),
                                        'pair2':    pair2,
                                        'outfile':  os.path.join(self.baseDir,\
                                                                 pairPrefix+\
                                                                 outfileSuffix)})
            j='remap' if remap else 'map'
            jobname=j
            stdout=os.path.join(self.baseDir, j+'.out')
            stderr=os.path.join(self.baseDir, j+'.err')
            logging.debug('Executing {}'.format(cmd))
            Command(cmd, jobname, stdout, stderr)
            retCmd.append(Command(cmd, j, stdout, stderr))
        
        return retCmd, self.baseDir
    
    def merge_bam(self, iteration, remerge=False, piped=False):
        '''
        construct command to merge bamfiles
        *_map.bam -> first mapping
        *_remap.bam -> second mapping for coverage
        also create index
        '''
        retCmd=list()
        if not remerge:
            bamFilesList=glob.glob(os.path.join(self.baseDir,'*_map.bam'))
            bamout=os.path.join(self.baseDir,'merged.map.bam')
        else:
            bamFilesList=glob.glob(os.path.join(self.baseDir,'*_remap.bam'))
            bamout=os.path.join(self.baseDir,'merged.remap.bam')
            
        bamFiles=" ".join(bamFilesList)
        cmdTemplate=string.Template('${samtools} merge -@ ${threads} ${bamout} ${bamFiles};'\
                            '\n${samtools} index ${bamout};\n')
        cmd=cmdTemplate.substitute({'samtools': self.MyProtocol.samtools,\
                        'bamFiles': bamFiles,\
                        'bamout':   bamout,\
                        'threads':  self.MyProtocol.nThreads})
        
        if piped:
            return cmd, None
        
        j='remerge' if remerge else 'merge'
        stdout=os.path.join(self.baseDir, j+'.out')
        stderr=os.path.join(self.baseDir, j+'err')
        retCmd.append(Command(cmd, j, stdout, stderr))
        
        return retCmd, self.baseDir
    
    def prep_varcall(self, iteration, piped=False):
        '''
        run realignerTargetcreator and indelRealigner on merged bam file
        '''
        retCmd=list()
        reference=os.path.join(self.baseDir,'reference.fa')
        inputBam=os.path.join(self.baseDir,'merged.map.bam')
        outputRTC=os.path.join(self.baseDir,'realingerTargetCreator.intervals')
        realignedBam=os.path.join(self.baseDir, 'merged.map.indelrealigned.bam')
        cmdTemplate=string.Template('${java} -jar ${gatk} -nt ${threads}'\
                                    ' -T RealignerTargetCreator -R ${reference}'\
                                    ' -I ${inputBam} -o ${outputRTC};\n'\
                                    '${java} -jar ${gatk} -T IndelRealigner'\
                                    ' -R ${reference} -I ${inputBam}'\
                                    ' -targetIntervals ${outputRTC}'\
                                    ' -o ${realignedBam};\n')
        
        cmd=cmdTemplate.substitute({'java':         self.MyProtocol.java,\
                                    'gatk':         self.MyProtocol.gatk,\
                                    'threads':      self.MyProtocol.nThreads,\
                                    'reference':    reference,\
                                    'inputBam':     inputBam,\
                                    'outputRTC':    outputRTC,\
                                    'realignedBam': realignedBam})
        
        if piped:
            return cmd, None
        
        j='prep_varcall'
        stdout=os.path.join(self.baseDir, j+'.out')
        stderr=os.path.join(self.baseDir, j+'.err')
        retCmd.append(Command(cmd, j, stdout, stderr))
        
        return retCmd, self.baseDir
    
    def varcall(self, iteration):
        '''
        run variant calls
        '''
        retCmd=list()
        reference=os.path.join(self.baseDir,'reference.fa')
        inputBam=os.path.join(self.baseDir, 'merged.map.indelrealigned.bam')
        nChunks=self.MyProtocol.nJobs
        
        #read in reference.fa.fai to get length index
        lengthIndex=[]
        n=0 #for naming the output
        with open(reference+'.fai', 'r') as f:
            for line in f:
                pack=line.rstrip().split()
                name, length=pack[0:2]
                lengthIndex.append([name, int(length)])
        #chunks is a lol of name-length tuples, each list is roughly
        #equal in terms of sequence content
        chunks=partition_to_re_chunks(lengthIndex, nChunks)
        cmdTemplate=string.Template('${java} -jar ${gatk} -nt ${threads}'\
                                    ' -T UnifiedGenotyper -R ${reference}'\
                                    ' -I ${inputBam} -L ${intervals}'\
                                    ' -glm BOTH -ploidy ${ploidy} '\
                                    '-o ${varcalls};\n')
        
        for chunk in chunks:
            intervalFile=os.path.join(self.baseDir, 'intervals.{}.list'.format(n))
            with open(intervalFile,'w') as intervals:
                print("\n".join(x[0] for x in chunk), file=intervals)
                
            varcalls=os.path.join(self.baseDir, 'varcalls.raw.{}.vcf'.format(n))
            cmd=cmdTemplate.substitute({'java':         self.MyProtocol.java,\
                                        'gatk':         self.MyProtocol.gatk,\
                                        'threads':      self.MyProtocol.nThreads,\
                                        'reference':    reference,\
                                        'inputBam':     inputBam,\
                                        'varcalls':     varcalls,\
                                        'intervals':    intervalFile,\
                                        'ploidy':       self.MyProtocol.ploidy})
            
            j='varcall'
            jobname=j
            stdout=os.path.join(self.baseDir, j+'.out')
            stderr=os.path.join(self.baseDir, j+'.err')
            logging.debug('Executing {}'.format(cmd))
            Command(cmd, jobname, stdout, stderr)
            retCmd.append(Command(cmd, j, stdout, stderr))
            n+=1
        
        return retCmd, self.baseDir
    
    def varintegration(self, iteration, piped=False):
        """
        Construction of call to ECintegrateVars
        """
        retCmd=list()
        
        integrateScript=os.path.join(self.MyProtocol.scriptBase, \
                                     'ECintegrateVars.py')
        vcfs=",".join(glob.glob(os.path.join(self.baseDir, '*.vcf')))
        cmdTemplate=string.Template('${python3} ${integrateScript} ${reference}'\
                                    ' ${vcfs} ${alignments} ${outFolder};\n')
        reference=os.path.join(self.baseDir,'reference.fa')
        bam=os.path.join(self.baseDir,'merged.map.indelrealigned.bam')
        
        cmd=cmdTemplate.substitute({'python3':      self.MyProtocol.python3,\
                                    'integrateScript':  integrateScript,\
                                    'reference':        reference,\
                                    'vcfs':             vcfs,\
                                    'alignments':       bam,\
                                    'outFolder':        self.baseDir})
        
        if piped:
            return cmd, None
        
        j='var_integration'
        stdout=os.path.join(self.baseDir, j+'.out')
        stderr=os.path.join(self.baseDir, j+'.err')
        retCmd.append(Command(cmd, j, stdout, stderr))
        
        return retCmd, self.baseDir
    
    def correction(self, iteration, piped=False):
        """
        Construct call to correction script
        """
        #first, create the folder for the next iteration so when we autocall
        #the index stage it knows what to do..
        os.mkdir(os.path.join(self.MyProtocol.outDir,\
                                  "ITERATION_{}".format(self.iteration+1)))
        
        
        retCmd=list()
        
        sanitizeScript=os.path.join(self.MyProtocol.scriptBase,\
                                    'ECsanitize.py')
        integratedVars=os.path.join(self.baseDir, \
                                    'reference.varcall.integrated.fa')
        bam=os.path.join(self.baseDir, 'merged.remap.bam')
        varTrack=os.path.join(self.baseDir, 'varTrack.tsv')

        cmdTemplate=string.Template('${python3} ${sanitizeScript} ${integratedVars}'\
                                    ' ${alignments} ${varTrack} ${outFolder};\n')
        
        cmd=cmdTemplate.substitute({'python3':      self.MyProtocol.python3,\
                                    'sanitizeScript':   sanitizeScript,\
                                    'integratedVars':   integratedVars,\
                                    'alignments':       bam,\
                                    'varTrack':         varTrack,\
                                    'outFolder':        self.baseDir})
        if piped:
            return cmd, None
        
        j='sanitizing'
        stdout=os.path.join(self.baseDir, j+'.out')
        stderr=os.path.join(self.baseDir, j+'.err')
        retCmd.append(Command(cmd, j, stdout, stderr))
        
        return retCmd, self.baseDir

    
    ########################
    ####HELPER FUNCTIONS####
    ########################
    
    def selfCall(self, stage, piped=False):
        cmdTemplate=string.Template('${python3} ${ourName} ${stage} ${protocol};\n')
        cmd=cmdTemplate.substitute({'python3':      self.MyProtocol.python3,\
                                    'ourName':  self.name,\
                                    'stage':    stage,\
                                    'protocol': self.protocolName})
        
        return cmd
    
    def constructCommand(self, jobname, cmdString):
        '''
        litte helper function to construct instances of "Command"
        '''
        
        retCmd=[]
        
        stdout=os.path.join(self.baseDir, jobname+'.out')
        stderr=os.path.join(self.baseDir, jobname+'.err')
        retCmd.append(Command(cmdString, jobname, stdout, stderr))
        return retCmd
        

if __name__ == '__main__':
    #format logger
    logFormat = "%(asctime)s [%(levelname)s] %(message)s"
    logging.basicConfig( stream=sys.stderr, format=logFormat)

    #check python version:
    if sys.version_info[0]<3:
        logging.critical("This program is not compatible with python 2.x or lower."\
                         "Please run on python3!".format(sys.version_info))
        sys.exit(1)

    #parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('stage', type=str, choices=['setup', 'index', 'map',\
                                                    'merge',  'prepvarcall', \
                                                    'varcall', 'varintegration',\
                                                    'reindex', 'remap', 'remerge',\
                                                    'correction', 'mapprep'])
    parser.add_argument('protocol', type=str)
    parser.add_argument('-p', '--piped', action='store_false')
    parser.add_argument('-v', '--verbosity', type=str,\
                        choices=['info','warning','error','debug'], default='info')
    args=parser.parse_args()
    #by default, we're in autorunning mode
    piped=args.piped
    try:
        i=StageDriver(args.protocol, args.stage, piped)
        #sys.exit([0])
    except:
        logging.error("Oh well..something went terribly wrong along the way..\n"\
                      "\t\t\t\tDid you run the previous stage yet?"
                      "\n\n\n", exc_info=True)
        

"""
This module was taken from the PBSUITE (v 14.9.9) available at http://sourceforge.net/projects/pb-jelly/
It has been published with the follwing licencsing:

##################################################



Copyright (c) '2013 Baylor College of Medicine
               Contributors: Adam English (english@bcm.edu)
                Affiliation: Human Genome Sequencing Center
                        URL: http://www.hgsc.bcm.tmc.edu/
                             https://sourceforge.net/projects/pb-jelly/
                             http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0047768
                             http://www.biomedcentral.com/1471-2105/15/180
                   Citation: English, Adam C., Stephen Richards, Yi Han, Min Wang,
                             Vanesa Vee, Jiaxin Qu, Xiang Qin, et al. "Mind the
                             Gap: Upgrading Genomes with Pacific Biosciences RS
                             Long-Read Sequencing Technology." PLoS ONE 7, no. 11
                             (November 21, 2012): e47768.
                             doi:10.1371/journal.pone.0047768.
                   Citation: English, Adam C., William J. Salerno, Jeffery G.
                             Reid. "PBHoney: identyfying genomic variants via
                             long-read discordance and interrupted mapping."
                             BMC Bioinformatics 2014, 15:180 (June 10, 2014).
                             doi:10.1186/1471-2105-15-180
                '

All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met: 

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer. 
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution. 

        THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
        ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
        WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
        DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
        ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
        (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
        LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
        ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
        (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
        SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those
of the authors and should not be interpreted as representing official policies, 
either expressed or implied, of the FreeBSD Project.


##################################################
"""

from string import Template
import tempfile
import subprocess, signal, logging, os, stat, sys

class Alarm(Exception):
    pass

def alarm_handler(signum, frame):
    raise Alarm


def exe(cmd, timeout=-1):
    """
    Executes a command through the shell.
    timeout in minutes! so 1440 mean is 24 hours.
    -1 means never
    """
    proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, \
                            stderr=subprocess.STDOUT, close_fds=True)
    signal.signal(signal.SIGALRM, alarm_handler)
    if timeout > 0:
        signal.alarm(int(timeout*60))  
    try:
        stdoutVal, stderrVal =  proc.communicate()
        print(cmd)
        logging.debug("Executing {}".format(cmd))
        logging.debug("STDERR:\n{}".format(stderrVal))
        logging.debug("STDOUT:\n{}".format(stdoutVal))
        signal.alarm(0)  # reset the alarm
    except Alarm:
        logging.error(("Command was taking too long. "
                       "Automatic Timeout Initiated after %d minutes") \
                       % (timeout))
        proc.kill()
        return 214,None,None
    
    retCode = proc.returncode
    return retCode,stdoutVal,stderrVal

class Command():
    def __init__(self, cmd, jobname, stdout, stderr):
        self.cmd = cmd
        self.jobname = jobname
        self.stdout = stdout
        self.stderr = stderr
    
    def asDict(self):
        return {"CMD":self.cmd, "JOBNAME":self.jobname, \
                "STDOUT":self.stdout, "STDERR":self.stderr}
    
class CommandRunner():
    """
    Uses a command template to run stuff. This is helpful for cluster commands
    and chunking several commands together
    """
    def __init__(self, template=None, njobs=0):
        """
        template: a string that will become the template for submitting to your cluster:
            #you can also go ahead and specify a string.Template
            default is to not submit to your cluster
            ${CMD} > ${STDOUT} 2> ${STDERR}
        njobs: (0)
            for clumping commands together and submitting them in a script
        """
        if template is None:
            template = "${CMD} > ${STDOUT} 2> ${STDERR}"
            self.runType = "Running"
        else:
            self.runType = "Submitting"
        self.template = Template(template)
        self.njobs = njobs
    
    def __call__(self, cmds, wDir = None, id = None):
        """
        Executes Commands - can either be a list or a single Command
        wDir is the working directory where chunk scripts will be written
        if id is None a random identifier will be applied when chunking
        """
        if wDir is None:
            wDir = "./"
        
        if type(cmds) != list:
            cmd = self.buildCommand(cmds)
            return exe(cmd)
        
        if self.njobs == 0:
            outRet = []
            for c in cmds:
                outRet.append(exe(self.buildCommand(c)))
            return outRet
        
        if id is None:
            id = tempfile.mkstemp(dir=wDir)[1]
        
        outputRet =[]
        for chunk, commands in enumerate( partition(cmds, self.njobs) ):
            outScript = open(os.path.join(wDir, "%s_chunk%d.sh" % (id, chunk)),'w')
            outScript.write("#!/bin/bash\n\n")
            for c in commands:
                outScript.write(c.cmd+"\n")
            outScript.close()
            #Add executeable 
            existing_permissions = stat.S_IMODE(os.stat(outScript.name).st_mode)
            if not os.access(outScript.name, os.X_OK):
                new_permissions = existing_permissions | stat.S_IXUSR
                os.chmod(outScript.name, new_permissions)
                
            submit = Command(outScript.name, \
                            id + "_chunk%d" % chunk, \
                            os.path.join(wDir, id + ("_chunk%d.out" % chunk)), \
                            os.path.join(wDir, id + ("_chunk%d.err" % chunk)))
            cmd = self.buildCommand(submit)
            outputRet.append(exe(cmd))
            
        return outputRet
        
    def checkTemplate(self):
        """
        Checks that my template works okay
        """
        temp.update({"CMD":"test", \
                     "STDOUT":"testo", \
                     "STDERR":"teste", \
                     "JOBNAME":"testn"})
        try:
            w = self.template.substitute(temp)
        except KeyError:
            logging.error("Your submission template is invalid ")
            sys.exit(1)

    def buildCommand(self, cmdSetup):
        """
        substitutes a template with a Command
        """
        return self.template.substitute(cmdSetup.asDict())

def partition(n,m):
    """
    Helper function. splits list n into m partitions
    """
    p = map(lambda x: list(), range(m))
    p=list(p)
    index = 0
    for item in n:
        p[index].append(item)
        if index < m-1:
            index += 1
        else:
            index = 0
    return filter(lambda x: len(x)>0, p)

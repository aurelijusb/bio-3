# Bioinformatics, task 3
#
# @author: Aurelijus Banelis 
from Bio.Align.Applications import MafftCommandline
from Bio import Entrez
import logging
import itertools
import re
import textwrap
import subprocess
import os
import exceptions
import time
import glob

class HpvAnalyser:
    # Configuration
    gene = 'L1'
    harmfull = [16, 18, 31, 33, 35, 51, 52]
    safe = [6, 11, 40, 42, 43, 44, 57, 81]
    email = 'aurelijus@banelis.lt'
    cdhit = '/usr/bin/cdhit-est'
    cdhitFastaExtention = '-cdhit'
    batch_size = 20
    retMax = 700
    log = logging.getLogger("hpv")
    mafft_exe = "/usr/bin/mafft --localpair --maxiterate 1000"
    mafft_lib = "/usr/lib/mafft/lib/mafft"
    
    # Properties
    query = '';
    ids = []
    count = 0
    webenv = None
    queryKey = None
    commonFile = ''
    commonAlignedFile = ''
    
    def generateEntrezQuery(self):
        """ Saves and returns query to get information about all HPV types. """
        ending = ') AND ' + self.gene + '[Gene Name] ' + \
                 'AND 1000:8000[Sequence Length]'
        types = []
        for typeNr in itertools.chain(self.harmfull, self.safe):
            types.append('"Human papillomavirus type %d"[Organism]' % typeNr)
        self.query = '(' + ' OR '.join(types) + ending
        return self.query
    
    def retrieveIds(self, query=None):
        """ Saves and returns Ids for HPV types. """
        if (query == None):
            query = self.query
        Entrez.email = self.email
        handle = Entrez.esearch(db="nucleotide", term=query, retmax=self.retMax,
                                usehistory="y")
        record = Entrez.read(handle)
        self.count = int(record["Count"])
        self.ids = record["IdList"]
        self.webenv = record["WebEnv"]
        self.queryKey = record["QueryKey"] 
        handle.close()
        
    def retrieveData(self):
        """ Get Full sequences from database and convert them to Fasta """
        self._emptyCachedFiles()
        count = len(self.ids)
        for start in range(0, count, self.batch_size):
            end = min(count, start+self.batch_size)
            self.log.info("Downloading record %i-%i/%d" %
                         (start+1, end, self.count))
            fetch_handle = Entrez.efetch(db="nucleotide", rettype="gbwithparts",
                                         retmode="text", retstart=start,
                                         retmax=self.batch_size,
                                         webenv=self.webenv,
                                         query_key=self.queryKey)
            combinedGenBank = fetch_handle.read()
            fetch_handle.close()
            for i, genBankText in enumerate(combinedGenBank.split("\n//\n")):
                if (len(genBankText) > 100):
                    self.log.info("Converting to Fasta %i/%i" %
                                   (start+i, self.count))
                    self._saveToFasta(genBankText, self.ids[start + i])
        self._removeCdHitTempFiles()
        
    def _emptyCachedFiles(self):
        """ Creates empty Fasta files for each HPV type. """ 
        for typeNr in itertools.chain(self.harmfull, self.safe):
            open(self._fastaName(typeNr), 'w').close()
        
    def _fastaName(self, typeNr, suffix = ''):
        """ Generates fasta file name by type number """
        return 'HPV' + str(typeNr) + suffix + '.fa'

    def _saveToFasta(self, genBankText, id):
        """Converts GenBank format to Fasta. 
        
        Takes only gene part and save to file as per HPV type.
        """
       
        # Parsing GenBank to Fasta sequence
        try:
            parts = genBankText.split('/gene="' + self.gene)
            boundaries = re.search('(\d+)[\.|>]+(\d+)', parts[1])
            start = int(boundaries.group(1))
            end = int(boundaries.group(2))
            origin = parts[-1].split("\nORIGIN      ")[1]
            fullFasta = ''
            for row in origin.split('\n'):
                fullFasta += row[11:].replace(' ', '').upper()
            usefullFasta = fullFasta[start:end]
            usefullFasta = "\n".join(textwrap.wrap(usefullFasta, 70))

            # Parsing GenBank to Fasta header
            title = re.search('DEFINITION  (.+)\n', parts[0]).group(1)
            id = re.search('  GI:(.+)\n', parts[0]).group(1)
            type = re.search('Human papillomavirus type (\d+)', parts[0])
            type = int(type.group(1))
            header = '>gi|' + id + ':' + str(start) + '-' + str(end) + ' ' + \
                     title
            self.log.info('Converted: ' + header)
                     
            # Saving to file
            data = header + "\n" + usefullFasta + "\n"
            file = open(self._fastaName(type), 'a')
            file.write(data)
            file.close()
            return data
        except:
            self.log.warn("Not converted to Fasta: " + id)
            print genBankText
        
    def removeDuplicates(self):
        """ Use CD-HIT to remove duplicates """
        for typeNr in itertools.chain(self.harmfull, self.safe):
            input = self._fastaName(typeNr)
            output = self._fastaName(typeNr, self.cdhitFastaExtention)
            self.log.info("CD-HIT for " + input)
            return self._run(self.cdhit + ' -i ' + input + ' -o ' + output)
        
    def _run(self, program, arguments=''): 
        """ Run program and return output of it """
        proc = subprocess.Popen([program, arguments], stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE, shell=True)
        out, err = proc.communicate()
        proc.wait()
        return out + err
    
    def _read(self, file):
        """ Reads file and returns its content """
        if (exists(file)):
            f = open(file, 'r')
            data = f.read()
            f.close()
            return data
        else:
            return None
        
    def _removeCdHitTempFiles(self):
        """ Remove .clstr files """
        for file in glob.glob('HPV*.clstr'):
            os.remove(file)

    def saveToCommonFile(self, fileName="HPV-all.fa"):
        """ Save all filtered sequences to common file """
        output = open(fileName, 'w')
        for type in self.harmfull:
            input = open(self._fastaName(type, self.cdhitFastaExtention))
            output.write(input.read())
            input.close()
        for type in self.safe:
            input = open(self._fastaName(type, self.cdhitFastaExtention))
            output.write(input.read())
            input.close()
        output.close()
        self.commonFile = fileName

    def align(self, resultFile="HPV-aligned.fa"):
        """ Use Mafft to align all sequences saved to common file """
        os.environ['MAFFT_BINARIES'] = self.mafft_lib
        mafft_cline = MafftCommandline(self.mafft_exe, input=self.commonFile)
        stdout, stderr = mafft_cline()
        handle = open(resultFile, "w")
        handle.write(stdout)
        handle.close()
        self.commonAlignedFile = resultFile
    

# Running
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)
log = logging.getLogger("hpv")
startTime = time.time()
log.info("Starting")
hpvAnalyser = HpvAnalyser()

log.info("Generating query")
hpvAnalyser.generateEntrezQuery()

log.info("Retrieving Ids")
hpvAnalyser.retrieveIds()

log.info("Retrieving Corrfinates for gene")
hpvAnalyser.retrieveData()

hpvAnalyser._removeCdHitTempFiles()

log.info("Removing duplicates")
hpvAnalyser.removeDuplicates()

log.info("Saving all to common file")
hpvAnalyser.saveToCommonFile()

log.info("Aligning")
hpvAnalyser.align()

log.info("Finished: " + str(round(time.time() - startTime)) + " seconds") 
# Bioinformatics, task 3
#
# @author: Aurelijus Banelis 
from genericpath import exists
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
import Bio.Blast.NCBIWWW
import Bio.Blast.NCBIXML

from Bio import Entrez
import logging
import itertools
import re
import textwrap


class HpvAnalyser:
    # Configuration
    gene = 'L1'
    harmfull = [16, 18, 31, 33, 35, 51, 52]
    safe = [6, 11, 40, 42, 43, 44, 57, 81]
    email = 'aurelijus@banelis.lt'
    batch_size = 3
    
    # Properties
    query = '';
    ids = []
    webenv = None
    queryKey = None
    
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
        handle = Entrez.esearch(db="nucleotide", term=query, usehistory="y")
        record = Entrez.read(handle)
        for key, value in record.items():
            print key, value
        self.ids = record["IdList"]
        self.webenv = record["WebEnv"]
        self.queryKey = record["QueryKey"] 
        handle.close()
        
    def retreveData(self):
        """ Get Full sequences from database and convert them to Fasta """
        self._emptyCachedFiles()
        for start in range(0,len(self.ids),batch_size):
            end = min(count, start+batch_size)
            logger.info("Going to download record %i to %i" % (start+1, end))
            fetch_handle = Entrez.efetch(db="nucleotide", rettype="gbwithparts",
                                         retmode="text", retstart=start,
                                         retmax=self.batch_size,
                                         webenv=self.webenv,
                                         query_key=self.queryKey)
            combinedGenBank = fetch_handle.read()
            fetch_handle.close()
            for i, genBankText in enumerate(combinedGenBank.split("\n//\n")):
                logger.info("Converting to Fasta %i/%i" % (start+i, end))
                self._saveToFasta(genBankText)
        
    def _emptyCachedFiles(self):
        """ Creates empty Fasta files for each HPV type. """ 
        for typeNr in itertools.chain(self.harmfull, self.safe):
            open(self._fastaName(typeNr), 'w').close()
        
    def _fastaName(self, type):
        """ Generates fasta file name by type number """
        return 'HPV' + str(type) + '.fa'

    def _saveToFasta(self, genBankText):
        """Converts GenBank format to Fasta. 
        
        Takes only gene part and save to file as per HPV type.
        """
       
        # Parsing GenBank to Fasta sequence
        parts = genBankText.split('/gene="' + self.gene)
        if (len(parts) > 1):
            boundaries = re.search('(\d+)[\.|>]+(\d+)', parts[1])
            start = int(boundaries.group(1))
            end = int(boundaries.group(2))
            origin = parts[2].split("\nORIGIN      \n")[1]
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
            header = '>gi|' + id + ':' + str(start) + '-' + str(end) + ' ' + title
                     
            # Saving to file
            data = header + "\n" + usefullFasta + "\n"
            file = open(self._fastaName(type), 'a')
            file.write(data)
            file.close()
            return data


# Running
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)
logging.info("Starting")
hpvAnalyser = HpvAnalyser()

logging.info("Generating query")
hpvAnalyser.generateEntrezQuery()

logging.info("Retrieving Ids")
hpvAnalyser.retrieveIds()

logging.info("Retrieving Corrfinates for gene")
hpvAnalyser.retrieveData()
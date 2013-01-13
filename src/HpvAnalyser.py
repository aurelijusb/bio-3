# Bioinformatics, task 3.2
#
# @author: Aurelijus Banelis 
from Bio import AlignIO, Entrez
from Bio.Align.Applications import MafftCommandline
from genericpath import exists
import glob
import itertools
import logging
import os
import re
import subprocess
import textwrap
import time
import sys

class HpvAnalyser:
    # Configuration
    gene = 'L1'
    harmfull = [16, 18, 31, 33, 35, 51, 52]
    safe = [6, 11, 40, 42, 43, 44, 57, 81]
    email = 'aurelijus@banelis.lt'
    cdhit = '/usr/bin/cdhit'
    cd_hit_extension = '-cdhit'
    batch_size = 20
    max_amount_to_retrieve = 700
    log = logging.getLogger("hpv")
    mafft_exe = "/usr/bin/mafft --localpair --maxiterate 1000"
    mafft_lib = "/usr/lib/mafft/lib/mafft"
    probe_size = 30
    probe_approximation = 5
    bases = ['A', 'C', 'T', 'G']
    
    # Properties
    _query = '';
    _ids = []
    _count = 0
    _web_env = None
    _query_key = None
    _unaligned_file = ''
    _aligned_file = 'HPV-aligned.fa'
    _sequences = []
    _conservations = {                      # Data structure:
                      'all':{},             #  _var[type][x][bp] = count
                      'harmfull':{},        #  type = all | harmfull | safe 
                      'safe':{}             #  bp = A | C | T | G
                     }                      #  x - column in sequnce
    _gaps = {       
             'all':{},                      # Data structure:
             'harmfull':{},                 #  _var[type][x] = count
             'safe':{}                      #  x - column in sequnce
            }
    _categories = ['harmful', 'safe', 'gaps_harmful', 'gaps_safe']
    
    _similarity = {}                        # {_categories}[x][y] = similarity
    
    # Main routine
    def main(self):
        """ Runs all methods to gather, align, analyse and save results """
        startTime = time.time()
        self.log.info("Starting")

        if (not exists('HPV-aligned.fa')):
            self.log.info("Generating _query")
            self.generate_entrez_query()
            
            self.log.info("Retrieving Ids")
            self.retrieve_ids()
            
            self.log.info("Retrieving Corrfinates for gene")
            self.retrieve_data()
            
            self.log.info("Removing duplicates")
            self.remove_duplicates()
            
            self.log.info("Saving all to common file")
            self.save_common_unaligned()
            
            self.log.info("Aligning")
            self.align()
        else:
            self.log.info("Using cached alignment")

        self.log.info("Loading sequences for analysis")
        self.load_aligned()
        self.show_sequences()
        
        self.log.info("Calculating similarity")
        self.update_conservation()
        self.show_conservations()
        self.consensus()
        
        difference = round(time.time() - startTime, 3)            
        self.log.info("Finished: " + str(difference) + " seconds") 
    
    
    # Methods
    def generate_entrez_query(self):
        """ Saves and returns _query to get information about all HPV types.
        """
        ending = ') AND ' + self.gene + '[Gene Name] ' + \
                 'AND 1000:8000[Sequence Length]'
        types = []
        for hpv_type in itertools.chain(self.harmfull, self.safe):
            types.append('"Human papillomavirus type %d"[Organism]' % hpv_type)
        self._query = '(' + ' OR '.join(types) + ending
        return self._query
    
    def retrieve_ids(self, query=None):
        """ Saves and returns Ids for HPV types.
        """
        if (query == None):
            query = self._query
        Entrez.email = self.email
        handle = Entrez.esearch(db="nucleotide", term=query,
                                retmax=self.max_amount_to_retrieve,
                                usehistory="y")
        record = Entrez.read(handle)
        self._count = int(record["Count"])
        self._ids = record["IdList"]
        self._web_env = record["WebEnv"]
        self._query_key = record["QueryKey"] 
        handle.close()
        
    def retrieve_data(self):
        """ Get Full sequences from database and convert them to Fasta.
        """
        self._empty_cached_files()
        count = len(self._ids)
        for start in range(0, count, self.batch_size):
            end = min(count, start + self.batch_size)
            self.log.info("Downloading record %i-%i/%d" % 
                          (start + 1, end, self._count))
            fetch_handle = Entrez.efetch(db="nucleotide", rettype="gbwithparts",
                                         retmode="text", retstart=start,
                                         retmax=self.batch_size,
                                         webenv=self._web_env,
                                         query_key=self._query_key)
            combined_gen_bank = fetch_handle.read()
            fetch_handle.close()
            for i, text in enumerate(combined_gen_bank.split("\n//\n")):
                if (len(text) > 100):
                    self.log.info("Converting to Fasta %i/%i" % 
                                   (start + i, self._count))
                    self._save_to_fasta(text, self._ids[start + i])
        
    def _empty_cached_files(self):
        """ Creates empty Fasta files for each HPV type.
        """ 
        for hpv_type in itertools.chain(self.harmfull, self.safe):
            open(self._fasta_name(hpv_type), 'w').close()
        
    def _fasta_name(self, typeNr, suffix=''):
        """ Generates fasta file name by type number """
        return 'HPV' + str(typeNr) + suffix + '.fa'

    def _save_to_fasta(self, gen_bank_text, sequence_id):
        """Converts GenBank format to Fasta. 
        
        Takes only gene part and save to handler as per HPV hpv_type.
        """
       
        # Parsing GenBank to Fasta sequence
        try:
            parts = gen_bank_text.split('/gene="' + self.gene)
            boundaries = re.search('(\d+)[\.|>]+(\d+)', parts[1])
            start = int(boundaries.group(1))
            end = int(boundaries.group(2))
            origin = parts[-1].split("\nORIGIN      ")[1]
            full_fasta = ''
            for row in origin.split('\n'):
                full_fasta += row[11:].replace(' ', '').upper()
            usefull_fasta = full_fasta[start:end]
            usefull_fasta = "\n".join(textwrap.wrap(usefull_fasta, 70))

            # Parsing GenBank to Fasta header
            title = re.search('DEFINITION  (.+)\n', parts[0]).group(1)
            sequence_id = re.search('  GI:(.+)\n', parts[0]).group(1)
            hpv_type = self._get_type(parts[0])
            header = '>gi|' + sequence_id + ':' + str(start) + '-' + \
                      str(end) + ' ' + title
            self.log.info('Converted: ' + header)
                     
            # Saving to handler
            data = header + "\n" + usefull_fasta + "\n"
            handler = open(self._fasta_name(hpv_type), 'a')
            handler.write(data)
            handler.close()
            return data
        except:
            self.log.warn("Not converted to Fasta: " + sequence_id)
            self.log.info(gen_bank_text)
            
    def _get_type(self, hpv_description):
        """ Extracts type number from sequence description.
        """
        hpv_description = hpv_description.replace('-', 'type')
        regex_result = re.search('Human papillomavirus type (\d+)',
                                 hpv_description)
        try:
            result = int(regex_result.group(1))
        except:
            print hpv_description
            print regex_result
            print result
        else:
            return result
        
    def remove_duplicates(self):
        """ Use CD-HIT to remove duplicates.
        """
        for hpv_type in itertools.chain(self.harmfull, self.safe):
            input_file = self._fasta_name(hpv_type)
            output_file = self._fasta_name(hpv_type, self.cd_hit_extension)
            self.log.info("CD-HIT for " + input_file)
            return self._run(self.cdhit + ' -i ' + input_file + ' -o ' + \
                              output_file)
        self._removeCdHitTempFiles()
        
    def _run(self, program, arguments=''): 
        """ Run program and return output of it.
        """
        proc = subprocess.Popen([program, arguments], stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE, shell=True)
        out, err = proc.communicate()
        proc.wait()
        return out + err
    
    def _read(self, file_name):
        """ Reads file_name and returns its content.
        """
        if (exists(file_name)):
            f = open(file_name, 'r')
            data = f.read()
            f.close()
            return data
        else:
            return None
        
    def _removeCdHitTempFiles(self):
        """ Remove .clstr files.
        """
        for file_path in glob.glob('HPV*.clstr'):
            os.remove(file_path)

    def save_common_unaligned(self, file_name="HPV-all.fa"):
        """ Save all filtered sequences to common file.
        """
        output_file = open(file_name, 'w')
        for hpv_type in self.harmfull:
            input_file = open(self._fasta_name(hpv_type,
                                               self.cd_hit_extension))
            output_file.write(input_file.read())
            input_file.close()
        for hpv_type in self.safe:
            input_file = open(self._fasta_name(hpv_type,
                                               self.cd_hit_extension))
            output_file.write(input_file.read())
            input_file.close()
        output_file.close()
        self._unaligned_file = file_name

    def align(self, result_file="HPV-aligned.fa"):
        """ Use Mafft to align all sequences saved to common file.
        """
        os.environ['MAFFT_BINARIES'] = self.mafft_lib
        mafft_cline = MafftCommandline(self.mafft_exe,
                                       input=self._unaligned_file)
        stdout, stderr = mafft_cline()
        handle = open(result_file, "w")
        handle.write(stdout + stderr)
        handle.close()
        self._aligned_file = result_file
    
    def load_aligned(self):
        """ Loads aligned sequence to futher analysis.
        """
        self._sequences = AlignIO.read(self._aligned_file, "fasta")

    def show_sequences(self):
        """ Showing aligned sequences with their HPV type.
        """
        for record in self._sequences:
            hpv_type = self._get_type(record.description)
            print "[%02d] %s\t%s\t%s" % (hpv_type, record.seq,
                           record.id, record.description)
            
    def _OLD_update_conservation(self):
        """ Caclulates similarity and gaps over all, harmfull and safe
            sequences. Resultas are saved to _conservations and _gaps.
        """
        # Default values
        categories = ['harmfull', 'safe']
        for category in categories + ['all']:
            self._conservations[category] = {}
            self._gaps[category] = {}
            for x in range(0, len(self._sequences[0])):
                self._conservations[category][x] = {}
                for base in self.bases:
                    self._conservations[category][x][base] = 0
                self._gaps[category][x] = 0
        
        # Calculating and saving
        for x in range(0, len(self._sequences[0])):
            for y in range(0, len(self._sequences)):
                hpv_type = self._get_type(self._sequences[y].description)
                if (hpv_type in self.harmfull):
                    category = 'harmfull'
                else:
                    category = 'safe'
                base = self._sequences[y,x].upper()
                if (base == '-'):
                    self._gaps[category][x] += 1
                    self._gaps['all'][x] += 1
                else:
                    self._conservations[category][x][base] += 1
                    self._conservations['all'][x][base] += 1
                    
    def _OLD_show_conservations1(self):
        """ Showing conservations horizontally bellow printed sequences.
        """
        print "Conservations"
        categories = ['harmfull', 'safe', 'all']
        for category in categories:
            for base in self.bases:
                for part in range(0, 2):
                    sys.stdout.write('     ')
                    for x in range(0, len(self._sequences[0])):
                        number = self._conservations[category][x][base]
                        symbol = self._get_numeral(number, 1 - part)
                        sys.stdout.write(symbol)
                    sys.stdout.write("\n")
                n = len(self._sequences[0]) / (len(category) + 6)
                for x in range(0, n):
                    sys.stdout.write(category + '##' + base + '###')
                sys.stdout.write("\n")

        print "Gaps"
        for category in categories:
            for part in range(0, 2):
                sys.stdout.write('     ')
                for x in range(0, len(self._sequences[0])):
                    number = self._gaps[category][x]
                    symbol = self._get_numeral(number, 1 - part)
                    sys.stdout.write(symbol)
                sys.stdout.write("\n")
            n = len(self._sequences[0]) / (len(category) + 3)
            for x in range(0, n):
                sys.stdout.write(category + '###')
            sys.stdout.write("\n")

    def update_conservation(self):
        """ Caclulates similarity each sequence symbol with safe and harmful.
        """
        # Default vaules
        n_y = len(self._sequences)
        n_x = len(self._sequences[0])
        for category in self._categories:
            self._similarity[category] = {}
            for x in range(0, n_x):
                self._similarity[category][x] = {}
                for y in range(0, n_y):
                    self._similarity[category][x][y] = 0
        
        # Filling with calculated
        for x in range(0, n_x):
            if (x % (n_x / 10) == 0):
                percent = x / (n_x / 10) * 10
                self.log.info("Calculatioing similarity: %d%%" % percent)
            for y in range(0, n_y):
                base = self._sequences[y,x].upper()
                if (base != '-'):
                    for y2 in range(0, n_y):
                        if (y != y2):
                            symbol = self._sequences[y2,x].upper()
                            hpv_type = self._get_type(self._sequences[y]
                                                            .description)
                            gap = (symbol == '-')
                            same = (base == symbol)
                            harmful = (hpv_type in self.harmfull) 
                            if (gap and harmful):
                                category = 'gaps_harmful'
                            elif (gap and not harmful):
                                category = 'gaps_safe'
                            elif (same and harmful):
                                category = 'harmful'
                            elif (same and not harmful):
                                category = 'safe'
                            else:
                                category = None
                                
                            if (category != None):
                                self._similarity[category][x][y] += 1
        
    def show_conservations(self):
        """ Showing conservations horizontally bellow printed sequences.
        """
        n_y = len(self._sequences)
        n_x = len(self._sequences[0])
        for y in range(0, n_y):
            for category in self._categories:
                for part in range(0, 2):
                    for x in range(0, n_x):
                        if (x == 0):
                            sys.stdout.write(self._short_category(category))
                        number = self._similarity[category][x][y]
                        symbol = self._get_numeral(number, 1 - part)
                        sys.stdout.write(symbol)
                    sys.stdout.write("\n")
                
    def _short_category(self, category):
        """ 4 symbol short name for category
        """
        if ("gaps" in category):
            return 'Ga_' + category[5:6] 
        else:
            return category[0:4]
        
    def _get_numeral(self, number, part):
        """ Gets symbol from number decimal representation.
        """
        if (part == 0):
            symbol = str(number % 10)
            min = 10
        elif (part == 1):
            symbol = str(number / 10 % 10)
            min = 100
        elif (part == 2):
            symbol = str(number / 100 % 10)
            min = 1000
        elif (part == 3):
            symbol = str(number / 1000 % 10)
            min = 10000
        else:
            symbol = '?'
        if (symbol == '0' and number < min):
            symbol = ' '
        return symbol
        
#    def show_conservations_png(self, file):
#        import pylab
#        categories = ['harmfull', 'safe', 'all']
#        for category in categories:
#            for base in self.bases:
#                numbers = []
#                for x in range(0, len(self._sequences[0])):
#                    number = self._conservations[category][x][base]
#                    numbers.append(number)
#                pylab.plot(numbers)
#        pylab.savefig(file)
        
    def consensus(self):
        from Bio.Align import AlignInfo
        summary_align = AlignInfo.SummaryInfo(self._sequences)
        consensus = summary_align.dumb_consensus(ambiguous = '-')
        print type(consensus)
        print '   ', consensus
        

# Running from command line
if __name__ == '__main__':
    logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)
    hpv_analyser = HpvAnalyser()
    hpv_analyser.main()

#FIXME: implement via features:
#from Bio import SeqIO
#Entrez.email = "aurelijus@banelis.lt"     # Always tell NCBI who you are
#handle = Entrez.efetch(db="nucleotide", id="186972394", rettype="gb", retmode="text")
#record = SeqIO.read(handle, "gb")
#print "%s, length %i, with %i features" % (record.name, len(record), len(record.features))
#for feature in record.features:
#    print feature   
#Bio.SeqFeature.SeqFeature

#import gtk
#import math
#
#class PyApp(gtk.Window):
#
#    def __init__(self):
#        super(PyApp, self).__init__()
#
#        self.set_title("Simple drawing")
#        self.resize(230, 150)
#        self.set_position(gtk.WIN_POS_CENTER)
#
#        self.connect("destroy", gtk.main_quit)
#
#        darea = gtk.DrawingArea()
#        darea.connect("expose-event", self.expose)
#        self.add(darea)
#
#        self.show_all()
#    
#    def expose(self, widget, event):
#        cr = widget.window.cairo_create()
##        cr.set_source_rgb(0.7, 0.2, 0.0)
##        for i in range(0, 200, 20):
##            cr.rectangle(i, 0, 19, 19)
##        cr.fill() 
#        cr.set_source_rgb(0.0, 0.0, 0.0)
#        cr.set_font_size(1.2)
#        x_bearing, y_bearing, width, height = cr.text_extents("a")[:4]
#        cr.move_to(0.5 - width / 2 - x_bearing, 0.5 - height / 2 - y_bearing)
#        cr.show_text("a")
#        
##        cr = widget.window.cairo_create()
##
##        cr.set_line_width(9)
##        cr.set_source_rgb(0.7, 0.2, 0.0)
##                
##        w = self.allocation.width
##        h = self.allocation.height
##
##        cr.translate(w/2, h/2)
##        cr.arc(0, 0, 50, 0, 2*math.pi)
##        cr.stroke_preserve()
##        
##        cr.set_source_rgb(0.3, 0.4, 0.6)
##        cr.fill()
#    
#
#PyApp()
#gtk.main()
    
    
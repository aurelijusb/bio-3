# Bioinformatics, task 3.2
#
# @author: Aurelijus Banelis 
from Bio import AlignIO, Entrez, SeqIO
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
import pickle

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
    _categories = ['harmful', 'safe', 'gaps_harmful', 'gaps_safe']
    
    _similarity = {}                        # {_categories}[x][y] = similarity
    _counts = {'harmful':0, 'safe':0}
    _candidates = {}                        # [] = [y, x1, x2]
    
    # Main routine
    def main(self):
        """ Runs all methods to gather, align, analyse and save results """
        startTime = time.time()
        self.log.info("Starting")

        if (not exists('HPV81.fa')):
            self.log.info("Generating _query")
            self.generate_entrez_query()
            
            self.log.info("Retrieving Ids")
            self.retrieve_ids()
            
            self.log.info("Retrieving Corrfinates for gene")
            self.retrieve_data()
        
        if (not exists('HPV-aligned.fa')):    
            self.log.info("Removing duplicates")
            self.remove_duplicates()
            
            self.log.info("Saving all to common file")
            self.save_common_unaligned()
            
            self.log.info("Aligning")
            self.align()
        else:
            self.log.info("Using cached alignment")

        # FIXME: debug mode
        self.brute_force()
        
        difference = round(time.time() - startTime, 3)            
        self.log.info("Finished: " + str(difference) + " seconds") 
        
        return
        
        self.log.info("Loading sequences for analysis")
        self.load_aligned()
        self.show_sequences()

        self.log.info("Calculating similarity")
        self.update_conservation()
        
        self.log.info("Candidates")
        self.update_candidates()
        
        self.show_gui()
        

    
    
    # Methods
    def generate_entrez_query(self):
        """ Saves and returns _query to get information about all HPV types.
        
            FIXME: rezult example:
            "Human papillomavirus type 18"[Organism] AND L1[Gene Name] AND 7000:8000[Sequence Length]
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
            self._run(self.cdhit + ' -i ' + input_file + ' -o ' + \
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
        handle.write(stdout)
        handle.close()
        if (stderr):
            print stderr    
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

    def update_conservation(self):
        """ Caclulates similarity each sequence symbol with safe and harmful.
            TODO: only harmfull
        """
        # Harful / safe counts
        n_y = len(self._sequences)
        n_x = len(self._sequences[0])
        self._counts['safe'] = 0
        self._counts['harmful'] = 0
        for y in range(0, n_y):
            hpv_type = self._get_type(self._sequences[y].description)
            if (hpv_type in self.safe):
                self._counts['safe'] += 1
            else:
                self._counts['harmful'] += 1
        
        # Loading from catch
        cache = 'conservations.pkl' 
        if (exists(cache)):
            self.log.info("Loading conservations from cache")
            pkl_file = open(cache, 'rb')
            self._similarity = pickle.load(pkl_file)
            pkl_file.close()
            return;
        
        # Default vaules
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
                            hpv_type = self._get_type(self._sequences[y2]
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
        
        # Caching
        self.log.info("Saving conservations to cache")
        output = open(cache, 'wb')
        pickle.dump(self._similarity, output)
        output.close()
        
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
                    
    def show_gui(self):
        """ Showing conservation in GUI
        """
        self.log.info("Launching GUI")
        gui = PyApp()
        gui.set_sequences(self._sequences)
        gui.set_conservation(self._similarity)
        gui.set_count(self._counts['harmful'])
        gui.set_candidates(self._candidates)
        gui.run();
                
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
        
    def update_candidates(self):
        """ Returning best candidates for probes
            FIXME: not finished
        """
        n_harmful = self._counts['harmful']
        n_x = len(self._sequences[0])
        ramp_harmful = n_harmful - 7
        ramp_safe = 5
        self._candidates = []
        for y in range(0, n_harmful):
            print str(y) + ': '
            last_good = 0
            gaps_passed = 0
            for x in range(0, n_x):
                gap = (self._sequences[y,x] == '-')
                if (not gap): 
                    harmful = self._similarity['harmful'][x][y] >= ramp_harmful
                    not_safe = self._similarity['safe'][x][y] <= ramp_safe
                    if (harmful and not_safe):
                        last_good += 1
                    elif (last_good > 25):
                        if (last_good > 35):
                            last_good = 35
                        print '[',
                        for i in range(x - last_good - gaps_passed, x):
                            sys.stdout.write(self._sequences[y,i])
                        print ']',
                        candidate = [y, x - last_good - gaps_passed, x]
                        self._candidates.append(candidate)
                        last_good = 0
                        gaps_passed = 0
                    else:
                        last_good = 0
                        gaps_passed = 0
                else:
                    gaps_passed += 1
            print ''
            
    def brute_force(self):
        self._load_unaligned()
#        self._print_unaligned()

    def _load_unaligned(self):
        self._unaligned = list(SeqIO.parse('HPV-all.fa', "fasta"))
        self._update_types()
        self.log.info("Unaligned:")
        self._print_unaligned()
        self._results = []
        self._size = 15
        self._errors = 2
        self.log.info("Rcursion")
        self._n_harful = 7
        for i in range(0, 7):
            self.log.info("==== %d =====" % i);
            self._first = i
            self._results = []
            self._brute_force_recurce()
            self.log.info("Results:")
            self._show_results()
            

    def _print_unaligned(self):
        for i, sequence in enumerate(self._unaligned):
            hpv_type = self._types[i]
            harm = hpv_type in self.harmfull
            harm2 = i <= self._n_harful
            assert (harm == harm2)

    def _update_types(self):
        self._types = []
        self._n_harful = None
        for i, sequence in enumerate(self._unaligned):
            hpv_type = self._get_type(sequence.description)
            print hpv_type, sequence.description
            if (hpv_type in self.safe and self._n_harful == None):
                self._n_harful = i - 1
            self._types.append(hpv_type)

    def _test_data(self):
        self._unaligned = ["ABCDEFGHJKL",
                           "---aBcDEFGHJKL",
                           "EFGHJKL"
                           ]
            
    def _brute_force_recurce(self, x1 = 0, y2 = -1):
        """ Searching similarities amond all sequences
        """
        size = self._size
        errors = self._errors
        original = self._unaligned[self._first].seq
        current = self._unaligned[y2].seq
        parts = 10
       
        if (y2 == -1):
            # First sequence
            n_x = len(original) - size + 1
            for x1 in range(0, n_x):
                if (x1 % round(n_x / float(parts)) == 0):
                    self.log.info("Comparing %d%%" % (x1 / float(n_x) * 100))
                if (self._first == 0):
                    start = 1
                else:
                    start = 0
                self._brute_force_recurce(x1, start)
        else:
            # Following sequences
            for i in range(0, len(current) - size + 1):
                if (self._compare(original, x1, current, i, size, errors)):
                    self._results.append([x1, y2, i])
                    if (y2 + 1 == self._first):
                        start = y2 + 2
                    else:
                        start = y2 + 1 
                    if (start < len(self._unaligned)):
                        self._brute_force_recurce(x1, start)
            
            
    def _compare(self, list1, x1, list2, x2, length, errors):
        """ Comapre 2 list with allowed error rate
        """
        if (len(list1) < x1 + length or len(list2) < x2 + length):
            return False
        difference = 0
        for i in range(0, length):
            if (list1[x1 + i] != list2[x2 + i]):
                difference += 1
            if (difference > errors):
                return False
        return True
    
    def _show_results(self):
        print len(self._results)
        for result in self._results:
            x1, y2, x2 = result
            source = self._unaligned[self._first].seq[x1:x1 + self._size]
            destination = self._unaligned[y2].seq[x2:x2 + self._size]
            if (y2 == len(self._unaligned)-1):
                arrow = '->  '
                print "%s(%d, %d) (%d, %d) %s %s" % (arrow, self._first, x1, y2,
                                                 x2, source, destination)
            else:
                arrow = '\t'
            


import gtk
import cairo
class PyApp(gtk.Window):
    """ Sequence alignemt visualization
    """
    
    sequences = None
    cr = None

    def __init__(self):
        super(PyApp, self).__init__()
        
        self.set_title("GUI")
        self.set_size_request(1500, 300)
        self.set_position(gtk.WIN_POS_CENTER)
        self.connect("destroy", gtk.main_quit)

    def run(self):
        darea = gtk.DrawingArea()
        darea.connect("expose-event", self.expose)
        darea.connect("motion_notify_event", self.motion_notify_event)
        darea.set_events(gtk.gdk.EXPOSURE_MASK \
                        | gtk.gdk.LEAVE_NOTIFY_MASK \
                        | gtk.gdk.POINTER_MOTION_MASK \
                        | gtk.gdk.POINTER_MOTION_HINT_MASK)
        school = gtk.ScrolledWindow()
        if (self.sequences):
            width = len(self.sequences[0]) * 10
            height = len(self.sequences) * 10
            darea.set_size_request(width, height)
        school.add_with_viewport(darea)
        self.add(school)

        self.show_all()
        gtk.main()

    def set_sequences(self, sequences):
        self.sequences = sequences
    
    def set_conservation(self, conservation):
        self.conservations = conservation
        
    def set_count(self, harmful):
        self.n_harmful = harmful
                
    def set_candidates(self, candidates):
        self.candidates = candidates
    
    def expose(self, widget, event):
        """ Painting
        """
        cr = widget.window.cairo_create()
        cr.set_font_size(10)
        colors = {'A':[1.0, 0.8, 0.8],
                 'C':[0.8, 1.0, 0.8],
                 'T':[0.8, 0.8, 1.0],
                 'G':[0.9, 0.6, 0.9],
                 '-':[1.0, 1.0, 1.0]
                 }

        n_x = len(self.sequences[0])
        n_y = len(self.sequences)
        for x in range(0, n_x):
            for y in range(0, n_y):
                # Fills by base
                color = colors[self.sequences[y,x].upper()]
                cr.set_source_rgb(color[0], color[1], color[2])
                cr.rectangle(x * 10, y * 10, 9, 9)
                cr.fill()
        
        cr.set_source_rgb(1.0, 0.0, 0.0)
        cr.rectangle(0, self.n_harmful * 10 - 1, n_x * 10, 2)
        cr.fill()

        # Bar by conservation
        cr.set_source_rgb(0.1, 0.0, 0.0)
        n_safe = n_y - self.n_harmful
        for x in range(0, n_x):
            for y in range(0, self.n_harmful):
                harm = self.conservations['harmful'][x][y]
                safe = self.conservations['safe'][x][y]
                gap = (self.sequences[y,x] == '-')
                candidate = (harm / float(self.n_harmful) >
                            safe / float(n_safe))
                if (not gap and candidate):
                    cr.rectangle(x * 10 - 2, y * 10, 10, 10)
        cr.set_line_width(1)
        cr.stroke()
        
        # Candidates
        cr.set_source_rgba(0.0, 0.4, 0.0, 0.3)
        for candidate in self.candidates:
            y, x1, x2 = candidate
            cr.rectangle(x1 * 10 - 1, y * 10, (x2 - x1) * 10, 10)
            cr.rectangle(x1 * 10 - 1, (n_y + 1) * 10, (x2 - x1) * 10, 4)
        cr.fill()
        
        
        # Text by base
        cr.set_source_rgb(0.1, 0.1, 0.1)
        for x in range(0, n_x):
            for y in range(0, n_y):
                cr.move_to(x * 10, y * 10 + 10)
                cr.show_text(str(self.sequences[y,x]))
                

    def motion_notify_event(self,widget, event):
        x, y, state = event.window.get_pointer()
        x = int(x / 10)
        y = int(y / 10)
        n_x = len(self.sequences[0])
        n_y = len(self.sequences)
        if (x < n_x and y < n_y and x >= 0 and y >= 0):
            symbol = self.sequences[y,x]
            harm = self.conservations['harmful'][x][y]
            safe = self.conservations['safe'][x][y]
            gap_harm = self.conservations['gaps_harmful'][x][y]
            gap_safe = self.conservations['gaps_safe'][x][y]
            if (y < self.n_harmful):
                type = "H"
            else:
                type = "S"
            self.set_title("%s %s (%d, %d): %d - %d Harm-safe | Gaps: %d %d" %
                           (type, symbol, x, y, harm, safe, gap_harm, gap_safe))
        else:
            self.set_title("GUI")
            
         
            
    
# Running from command line
if __name__ == '__main__':
    logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)
    hpv_analyser = HpvAnalyser()
    hpv_analyser.main()
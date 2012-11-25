# Bioinformatics 2
#
# @author: Aurelijus Banelis 
from genericpath import exists
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
import Bio.Blast.NCBIWWW
import Bio.Blast.NCBIXML


# Preparing input data
print "Started..."
hsa_sequence = SeqIO.read("data1.fa", "fasta");


# Gathering sequences
sequences = []
if (exists("cached-unaligned.fa")):
    # Loading from cache
    sequences = SeqIO.parse("cached-unaligned.fa", "fasta");
else:
    # Downloading sequences from online database
    print "Getting sequences from Swissport database..."
    ncbi = Bio.Blast.NCBIWWW.qblast(program="blastp" , database="swissprot", sequence=hsa_sequence.seq)
        
    # Leaving sequences with good coverage
    for i, alignment in enumerate(Bio.Blast.NCBIXML.read(ncbi).alignments):
        hsp = alignment.hsps[0]        
        #FIXME: not same as in http://www.ncbi.nlm.nih.gov/blast/Blast.cgi
        coverage = float(hsp.query_end-hsp.query_start+1)/float(len(hsa_sequence))*100
        if (coverage > 80):
            name = alignment.title.split(' ')[0]
            sequence = Seq(hsp.sbjct)
            record = SeqRecord(id=name, name=name, description=alignment.title, seq=sequence)
            sequences.append(record)
    
    # Caching results
    print "Saving sequences to cache..."
    cache_handle = open("cached-unaligned.fa", "w")
    SeqIO.write(sequences, cache_handle, "fasta")
    cache_handle.close()


# Printing sequences
for i, sequence in enumerate(sequences):
    print "______%d_________" % i
    print sequence.description
    print sequence.seq
    
    
    
    
    
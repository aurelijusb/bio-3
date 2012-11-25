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
    sequences = list(SeqIO.parse("cached-unaligned.fa", "fasta"));
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


# Printing unaligned sequences
print "Unaligned."
for i, sequence in enumerate(sequences):
    print "[%d] Unaligned %s\n\t\t\t\t\t\t\t%s" % ((i+1), sequence.id, sequence.seq)
    
    
#    
## Aligning
#if (not exists("cached-aligned.fa")):
#    print "Aligning..."
#    from Bio.Align.Applications import MafftCommandline
#    import os
#    os.environ['MAFFT_BINARIES'] = "/usr/lib/mafft/lib/mafft"
#    mafft_exe = "/usr/bin/mafft --localpair --maxiterate 1000"
#    in_file = "cached-unaligned.fa"
#    mafft_cline = MafftCommandline(mafft_exe, input=in_file)
#    stdout, stderr = mafft_cline()
#    handle = open("cached-aligned.fa", "w")
#    handle.write(stdout)
#    handle.close()
#
## Printing aligned sequences
#print "Aligned."
#from Bio import AlignIO
#multiple_seq_alignment = AlignIO.read("cached-aligned.fa", "fasta")
#for i, record in enumerate(multiple_seq_alignment):
#    print "[%d] Aligned %s\n\t\t\t\t\t\t\t%s" % ((i+1), record.id, record.seq)
    

# Calculating conservation:
conservations = []
n_y = len(sequences)
n_x = len(hsa_sequence)
for x in range(0, n_x):
    symbol = sequences[0][x]
    same = 0
    if (symbol != '-'):
        for y in range(1, n_y):
            if (len(sequences[y].seq) > x and sequences[y][x] == symbol):
                same += 1
        conservations.append(same)
    else:
        conservations.append(0)
    
# Printing conservations:
print "Conservation (index, symbol, row, block):"
for x in range(0, n_x):
    print "{0:3}".format(x),
print ''
for x in range(0, n_x):
    print "  %c" % hsa_sequence[x],
print ''
for x in range(0, n_x):
    print "{0:3}".format(conservations[x]),
print ''
    
# Calculating conservatio sums
conservation_sums = []
size = 15
current_sum = 0
for x in range(0, size - 1):
    current_sum += conservations[x]
conservation_sums.append(current_sum)

index_min = 0
min_value = conservation_sums[0]
index_max = 0
max_value = conservation_sums[0]

for x in range(size, n_x):
    current_sum += conservations[x]
    current_sum -= conservations[x - size]
    if (current_sum > max_value):
        max_value = current_sum
        index_max = x
    elif (current_sum < min_value):
        min_value = current_sum
        index_min = x
    conservation_sums.append(current_sum)

for x in range(0, n_x - size):
    print "{0:3}".format(conservation_sums[x]),
print ''

print "\nMin: [%d]=%d | Max: [%d]=%d" % (index_min, min_value, index_max, max_value)
import sys
print '\nMin: ',
for x in range(index_min, index_min + size):
    sys.stdout.write(hsa_sequence[x])
print '\nMax: ',
for x in range(index_max, index_max + size):
    sys.stdout.write(hsa_sequence[x])
print ''

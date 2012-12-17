# Bioinformatics 2
#
# @author: Aurelijus Banelis 
from genericpath import exists
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
import Bio.Blast.NCBIWWW
import Bio.Blast.NCBIXML

# Settings
fragment_size = 15



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
        coverage = float(hsp.query_end-hsp.query_start+1)/float(len(hsa_sequence))*100
        if (coverage > 80):
            name = alignment.title.split(' ')[0]
            sequence = Seq(hsp.sbjct)
            record = SeqRecord(id=name, name=name, description=alignment.title, seq=sequence)
            sequences.append(record)
    
    # Caching results
    print "Saving unaligned sequences to cache..."
    cache_handle = open("cached-unaligned.fa", "w")
    SeqIO.write(sequences, cache_handle, "fasta")
    cache_handle.close()


# Printing unaligned sequences
print "Unaligned sequences:"
for i, sequence in enumerate(sequences):
    print "[%d] Unaligned %s\n\t\t\t\t\t\t\t%s" % ((i+1), sequence.id, sequence.seq)
    
    
# Aligning
if (not exists("cached-aligned.fa")):
    print "Aligning..."
    try:
        from Bio.Align.Applications import MuscleCommandline
        align = MuscleCommandline(input="cached-unaligned.fa", out="cached-aligned.fa")
        align()
    except ImportError:
        print "[Error] Muscle not installed. "    

# Reading cached alignment
from Bio import AlignIO
cache_handle = open("cached-aligned.fa", "r")
multipleSeqAlignment = AlignIO.read(cache_handle, "fasta")
cache_handle.close()
alignmened = multipleSeqAlignment.get_all_seqs()

# Printing aligned sequences
print "Alaligned sequences:"
for i, sequence in enumerate(alignmened):
    print "[%d] Aligned %s\n\t\t\t\t\t\t\t%s" % ((i+1), sequence.id, sequence.seq)

       
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
print "\nSimilarity (index, symbol, column similarity, block similarity):"
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
current_sum = 0
for x in range(0, fragment_size):
    current_sum += conservations[x]
conservation_sums.append(current_sum)
index_min = 0
min_value = conservation_sums[0]
index_max = 0
max_value = conservation_sums[0]

for x in range(fragment_size, n_x):
    current_sum += conservations[x]
    current_sum -= conservations[x - fragment_size]
    if (current_sum > max_value):
        max_value = current_sum
        index_max = x - fragment_size + 1
    elif (current_sum < min_value):
        min_value = current_sum
        index_min = x - fragment_size + 1
    conservation_sums.append(current_sum)

for x in range(0, n_x - fragment_size):
    print "{0:3}".format(conservation_sums[x]),
    
print ''

print "\nFormat: [index]=value | Least expected: [%d]=%d | Most expected: [%d]=%d" % (index_min, min_value, index_max, max_value)
import sys
print '\nLeast expected sequence: ',
for x in range(index_min, index_min + fragment_size):
    sys.stdout.write(hsa_sequence[x])
print '\nMost expected sequence: ',
for x in range(index_max, index_max + fragment_size):
    sys.stdout.write(hsa_sequence[x])
print ''

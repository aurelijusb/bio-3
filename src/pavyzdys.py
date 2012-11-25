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


# Printing unaligned sequences
print "Unaligned."
for i, sequence in enumerate(sequences):
    print "[%d] Unaligned %s\n\t\t\t\t\t\t\t%s" % ((i+1), sequence.id, sequence.seq)
    
    
# Aligning
if (not exists("cached-aligned.fa")):
    print "Aligning..."
    from Bio.Align.Applications import MafftCommandline
    import os
    os.environ['MAFFT_BINARIES'] = "/usr/lib/mafft/lib/mafft"
    mafft_exe = "/usr/bin/mafft --localpair --maxiterate 1000"
    in_file = "cached-unaligned.fa"
    mafft_cline = MafftCommandline(mafft_exe, input=in_file)
    stdout, stderr = mafft_cline()
    handle = open("cached-aligned.fa", "w")
    handle.write(stdout)
    handle.close()

# Printing aligned sequences
print "Aligned."
from Bio import AlignIO
multiple_seq_alignment = AlignIO.read("cached-aligned.fa", "fasta")
for i, record in enumerate(multiple_seq_alignment):
    print "[%d] Aligned %s\n\t\t\t\t\t\t\t%s" % ((i+1), record.id, record.seq)
    
multiple_seq_alignment    

# Calculating conservation:
conservations = []
n_y = len(multiple_seq_alignment)
n_x = multiple_seq_alignment.get_alignment_length()
for x in range(0, n_x):
    symbol = multiple_seq_alignment[0,x]
    same = 0
    if (symbol != '-'):
        for y in range(1, n_y):
            if (multiple_seq_alignment[y,x] == symbol):
                same += 1
        conservations.append(same)
    else:
        conservations.append(0)
    
# Printing conservations:
print "Conservation:"
for x in range(0, n_x):
    print "  %c" % multiple_seq_alignment[0, x],
print ''
for x in range(0, n_x):
    print "{0:3}".format(conservations[x]),
print ''
    
# Calculating conservatio sums
conservation_sums = []
size = 15
sum = 0
for x in range(0, size - 1):
    sum += conservations[x]
conservation_sums.append(sum)

index_min = 0
min = conservation_sums[0]
index_max = 0
max = conservation_sums[0]

for x in range(size, n_x):
    sum += conservations[x]
    sum -= conservations[x - size]
    if (sum > max):
        max = sum
        index_max = x
    elif (sum < min):
        min = sum
        index_min = x
    conservation_sums.append(sum)

for x in range(0, n_x - size):
    print "{0:3}".format(conservation_sums[x]),
print ''

print "Min: [%d]=%d | Max: [%d]=%d" % (index_min, min, index_max, max)
print 'Min:',
for x in range(index_min, index_min + size):
    print multiple_seq_alignment[0,x],
print '\nMax:',
for x in range(index_max, index_max + size):
    print multiple_seq_alignment[0,x],
print ''
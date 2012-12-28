# Bioinformatics 3
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

###########################################
# Databases:
#    nucleotide
#    chromosome
#    gene
#    probe
#
# http://www.ncbi.nlm.nih.gov/sites/gquery
#
#
#
###########################################




# Preparing input data
print "Started..."
hsa_sequence = SeqIO.read("data3.fa", "fasta");

from Bio import Entrez
Entrez.email = 'aurelijus@banelis.lt'
handle = Entrez.egquery(term="Human papillomavirus type")
record = Entrez.read(handle)

print "-------------"

#for row in record["eGQueryResult"]:
#    for key in row.keys():
#        print "%s=%s" % (key, row[key])

for row in record["eGQueryResult"]:
    if (row['Status'] == 'Ok' and row['Count'] > 0):
        print row['Count'], row['MenuName'], row['DbName'] 
    


exit(0)





from Bio import Entrez
Entrez.email = 'aurelijus@banelis.lt'
handle = Entrez.esearch(db="nucleotide",term="Human papillomavirus")
record = Entrez.read(handle)

for key in record.keys():
    print "%s=%s" % (key, record[key])
    
    
#    
#    
    
##
## Gathering sequences
#sequences = []
#if (exists("cached-unaligned.fa")):
#    # Loading from cache
#    sequences = list(SeqIO.parse("cached-unaligned.fa", "fasta"));
#else:
#    # Downloading sequences from online database
#    print "Getting sequences from Swissport database..."
#    ncbi = Bio.Blast.NCBIWWW.qblast(program="blastp" , database="gene", sequence=hsa_sequence.seq)
#        
#    # Leaving sequences with good coverage
#    for i, alignment in enumerate(Bio.Blast.NCBIXML.read(ncbi).alignments):
#        hsp = alignment.hsps[0]
#        coverage = float(hsp.query_end-hsp.query_start+1)/float(len(hsa_sequence))*100
#        if (coverage > 80):
#            name = alignment.title.split(' ')[0]
#            sequence = Seq(hsp.sbjct)
#            record = SeqRecord(id=name, name=name, description=alignment.title, seq=sequence)
#            sequences.append(record)
#    
#    # Caching results
#    print "Saving unaligned sequences to cache..."
#    cache_handle = open("cached-unaligned.fa", "w")
#    SeqIO.write(sequences, cache_handle, "fasta")
#    cache_handle.close()
#

#
## Printing unaligned sequences
#print "Unaligned sequences:"
#for i, sequence in enumerate(sequences):
#    print "[%d] Unaligned %s\n\t\t\t\t\t\t\t%s" % ((i+1), sequence.id, sequence.seq)
#    
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
#
#       
## Calculating conservation:
#conservations = []
#n_y = len(sequences)
#n_x = len(hsa_sequence)
#for x in range(0, n_x):
#    symbol = sequences[0][x]
#    same = 0
#    if (symbol != '-'):
#        for y in range(1, n_y):
#            if (len(sequences[y].seq) > x and sequences[y][x] == symbol):
#                same += 1
#        conservations.append(same)
#    else:
#        conservations.append(0)
#    
## Printing conservations:
#print "\nSimilarity (index, symbol, column similarity, block similarity):"
#for x in range(0, n_x):
#    print "{0:3}".format(x),
#print ''
#for x in range(0, n_x):
#    print "  %c" % hsa_sequence[x],
#print ''
#for x in range(0, n_x):
#    print "{0:3}".format(conservations[x]),
#print ''
#    
## Calculating conservatio sums
#conservation_sums = []
#current_sum = 0
#for x in range(0, fragment_size):
#    current_sum += conservations[x]
#conservation_sums.append(current_sum)
#index_min = 0
#min_value = conservation_sums[0]
#index_max = 0
#max_value = conservation_sums[0]
#
#for x in range(fragment_size, n_x):
#    current_sum += conservations[x]
#    current_sum -= conservations[x - fragment_size]
#    if (current_sum > max_value):
#        max_value = current_sum
#        index_max = x - fragment_size + 1
#    elif (current_sum < min_value):
#        min_value = current_sum
#        index_min = x - fragment_size + 1
#    conservation_sums.append(current_sum)
#
#for x in range(0, n_x - fragment_size):
#    print "{0:3}".format(conservation_sums[x]),
#    
#print ''
#
#print "\nFormat: [index]=value | Least expected: [%d]=%d | Most expected: [%d]=%d" % (index_min, min_value, index_max, max_value)
#import sys
#print '\nLeast expected sequence: ',
#for x in range(index_min, index_min + fragment_size):
#    sys.stdout.write(hsa_sequence[x])
#print '\nMost expected sequence: ',
#for x in range(index_max, index_max + fragment_size):
#    sys.stdout.write(hsa_sequence[x])
#print ''

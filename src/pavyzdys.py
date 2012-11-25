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
    
# Phylogenetic trees
if (not exists("cached-aligned.aln")):
    print "Prepairing trees (1/3)..."
    from Bio.Align.Applications import MuscleCommandline
    cmdline = MuscleCommandline(input="cached-aligned.fa", out="cached-aligned.aln", clw=True)
    cmdline()
if (not exists("cached-aligned.phy")):
    print "Prepairing trees (2/3)..."
    AlignIO.convert("cached-aligned.aln", "clustal", "cached-aligned.phy", "phylip-relaxed")

if (not exists("cached-aligned.phy_phyml_tree.txt")):
    print "Prepairing trees (3/3)..."
    from Bio.Phylo.Applications import PhymlCommandline
    #TODO: optimization
    cmdline = PhymlCommandline(input='cached-aligned.phy', datatype='aa', model='LG')
    out_log, err_log = cmdline()

print "Painting phylogenetic tree..."
from Bio import Phylo
egfr_tree = Phylo.read("cached-aligned.phy_phyml_tree.txt", "newick")
Phylo.draw_ascii(egfr_tree)

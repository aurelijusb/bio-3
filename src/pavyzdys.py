#Sis skrtipukas vykdo NCBI blast uzklausa ir atspausdina rastus atitikmenis bei  sekas
from genericpath import exists
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


#enterez_query = '"serum albumin"[Protein name] AND mammals[Organism]'
#ncbi = Bio.Blast.NCBIWWW.qblast(program="blastp" , database="swissprot", sequence="4502027", entrez_query=enterez_query,  hitlist_size=500,  expect=100.0) #seka, pagal kuria ieskoma  galima nurodyti ir kodu

# Nuskaitom XML formatu parsiustus duomenis. Patartina siuos duomenis issisaugoti i faila, ir pakartotinai kreiptis tik esant butinybei.
#blast = Bio.Blast.NCBIXML.read(ncbi)
#for sequence in blast.alignments:
#    print '>%s'%sequence.title # rasto atitikmens pavadinimas fasta formatu
#    print sequence.hsps[0].sbjct # rasto atitikmens labiausiai patikimo sutampancio fragmento seka. Kiti maziau patikimi fragmentai [kuriu indeksai didesni nei 0] nedomina. 

print "Started..."


# Reading input data
from Bio import SeqIO
hsa_sequence = SeqIO.read("data1.fa", "fasta");

#Qriteria for reliable results
import math
size = len(hsa_sequence);
percent = 80
slen = "%d:%d[slen]" % (round(size * percent / 100), size)
print slen

#from Bio.Blast import NCBIWWW
#print "Gathering..."
#result_handle = NCBIWWW.qblast("blastn", "nt", "8332116")
#
#save_file = open("my_blast.xml", "w")
#save_file.write(result_handle.read())
#save_file.close()
#result_handle.close()

#Searching
sequences = []
if (not exists("sequence.fasta")):
    sequences = SeqIO.parse("cached-unaligned.fa", "fasta");
else:
    print "Getting sequences from Swissport database..."
    import Bio.Blast.NCBIWWW
    import Bio.Blast.NCBIXML
    ncbi = Bio.Blast.NCBIWWW.qblast(program="blastp" , database="swissprot", sequence=hsa_sequence.seq, entrez_query=slen)
    i = 0
    for sequence in Bio.Blast.NCBIXML.read(ncbi).alignments:
        i += 1
        name = sequence.title.split(' ')[0]
        sequences.append(SeqRecord(id=name, name=name, description=sequence.title, seq=Seq(sequence.hsps[0].sbjct)))
        if (i > 29):
            for j, hh in enumerate(sequence.hsps):
                print ">%d>%s" % (j, hh.sbjct)
    
    # Caching Fasta
    output_handle = open("example.fasta", "w")
    SeqIO.write(sequences, output_handle, "fasta")
    output_handle.close()
    
    # Caching XML
    ncbi.seek(0, 0)
    save_file = open("my_blast.xml", "w")
    save_file.write(ncbi.read())
    save_file.close()
    ncbi.close()

print "fiished"

#    sequences = blast.alignments
#      
for i, sequence in enumerate(sequences):
    print "______%d_________" % i
    print sequence
    
    
    
    
    
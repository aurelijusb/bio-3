f = open('../data/duomenys-aligned.fa', 'r')
# Reading FASTA file
names = []
sequences = []
index = -1;
for line in f:
    if (line[0] == '>'):
        names.append(line[1:].strip())
        sequences.append('')
        index += 1 
    else:
        sequences[index] += line.strip()
f.close()
        
# Sequences
for i, sequence in enumerate(sequences):
    print "%s > %s" % (sequence, names[i])
    
print "_____________________________________________"
    
# Remove gaps
mainSequence = sequences[0];
n = len(mainSequence)
i = 0
while (i < n):
    if (mainSequence[i] == '-'):
        mainSequence = mainSequence[0:i] + mainSequence[i+1:]
        n -= 1 
        for j, sequence in enumerate(sequences):
            sequences[j] = sequences[j][0:i] + sequences[j][i+1:]
    else:
        i += 1

print "_____________________________________________"
        
# Sequences
for i, sequence in enumerate(sequences):
    print "%s > %s" % (sequence, names[i])
    
f = open('../data/duomenys-trimmed.fa', 'w')
# TO FASTA
for i, sequence in enumerate(sequences):
    f.write(">%s\n" % names[i])
    for j in range(0, len(sequence)-1, 72):
        f.write("%s\n" % sequence[j:j+72])
f.close()
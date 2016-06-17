mapFI = open('/home/fs01/yw595/Entrez_and_IPI_unique.csv')
IPIToEntrez = {}
for line in mapFI:
    line = line.strip()
    words = line.split('\t')
    IPIToEntrez[words[1]] = words[0]
mapFI.close()

protFI = open('/home/fs01/yw595/UACC_257_protexp.csv')
entrezToExp = {}
for line in protFI:
    line = line.strip()
    words = line.split(',')
    if words[0] != 'NA' and words[0] in IPIToEntrez:
        entrezToExp[IPIToEntrez[words[0]]] = pow(2,float(words[1]))
protFI.close()

#print(IPIToEntrez)

outFI = open('/home/fs01/yw595/MATLAB/FALCON/input/NCI60Sims/nci60mRNA/UACC_257_test.csv','w')
outFI.write('gene\tmean\tvar\n')
for entrez in entrezToExp:
    outFI.write(entrez + '\t' + str(entrezToExp[entrez]) + '\t' + '1\n')
outFI.close()

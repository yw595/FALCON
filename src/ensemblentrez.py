inputDir = '/mnt/vdb/home/ubuntu2/MATLAB/FALCON/input/'
inFI = open(inputDir+'ensemblentrez.txt')
outFI = open(inputDir+'ensemblentrez2.txt','w')
genemap = {}
currentkey = ''
for line in inFI:
    geneid = line[line.index('>')+1:line.index('<')]
    if geneid!='':
        if geneid.startswith('ENSG'):
            genemap[geneid] = []
            currentkey = geneid
        else:
            genemap[currentkey].append(geneid)
inFI.close()
for currentkey in genemap:
    if len(genemap[currentkey])>1:
        outFI.write(currentkey+'\t')
        for i in range(len(genemap[currentkey])-1):
            outFI.write(genemap[currentkey][i]+'\t')
        outFI.write(genemap[currentkey][len(genemap[currentkey])-1]+'\n')
outFI.close()

use1051 = False
useGTEx = True
if use1051:
    rnaFI = open(inputDir+'rnaseq1051.txt')
    line = rnaFI.readline().strip()
    allIDs = line.split('\t')[4:]
    edgeIDs = line.split('\t')[4:13]
    coreIDs = line.split('\t')[13:]
    expMap = {}
    for i in range(len(edgeIDs)):
        expMap[edgeIDs[i]] = {}
    for i in range(len(coreIDs)):
        expMap[coreIDs[i]] = {}
    for line in rnaFI:
        words = line.strip().split('\t')
        for i in range(4,len(words)):
            if words[0] in genemap and len(genemap[words[0]]) > 1:
                entrezID = genemap[words[0]][1]
                if entrezID not in expMap[allIDs[i-4]]:
                    expMap[allIDs[i-4]][entrezID] = 0
                if words[i]=='':
                    words[i]='0'
                expMap[allIDs[i-4]][entrezID] = expMap[allIDs[i-4]][entrezID] + float(words[i])
    rnaFI.close()

    for expID in expMap:
        outFI = open(inputDir+'/rnaseq1051/'+expID+'.txt','w')
        outFI.write('gene\tmean\tvar\n')
        for geneid in expMap[expID]:
            outFI.write(geneid+'\t'+str(expMap[expID][geneid])+'\t'+'1.0'+'\n')
        outFI.close()

if useGTEx:
    firstly = False
    rnaFI = open(inputDir+'gtex/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_reads.gct')
    rnaFI.readline()
    rnaFI.readline()
    line = rnaFI.readline().strip()
    allIDs = line.split('\t')[2:]
    expMap = {}
    for i in range(len(allIDs)):
        expMap[allIDs[i]] = {}
    linecount = 0
    for line in rnaFI:
        linecount += 1
        if (firstly and linecount < 25000) or (not firstly and linecount>=25000):
            print(linecount)
            words = line.strip().split('\t')
            for i in range(2,len(words)):
                genename = words[0]
                if genename.index('.')!=-1:
                    genename = genename[:genename.index('.')]
                if genename in genemap and len(genemap[genename]) > 1:
                    entrezID = genemap[genename][1]
                    if entrezID not in expMap[allIDs[i-2]]:
                        expMap[allIDs[i-2]][entrezID] = 0
                    if words[i]=='':
                        words[i]='0'
                    expMap[allIDs[i-2]][entrezID] = expMap[allIDs[i-2]][entrezID] + float(words[i])
    rnaFI.close()

    for expID in expMap:
        if firstly:
            outFI = open(inputDir+'/gtex/'+expID+'.txt','w')
        else:
            outFI = open(inputDir+'/gtex/'+expID+'.txt','a')
        if firstly:
            outFI.write('gene\tmean\tvar\n')
        for geneid in expMap[expID]:
            outFI.write(geneid+'\t'+str(expMap[expID][geneid])+'\t'+'1.0'+'\n')
        outFI.close()

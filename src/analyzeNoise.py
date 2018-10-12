import os
import scipy.stats
import numpy as np
import sys
sys.path.append('/mnt/vdb/home/ubuntu2/MATLAB/SEED/src')
from bootstrapAGORAIBD import writeData

inDir = '/mnt/vdb/home/ubuntu2/MATLAB/FALCON/output/gbmscript/gtex/'
ages = ['young','old']
fullmap = {}
for age in ages:
    inDir1 = inDir+age+'/'
    files = os.listdir(inDir1)
    agemap = {}
    for afile in files:
        prefix = afile[:afile.index('_')]
        inFI = open(inDir1+afile)
        inFI.readline()
        if afile.endswith('noise1'):
            rxnstoexprratios = {}
            rxnstoflux1 = {}
        else:
            rxnstoflux0 = {}
        lineCount = 0
        for line in inFI:
            lineCount += 1
            words = line.strip().split('\t')
            if words[1]!='NaN' and words[2]!='NaN' and float(words[2])!=0:
                if afile.endswith('noise1'):
                    exprratio = abs(float(words[1])-float(words[2]))/float(words[2])
                    rxnstoexprratios[lineCount] = exprratio
                    rxnstoflux1[lineCount] = float(words[0])
                else:
                    rxnstoflux0[lineCount] = float(words[0])
        inFI.close()
        if prefix not in agemap:
            agemap[prefix] = {}
        if afile.endswith('noise1'):
            agemap[prefix]['exprratio'] = rxnstoexprratios
            agemap[prefix]['flux1'] = rxnstoflux1
        else:
            agemap[prefix]['flux0'] = rxnstoflux0
    fullmap[age] = agemap

#nonsense = nonsense+1
wilcarrs = {}
printarr1 = []
printarr2 = []
printarr3 = []
printarr4 = []
stdcutarrs = {}
for age in ages:
    agemap = fullmap[age]
    wilcarr = []
    stdcutarr = []
    print(age)
    for prefix in agemap:
        prefixkeys = agemap[prefix]['exprratio'].keys()
        prefixexprratios = []
        prefixfluxratios = []
        prefixcompratios = []
        for i in range(len(prefixkeys)):
            if prefixkeys[i] in agemap[prefix]['flux0'] and prefixkeys[i] in agemap[prefix]['flux1']:
                flux0 = agemap[prefix]['flux0'][prefixkeys[i]]
                flux1 = agemap[prefix]['flux1'][prefixkeys[i]]
                exprratio = agemap[prefix]['exprratio'][prefixkeys[i]]
                if flux0!=0:
                    fluxratio = abs(flux1-flux0)/flux0
                    prefixfluxratios.append(fluxratio)
                    prefixexprratios.append(exprratio)
                    if exprratio!=0:
                        prefixcompratios.append(fluxratio/exprratio)
                        wilcarr.append(fluxratio/exprratio)
                        printarr1.append(fluxratio/exprratio)
                        printarr2.append(age+' '+prefix)
                        printarr3.append(age)
                        printarr4.append(prefix)
        [rho, pval] = scipy.stats.spearmanr(prefixexprratios,prefixfluxratios)
        print(prefix)
        print(rho)
        print(pval)
        print(np.median(prefixcompratios))
        print(max(prefixcompratios))
        print(min(prefixcompratios))
        print(np.std(prefixcompratios))
        prefixcompratios = np.array(prefixcompratios)
        stdcutarrtemp = prefixcompratios[abs(prefixcompratios-np.mean(prefixcompratios)) < 2*np.std(prefixcompratios)]
        stdcutarr.extend(stdcutarrtemp)
    wilcarrs[age] = wilcarr
    stdcutarrs[age] = stdcutarr
        #nonsense = nonsense+1
#print(wilcarrs['young'])
#print(wilcarrs['old'])
[stat, pval] = scipy.stats.mannwhitneyu(wilcarrs['young'],wilcarrs['old'])
print(stat)
print(pval)
print(np.mean(wilcarrs['young'])-np.mean(wilcarrs['old']))
[stat, pval] = scipy.stats.mannwhitneyu(stdcutarrs['young'],stdcutarrs['old'])
print(stat)
print(pval)
print(np.mean(wilcarrs['young'])-np.mean(wilcarrs['old']))
print(np.std(wilcarrs['young']))
print(np.std(wilcarrs['old']))
writeData([printarr1,printarr2],'/mnt/vdb/home/ubuntu2/compareNoise.txt',delimiter='\t',headers=['ratio','label'])
writeData([printarr1,printarr3,printarr4],'/mnt/vdb/home/ubuntu2/compareNoiseANOVA.txt',delimiter='\t',headers=['ratio','age','sample'])

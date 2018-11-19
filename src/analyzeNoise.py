import os
import scipy.stats
from scipy.stats import hypergeom
import numpy as np
import sys
import copy
sys.path.append('/mnt/vdb/home/ubuntu2/MATLAB/SEED/src')
from bootstrapAGORAIBD import writeData
from bootstrapAGORAIBD import segmentFluxBySubsystemGroup
from bootstrapAGORAIBD import sortByIdx
from bootstrapAGORAIBD import addIdxStrings
from bootstrapAGORAIBD import absArr

inDir = '/mnt/vdb/home/ubuntu2/MATLAB/FALCON/output/gbmscript/gtex/'
ages = ['young','old']
checkedprefixes = []
fullmap = {}
samplelim = 3
for age in ages:
    inputDir1 = inDir+age+'/'
    files = os.listdir(inputDir1)
    agemap = {}
    for afile in files:
        testFI = open(inputDir1+afile)
        linecount = 0
        for line in testFI:
            linecount += 1
        testFI.close()
        if afile.find('natural')!=-1 and afile.find('noiseall')!=-1 and linecount>10:
            prefix = afile[:afile.index('_')]
            commonfilepart = afile[:afile.index('natural')-1]
            if prefix not in checkedprefixes:
                checkedprefixes.append(prefix)
                tissue = afile.split('_')[1][6:]
                if tissue not in agemap:
                    agemap[tissue] = {}
                agemap[tissue][prefix] = {}
                if os.path.exists(inputDir1+commonfilepart):
                    inFI = open(inputDir1+commonfilepart)
                    inFI.readline()
                    triplelist = [[],[],[]]
                    for line in inFI:
                        words = line.strip().split('\t')
                        for j in range(len(words)):
                            if words[j]=='NaN':
                                appendval = float('nan')
                            else:
                                appendval = float(words[j])
                            triplelist[j].append(appendval)
                    agemap[tissue][prefix]['nonoise'] = triplelist
                for k in range(1,samplelim+1,1):
                    if os.path.exists(inputDir1+commonfilepart+'_naturalnoiseall'+str(k)):
                        inFI = open(inputDir1+commonfilepart+'_naturalnoiseall'+str(k))
                        inFI.readline()
                        triplelist = [[],[],[]]
                        for line in inFI:
                            words = line.strip().split('\t')
                            for j in range(len(words)):
                                if words[j]=='NaN':
                                    appendval = float('nan')
                                else:
                                    appendval = float(words[j])
                                triplelist[j].append(appendval)
                        agemap[tissue][prefix]['naturalnoise'+str(k)] = triplelist
    fullmap[age] = agemap

tissuescommon = []
for tissue in fullmap['young'].keys():
    if tissue not in tissuescommon:
        tissuescommon.append(tissue)
for tissue in fullmap['old'].keys():
    if tissue not in tissuescommon:
        tissuescommon.append(tissue)
tissuessubstable = []
tissuessubsheader = []
for tissue in tissuescommon:
    totalvar = {}
    totalvar['old'] = {}
    totalvar['young'] = {}
    totalvarrxn = {}
    totalvarrxn['old'] = {}
    totalvarrxn['young'] = {}
    totalvarrxn2 = {}
    totalvarrxn2['old'] = {}
    totalvarrxn2['young'] = {}
    rxnpvals = [1 for i in range(7440)]
    rxnrangediffs = [0 for i in range(7440)]
    for age in fullmap:
        for prefix in fullmap[age][tissue]:
            totalvarrxn2[age][prefix] = [['nan' for nanidx2 in range(samplelim)] for nanidx in range(7440)]
            rxnminarr = ['nan' for i in range(7440)]
            rxnmaxarr = ['nan' for i in range(7440)]
            print(age)
            print(prefix)
            rxnvalarr = []
            for k in range(1,samplelim+1,1):
                print(k)
                for l in range(7440):
                    rxnval = fullmap[age][tissue][prefix]['naturalnoise'+str(k)][0][l]
                    if rxnval!='nan':
                        totalvarrxn2[age][prefix][l][k-1] = rxnval
                    if rxnminarr[l]=='nan' or rxnval < rxnminarr[l]:
                        rxnminarr[l] = rxnval
                    if rxnmaxarr[l]=='nan' or rxnval > rxnmaxarr[l]:
                        rxnmaxarr[l] = rxnval

            totalvartemp = 0
            totalvarrxn[age][prefix] = ['nan' for nanidx in range(7440)]
            for i in range(len(rxnminarr)):
                if rxnminarr[i]!='nan' and rxnmaxarr[i]!='nan':
                    totalvartemp += (rxnmaxarr[i]-rxnminarr[i])
                    totalvarrxn[age][prefix][i] = (rxnmaxarr[i]-rxnminarr[i])
            totalvar[age][prefix] = totalvartemp
    [totalstat, totalpval] = scipy.stats.mannwhitneyu(totalvar['young'].values(),totalvar['old'].values())
    totalvararr = [[],[]]
    testrealall = True
    for i in range(7440):
        rxnfluxlist = [[],[]]
        ages = ['young','old']
        allnum = True
        for j in range(len(ages)):
            for prefix in totalvarrxn[ages[j]]:
                if testrealall:
                    for k in range(samplelim):
                        if totalvarrxn2[ages[j]][prefix][i][k]=='nan':
                            allnum = False
        for j in range(len(ages)):
            for prefix in totalvarrxn[ages[j]]:
                if testrealall:
                    for k in range(samplelim):
                        if not allnum:
                            #print('HERE')
                            rxnfluxlist[j].append(0)
                        else:
                            rxnfluxlist[j].append(totalvarrxn2[ages[j]][prefix][i][k]+1000)
                else:
                    rxnfluxlist[j].append(totalvarrxn[ages[j]][prefix][i]+1000)
        allrxnfluxmean = copy.deepcopy(rxnfluxlist[0])
        allrxnfluxmean.extend(rxnfluxlist[1])
        allrxnfluxmean = np.mean(allrxnfluxmean)
        for j in range(len(ages)):
            for k in range(len(rxnfluxlist[j])):
                rxnfluxlist[j][k] = abs((rxnfluxlist[j][k]-allrxnfluxmean)/allrxnfluxmean)
        if max(rxnfluxlist[0])!=max(rxnfluxlist[1]) and min(rxnfluxlist[0])!=min(rxnfluxlist[1]) and allrxnfluxmean!=0:
            [rxnstat, rxnpval] = scipy.stats.mannwhitneyu(rxnfluxlist[0],rxnfluxlist[1])
            #nonsense = nonsense+1
            rxnpvals[i] = rxnpval
            rxnrangediffs[i] = np.mean(rxnfluxlist[0])-np.mean(rxnfluxlist[1])
    inFI = open('/mnt/vdb/home/ubuntu2/recon2.txt')
    recon2rxns = []
    recon2rxnnames = []
    recon2subs = []
    for line in inFI:
        words = line.strip().split('\t')
        recon2rxns.append(words[0])
        recon2rxnnames.append(words[1])
        recon2subs.append(words[2])
    inFI.close()

    rxnpvals = np.array(rxnpvals)
    recon2subs = np.array(recon2subs)
    rxnrangediffs = np.array(rxnrangediffs)
    uniqSubs = np.unique(recon2subs)
    numrxnspos = []
    numsigrxnspos = []
    hypergeomparrpos = []
    for sub in uniqSubs:
        N = len(recon2rxns)
        M = sum(recon2subs==sub)
        K = sum(np.logical_and(rxnpvals<.05,rxnrangediffs>0))
        x = sum(np.logical_and(np.logical_and(rxnpvals<.05,rxnrangediffs>0),recon2subs==sub))
        numrxnspos.append(M)
        numsigrxnspos.append(x)
        hypergeomparrpos.append(1-hypergeom.cdf(x-1,N,M,K))
    numrxnsneg = []
    numsigrxnsneg = []
    hypergeomparrneg = []
    for sub in uniqSubs:
        N = len(recon2rxns)
        M = sum(recon2subs==sub)
        K = sum(np.logical_and(rxnpvals<.05,rxnrangediffs<0))
        x = sum(np.logical_and(np.logical_and(rxnpvals<.05,rxnrangediffs<0),recon2subs==sub))
        numrxnsneg.append(M)
        numsigrxnsneg.append(x)
        hypergeomparrneg.append(1-hypergeom.cdf(x-1,N,M,K))
    labelarr = []
    valuesarr = []
    youngvals = totalvar['young'].values()
    for i in range(len(youngvals)):
        labelarr.append('young')
        valuesarr.append(youngvals[i])
    oldvals = totalvar['old'].values()
    for i in range(len(oldvals)):
        labelarr.append('old')
        valuesarr.append(oldvals[i])
    writeData([labelarr,valuesarr],'/mnt/vdb/home/ubuntu2/youngOldPopDiff'+tissue+'.txt',delimiter='\t',headers=['label','overallfluxdiff'])
    writeData([recon2rxns,recon2rxnnames,recon2subs,rxnpvals,rxnrangediffs],'/mnt/vdb/home/ubuntu2/youngOldStabilityDiffRxns'+tissue+'.txt',delimiter='\t',headers=['rxn','rxnname','sub','wilcoxon p-val','average flux range diff'])
    if len(tissuessubstable)==0:
        tissuessubstable.append(list(uniqSubs))
        tissuessubsheader.append('tissue')    
    writeData([uniqSubs,numrxnspos,numsigrxnspos,hypergeomparrpos],'/mnt/vdb/home/ubuntu2/youngOldStabilityDiffSubsPos'+tissue+'.txt',delimiter='\t',headers=['sub','num reactions','num sig reactions','hypergeometric p-val'])
    writeData([uniqSubs,numrxnsneg,numsigrxnsneg,hypergeomparrneg],'/mnt/vdb/home/ubuntu2/youngOldStabilityDiffSubsNeg'+tissue+'.txt',delimiter='\t',headers=['sub','num reactions','num sig reactions','hypergeometric p-val'])
    tissuessubstable.append(hypergeomparrpos)
    tissuessubstable.append(hypergeomparrneg)
    tissuessubsheader.append(tissue+' pos')
    tissuessubsheader.append(tissue+' neg')
writetable = []
sumcol = [0 for i in range(len(tissuessubstable[0]))]
for i in range(len(tissuessubstable[0])):
    for j in range(1,len(tissuessubstable)):
        if tissuessubstable[j][i]<.005:
            sumcol[i] += 1
tissuessubstable[0].append('num sig subsystems in tissue')
writetable.append(tissuessubstable[0])
#writetable.append('num sig subsystems in tissue')
for i in range(1,len(tissuessubstable)):
    tissuessubstable[i].append(sum(np.array(tissuessubstable[i])<.005))
    writetable.append(tissuessubstable[i])
writetable.append(sumcol)
tissuessubsheader.append('num sig occurrences in subsystem')
writeData(writetable,'/mnt/vdb/home/ubuntu2/youngOldAllTissues.txt',delimiter='\t',headers=tissuessubsheader)
nonsense = nonsense+1

inDir = '/mnt/vdb/home/ubuntu2/MATLAB/FALCON/output/gbmscript/gtex/'
ages = ['young','old']
fullmap = {}
for age in ages:
    inDir1 = inDir+age+'/'
    files = os.listdir(inDir1)
    agemap = {}
    for afile in files:
        if afile.find('tissue')==-1:
            prefix = afile[:afile.index('_')]
            inFI = open(inDir1+afile)
            inFI.readline()
            if afile.endswith('noise1'):
                rxnstoexprratios = {}
                rxnstoflux1 = {}
            else:
                rxnstoflux0 = {}
                rxnstoexpr0 = {}
            lineCount = 0
            for line in inFI:
                lineCount += 1
                words = line.strip().split('\t')
                if words[1]=='NaN':
                    firstword = words[1]
                else:
                    firstword = float(words[1])
                if words[2]=='NaN':
                    secondword = words[2]
                else:
                    secondword = float(words[2])
                if afile.endswith('noise1'):
                    if firstword!='NaN' and secondword!='NaN' and secondword!=0:
                        exprratio = abs(float(words[1])-float(words[2]))/float(words[2])
                    else:
                        exprratio = 'NaN'
                    rxnstoexprratios[lineCount] = exprratio
                    rxnstoflux1[lineCount] = float(words[0])
                else:
                    rxnstoflux0[lineCount] = float(words[0])
                    rxnstoexpr0[lineCount] = firstword
            inFI.close()
            if prefix not in agemap:
                agemap[prefix] = {}
            if afile.endswith('noise1'):
                agemap[prefix]['exprratio'] = rxnstoexprratios
                agemap[prefix]['flux1'] = rxnstoflux1
            else:
                agemap[prefix]['flux0'] = rxnstoflux0
                agemap[prefix]['expr0'] = rxnstoexpr0
        fullmap[age] = agemap

inFI = open('/mnt/vdb/home/ubuntu2/recon2.txt')
recon2rxns = []
recon2rxnnames = []
recon2subs = []
for line in inFI:
    words = line.strip().split('\t')
    recon2rxns.append(words[0])
    recon2rxnnames.append(words[1])
    recon2subs.append(words[2])
inFI.close()
        
#print(fullmap)
printarr1 = []
printarr2 = []
reshapedmapretain = {}
ratiosarrretain = [[],[]]
youngidxs = []
oldidxs = []
agearr = []
samplenamearr = []
allexprarr = []
offset = 0
for age in ['old','young']:
    agemap = fullmap[age]
    reshapedmap = {}
    for prefix in agemap.keys():
        agearr.append(age)
        samplenamearr.append(prefix)
        allexprarr.append([])
        for lineCount in agemap[prefix]['flux0'].keys():
            if lineCount not in reshapedmap:
                reshapedmap[lineCount] = {}
                reshapedmap[lineCount]['expr'] = []
                reshapedmap[lineCount]['flux'] = []
            reshapedmap[lineCount]['expr'].append(agemap[prefix]['expr0'][lineCount])
            reshapedmap[lineCount]['flux'].append(agemap[prefix]['flux0'][lineCount])
    reshapedmapretain[age] = reshapedmap
    ratioarr = []
    ratiosarrretaintemp = [0 for z in range(7440)]
    for lineCount in reshapedmap.keys():
        fluxarr = reshapedmap[lineCount]['flux']
        exprarr = reshapedmap[lineCount]['expr']
        for i in range(len(exprarr)):
            allexprarr[offset+i].append(exprarr[i])
        nonan = True
        for i in range(len(exprarr)):
            if exprarr[i]=='NaN':
                nonan = False
        if nonan and np.mean(exprarr)!=0:#max(exprarr)-min(exprarr)!=0:
            correctstat = np.std(exprarr)/np.mean(exprarr)#(max(fluxarr)-min(fluxarr))/(max(exprarr)-min(exprarr))
            ratioarr.append(correctstat)
            ratiosarrretaintemp[lineCount-1] = correctstat
            printarr2.append(age)
            printarr1.append(correctstat)
    print(np.median(ratioarr))
    print(np.mean(ratioarr))
    print(np.std(ratioarr))
    if age=='old':
        ratioarrold = ratioarr
        ratiosarrretain[1] = ratiosarrretaintemp
    else:
        ratioarryoung = ratioarr
        ratiosarrretain[0] = ratiosarrretaintemp
    offset += len(samplenamearr)

writedata = [recon2rxns,recon2rxnnames,recon2subs]
writeheaders = ['rxn','rxnname','sub']
for i in range(len(allexprarr)):
    writedata.append(allexprarr[i])
    writeheaders.append(agearr[i]+' '+samplenamearr[i])
writeData(writedata,'/mnt/vdb/home/ubuntu2/youngOldExprData.txt',delimiter='\t',headers=writeheaders)
#[stat, pval] = scipy.stats.mannwhitneyu(ratioarrold, ratioarryoung)
[stat, pval] = scipy.stats.wilcoxon(ratiosarrretain[0], ratiosarrretain[1])
printarr5 = []
printarr6 = []
for i in range(len(ratiosarrretain[0])):
    if ratiosarrretain[0][i]!=0 and ratiosarrretain[1][i]!=0:
        printarr5.append(ratiosarrretain[0][i])
        printarr6.append(ratiosarrretain[1][i])
writeData([printarr5,printarr6],'/mnt/vdb/home/ubuntu2/compareNoiseNaturalScatter.txt',delimiter='\t',headers=['youngratio','oldratio'])
print(stat)
print(pval)
print(len(fullmap['old']))
print(len(fullmap['young']))
writeData([printarr1,printarr2],'/mnt/vdb/home/ubuntu2/compareNoiseNatural.txt',delimiter='\t',headers=['ratio','label'])

inFI = open('/mnt/vdb/home/ubuntu2/pyruvatecidxs.txt')
cidxs = []
for line in inFI:
    cidxs.append(int(line.strip()))
inFI.close()
inFI = open('/mnt/vdb/home/ubuntu2/pyruvatemidxs.txt')
midxs = []
for line in inFI:
    midxs.append(int(line.strip()))
inFI.close()

# wilcpvals = [1 for i in range(7440)]
# noiseratiodiffs = [0 for i in range(7440)]
# for i in range(len(ratiosarrretain[0])):
#     noiseratiodiffs[i] = np.mean(ratiosarrretain[0][i])-np.mean(ratiosarrretain[1][i])
#     [stat, pval] = scipy.stats.mannwhitneyu(ratiosarrretain[0][i], ratiosarrretain[1][i])
#     wilcpvals[i] = pval
# writeData([recon2rxns,recon2rxnnames,recon2subs,noiseratiodiffs,wilcpvals],'/mnt/vdb/home/ubuntu2/recon2rxnnoisediffs.txt',headers=['rxn','rxnname','sub','noiseratiodiff','wilcpval'])
rxnfluxdiffs = [0 for i in range(7440)]
rxnpvals = [1 for i in range(7440)]
for cidx in range(1,7441):#cidxs:
    if cidx in reshapedmapretain['old'] and cidx in reshapedmapretain['young']:
        oldfluxes = reshapedmapretain['old'][cidx]['flux']
        youngfluxes = reshapedmapretain['young'][cidx]['flux']
        if max(oldfluxes)!=max(youngfluxes) and min(oldfluxes)!=min(youngfluxes):
            [stat, pval] = scipy.stats.mannwhitneyu(oldfluxes, youngfluxes)
            #print(cidx)
            #print(stat)
            rxnpvals[cidx-1] = pval
            rxnfluxdiffs[cidx-1] = np.mean(youngfluxes)-np.mean(oldfluxes)
            #noiseratiodiffs[cidx-1] = np.mean(youngfluxes)-np.mean(oldfluxes)
            #wilcpvals[cidx-1] = pval

writeData([recon2rxns,recon2rxnnames,recon2subs,rxnfluxdiffs,rxnpvals,ratiosarrretain[0],ratiosarrretain[1]],'/mnt/vdb/home/ubuntu2/recon2rxndiffs.txt',delimiter='\t',headers=['rxn','rxnname','sub','fluxdiff','wilcoxon p-val','young flux to expression ratio','old flux to expression ratio'])
            
[subDiffArr,subDiffArrNum,subDiffMapIndividual,rxnDiffArr,rxnDiffArrNum,rxnDiffMapIndividual] = segmentFluxBySubsystemGroup([0],[1],ratiosarrretain,recon2rxns,recon2rxnnames,recon2subs)

youngratios = np.array(ratiosarrretain[0])
oldratios = np.array(ratiosarrretain[1])
recon2subsTemp = np.array(recon2subs)
subDiffWilcPVals = []
for i in range(len(subDiffArr)):
    if max(youngratios[recon2subsTemp==subDiffArr[i]])==max(oldratios[recon2subsTemp==subDiffArr[i]]) and min(youngratios[recon2subsTemp==subDiffArr[i]])==min(oldratios[recon2subsTemp==subDiffArr[i]]):
        pval = 1
    else:
        [stat, pval] = scipy.stats.mannwhitneyu(youngratios[recon2subsTemp==subDiffArr[i]],oldratios[recon2subsTemp==subDiffArr[i]])
    subDiffWilcPVals.append(pval)

sortIdxs = zip(absArr(subDiffArrNum),range(len(subDiffArrNum)))
sortIdxs.sort(reverse=True)
sortIdxsTemp = []
for i in range(len(sortIdxs)):
    sortIdxsTemp.append(sortIdxs[i][1])
sortIdxs = sortIdxsTemp
subDiffArr = sortByIdx(subDiffArr,sortIdxs)
subDiffArrNum = sortByIdx(subDiffArrNum,sortIdxs)
subDiffWilcPVals = sortByIdx(subDiffWilcPVals,sortIdxs)

writeData([addIdxStrings(subDiffArr),subDiffArrNum,subDiffWilcPVals],'/mnt/vdb/home/ubuntu2/recon2noisediffs.txt',delimiter='\t',headers=['sub','average noise diff','wilcoxon p-val'])

if False:
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

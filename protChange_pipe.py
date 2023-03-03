#!/usr/bin/env python3
import numpy as np
import pandas as pd
from statistics import mean
import os
import sys
# /projectnb/evolution/zwakefield/Annotations/hg38_gencode/hg38_for_df.txt
# /projectnb2/evolution/zwakefield/proteinChange/full_sheet.xlsx
gtf_path = sys.argv[1]
exon_path = sys.argv[2]
outname = sys.argv[3]
ex_type = sys.argv[4]
# import gtf
# Need chr, type (gene, exon, transcript), start, stop, strand, geneID, transcriptID
gtf = pd.read_csv(gtf_path, delimiter="\t", header = None, index_col = False, names = ['chr', 'hav', 'type', 'start', 'stop', 'per', 'strand', 'per2', 'more'])
gtf[['geneID', 'transcriptID']] = gtf.more.str.split(" ", expand=True).loc[:,1:3:2].replace({'"':''}, regex=True).replace({';':''}, regex=True).replace({'transcribed_unprocessed_pseudogene':'gene'}, regex = True)
gtf['geneID'] = gtf.geneID.str.split('.').str[0]
gtf['transcriptID'] = gtf.transcriptID.str.split('.').str[0]
print("gtf loaded...")


ex_count = gtf.more.str.split("exon_id", expand=True).iloc[:,0].str.split("exon_number ", expand=True).iloc[:,1]
l_count_m = []
f_count_m = []
for i in ex_count:
    if i is None:
        l_count_m.append(-1)
#     elif isinstance(i, float):
#         l_count_m.append(int(i))
    else:
        l_count_m.append(int(i.strip("; ")))
for i in ex_count:
    if i is None:
        f_count_m.append(99999)
#     elif isinstance(i, float):
#         l_count_m.append(int(i))
    else:
        f_count_m.append(int(i.strip("; ")))       
        
transGTF = gtf[(gtf.type == 'transcript')]
lastEx = []
firstEx = []
for i in range(len(transGTF.index)):
    if i == (len(transGTF.index)-1):
        end = gtf.shape[0]
    else:
        end = (transGTF.index[i+1])
    lastEx.append(transGTF.index[i]+(np.argmax(l_count_m[transGTF.index[i]:end])))
    firstEx.append(transGTF.index[i]+(np.argmin(f_count_m[transGTF.index[i]:end])))




exon = pd.read_excel(exon_path)
exon[['chr', 'loc']] = exon.exon.str.split(":", expand=True)
exon[['start', 'stop']] = exon['loc'].str.split("-", expand=True)
exon[['geneR', 'rand']] = exon['gene'].str.split("[.]", expand=True)
redExon = exon[['geneR', 'chr', 'start', 'stop']]
redExon.dropna(
    axis=0,
    how='any',
    thresh=None,
    subset=None,
    inplace=True
)
print("exon loaded...")

def flatten(l):
    return [item for sublist in l for item in sublist]

def getTranscripts(gtf, redExon, ex_type):
    ind = []
    transID = []
    outGeneID = []
    counter = 0
    total = 0
    hybCapture = 0
    print(f"searching for {ex_type}...")
    if ex_type == 'AFE':
        pList = firstEx
    else:
        pList = lastEx
    for n, row in redExon.iterrows():
        if n % round(redExon.shape[0]/10, 0) == 0:
            print(str(n) + " exons matched, " + str(n/redExon.shape[0]) + " completed")
        currExon = []
        minExon = []
        for m, gtfRow in gtf[(gtf.geneID == row['geneR']) & (gtf.type == 'exon')].iterrows():
            if len(gtfRow) > 0:
                currExon.append(len(np.intersect1d(range(int(row['start']), (int(row['stop'])+1)), range(int(gtfRow['start']), (int(gtfRow['stop'])+1))))/len(np.union1d(range(int(row['start']), (int(row['stop'])+1)), range(int(gtfRow['start']), (int(gtfRow['stop'])+1)))))
                minExon.append(len(np.intersect1d(range(int(row['start']), (int(row['stop'])+1)), range(int(gtfRow['start']), (int(gtfRow['stop'])+1))))/np.minimum(np.absolute(int(row['start'])-(int(row['stop'])+1)), np.absolute(int(gtfRow['start'])-(int(gtfRow['stop'])+1))))
        if len(gtf[(gtf.geneID == row['geneR']) & (gtf.type == 'exon')]) > 0:
#             topMatch = gtf[(gtf.geneID == row['geneR']) & (gtf.type == 'exon')].iloc[[index for index in range(len(minExon)) if minExon[index] == max(minExon)]]
            topMatch = gtf[(gtf.geneID == row['geneR']) & (gtf.type == 'exon')].iloc[[index for index in range(len(currExon)) if currExon[index] == max(currExon)]]
            ## if matches one best
            if (sorted(currExon, key = lambda x:float(x))[::-1][0] > sorted(currExon, key = lambda x:float(x))[::-1][1]) and ((np.array(gtf[(gtf.geneID == row['geneR']) & (gtf.type == 'exon')].iloc[[int(np.argmax(currExon))]].index)) in pList):
                ind.append(np.argmax(currExon))
                transID.append(gtf[(gtf.geneID == row['geneR']) & (gtf.type == 'exon')]['transcriptID'].iloc[[int(np.argmax(currExon))][0]])
                counter += 1
                total += 1
            ## if matches multiple best, with all as first exons
            elif (sorted(currExon, key = lambda x:float(x))[::-1][0] == sorted(currExon, key = lambda x:float(x))[::-1][1] and sorted(currExon, key = lambda x:float(x))[::-1][0] > .8) and ((np.array(gtf[(gtf.geneID == row['geneR']) & (gtf.type == 'exon')].iloc[[int(np.argmax(currExon))]].index)) in pList):
                ind.append(np.argmax(currExon))
                transID.append(gtf[(gtf.geneID == row['geneR']) & (gtf.type == 'exon')]['transcriptID'].iloc[[int(np.argmax(currExon))][0]])
                counter += 1
                total += 1
            ## if matches a few, only one is a first exon 
            elif (sorted(currExon, key = lambda x:float(x))[::-1][0] == sorted(currExon, key = lambda x:float(x))[::-1][1] and sorted(currExon, key = lambda x:float(x))[::-1][0] > .8) and any(i in pList for i in (np.array(topMatch.index))):
    #                 gtf.iloc[[index for index in topMatch.index if index in np.array(gtf.index[(gtf.type == 'transcript')] + 1)]]['transcriptID'].values[0]
                ind.append(np.argmax(currExon))
                transID.append(gtf.iloc[[index for index in topMatch.index if index in pList]]['transcriptID'].values[0])
                counter += 1
                total += 1
                hybCapture += 1
            else:
                ind.append("None")
                transID.append("None")
                total += 1
                outGeneID.append(row['geneR'] + ";" + row['chr'] + ":" + row['start'] + "-" + row['stop'])
        else:
            ind.append("None")
            transID.append("None")
            outGeneID.append(row['geneR'] + ";" + row['chr'] + ":" + row['start'] + "-" + row['stop'])
            total += 1
    redExon[['transcript']] = transID
    outExons = redExon[(redExon.transcript != "None")].groupby('geneR').filter(lambda x: len(x) > 1)

    tID = []
    for x in outExons.transcript:
        if x != "None":
            tID.append(x)
        else:
            tID.append("None")

    eiID = []
    for n, row in outExons.iterrows():
        if x != "None":
            eiID.append(row['geneR'] + ";" + row['chr'] + ";" + row['start'] + ";" + row['stop'])
        else:
            eiID.append("None")

    print(str(counter) + " out of " + str(total) + " are one-to-one... = " + str(round(counter/total, 3)))
    print("Genes with poor alignment: ")
    for o in outGeneID:
        print(o)
    print("hybrid capture attempt count: " + str(hybCapture))
    return ind, transID, outExons, tID, eiID, outGeneID




def bedify(outExons, tID, eiID, gtf, saveBED = True, outname = "toBed"):
    up = 0
    eiIDa = [[eiID[up]] * gtf[(gtf.transcriptID == tID[0]) & (gtf.type == "exon")][['chr', 'start', 'stop', 'transcriptID', 'geneID', 'strand']].shape[0]]
    bed = gtf[(gtf.transcriptID == tID[0]) & (gtf.type == "exon")][['chr', 'start', 'stop', 'transcriptID', 'geneID', 'strand']]
    for b in tID[1:]:
        up += 1
        if b != "None":
            bed = pd.concat([bed, gtf[(gtf.transcriptID == b) & (gtf.type == "exon")][['chr', 'start', 'stop', 'transcriptID', 'geneID', 'strand']]])
            eiIDa.append([eiID[up]] * gtf[(gtf.transcriptID == b) & (gtf.type == "exon")][['chr', 'start', 'stop', 'transcriptID', 'geneID', 'strand']].shape[0])
    bed['eiID'] = flatten(eiIDa)
    bed['score'] = 0
    bed['squish'] = bed['transcriptID'] + "#" + bed['eiID']

    bed.rename(
        columns={"chr": "chrom", "start": "chromStart", "stop": "chromEnd", "squish": "name"},
        inplace=True,
    )

    toBed = bed[['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand']]

    toBed[['chromStart']] = toBed[['chromStart']] + 1
    toBed['chromStart'] = toBed['chromStart'].astype(int)
    toBed['chromEnd'] = toBed['chromEnd'].astype(int)

    toBed
    if saveBED == True:
        toBed.to_csv(f"./{outname}.bed", sep='\t', header=None, quoting=None, index=False)

    return toBed



ind, transID, outExons, tID, eiID, outGeneID = getTranscripts(gtf = gtf, redExon = redExon, ex_type = ex_type)
print("Transcripts extracted...")
toBed = bedify(outExons=outExons, tID=tID, eiID=eiID, gtf=gtf, saveBED = True, outname = outname)
print("Bed file written...")
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
# import gtf
# Need chr, type (gene, exon, transcript), start, stop, strand, geneID, transcriptID
gtf = pd.read_csv(gtf_path, delimiter="\t", header = None, index_col = False, names = ['chr', 'hav', 'type', 'start', 'stop', 'per', 'strand', 'per2', 'more'])
gtf[['geneID', 'transcriptID']] = gtf.more.str.split(" ", expand=True).loc[:,1:3:2].replace({'"':''}, regex=True).replace({';':''}, regex=True).replace({'transcribed_unprocessed_pseudogene':'gene'}, regex = True)
gtf['geneID'] = gtf.geneID.str.split('.').str[0]
gtf['transcriptID'] = gtf.transcriptID.str.split('.').str[0]
print("gtf loaded...")


exon = pd.read_excel("/projectnb2/evolution/zwakefield/proteinChange/full_sheet.xlsx")
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



def getTranscripts(gtf, redExon):
    ind = []
    transID = []
    counter = 0
    total = 0
    for n, row in redExon.iterrows():
        currExon = []
        for m, gtfRow in gtf[(gtf.geneID == row['geneR']) & (gtf.type == 'exon')].iterrows():
            if len(gtfRow) > 0:
                currExon.append(len(np.intersect1d(range(int(row['start']), (int(row['stop'])+1)), range(int(gtfRow['start']), (int(gtfRow['stop'])+1))))/len(np.union1d(range(int(row['start']), (int(row['stop'])+1)), range(int(gtfRow['start']), (int(gtfRow['stop'])+1)))))
        if len(gtf[(gtf.geneID == row['geneR']) & (gtf.type == 'exon')]) > 0:

            if sorted(currExon, key = lambda x:float(x))[::-1][0] > sorted(currExon, key = lambda x:float(x))[::-1][1]:
                ind.append(np.argmax(currExon))
                transID.append(gtf[(gtf.geneID == row['geneR']) & (gtf.type == 'exon')]['transcriptID'].iloc[[int(np.argmax(currExon))][0]])
            else:
                ind.append("None")
                transID.append("None")
        else:
            ind.append("None")
            transID.append("None")
    # print(ind)
    redExon[['transcript']] = transID
    outExons = redExon[(redExon.transcript != "None")].groupby('geneR').filter(lambda x: len(x) > 1)

    tID = []
    for x in outExons.transcript:
        if x != "None":
            tID.append(x)
        else:
            tID.append("None")
    return ind, transID, outExons, tID


def bedify(outExons, tID, gtf, saveBED = True, outname = "toBed"):
    bed = gtf[(gtf.transcriptID == tID[0]) & (gtf.type == "transcript")][['chr', 'start', 'stop', 'transcriptID', 'strand']]
    for b in tID[1:]:
        if b != "None":
            bed = pd.concat([bed, gtf[(gtf.transcriptID == b) & (gtf.type == "transcript")][['chr', 'start', 'stop', 'transcriptID', 'strand']]])

    bed['score'] = 0
    bed.rename(
        columns={"chr": "chrom", "start": "chromStart", "stop": "chromEnd", "transcriptID": "name"},
        inplace=True,
    )

    toBed = bed[['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand']]

    toBed[['chromEnd']] = toBed[['chromEnd']] + 1
    toBed['chromStart'] = toBed['chromStart'].astype(int)
    toBed['chromEnd'] = toBed['chromEnd'].astype(int)
    if saveBED == True:
        toBed.to_csv(f"./{outname}.bed", sep='\t', header=None, quoting=None, index=False)

    return toBed



ind, transID, outExons, tID = getTranscripts(gtf = gtf, redExon = redExon)
toBed = bedify(outExons=outExons, tID=tID, gtf=gtf, saveBED = True, outname = outname)

from RaceID import SCseq

def getExpData(SCseq_object, genes=None):
    expdata = SCseq_object.expdata
    if genes == None:
        genes = SCseq_object.genes
    return expdata[expdata.index.isin(genes)]
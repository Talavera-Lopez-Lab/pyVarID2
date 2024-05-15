from RaceID3_StemID2_class import SCseq

def getExpData(SCseq_object, genes=None):
    if genes == None:
        # the genes attribute must be assigned in the filterdata method however i cant find where
        genes = SCseq_object.genes
    return None
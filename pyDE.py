import pandas as pd
from pydeseq2.dds import DeseqDataSet

class DESeq2NormCounts:
    """
    Differential expression analysis using DESeq2
    """

    def __init__(self,rawCounts:str,metadat:str):
        self.rawCounts = rawCounts
        self.metadat = metadat
        self.mat = None

    def tidyCounts(self):
        df = pd.read_csv(self.rawCounts,skiprows=1,sep="\t")
        df = df.drop(columns=['Chr','Start','End','Strand','Length'])
        self.mat = df.T
        return self.mat
    
    def runDESeq2(self):
        dds = DeseqDataSet(counts = self.mat,
                            metadata = self.metadat,
                            design = "~cond")
        dds.deseq2()
        return dds
    

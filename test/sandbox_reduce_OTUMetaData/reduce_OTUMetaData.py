#!/home/sumin/anaconda3/bin/python


import pandas as pd
import numpy as np
import sys
import os



if __name__ == '__main__':
    """
    Take PI's OTUMetaData.tsv, which has lots of columns, and just extract 3 columns to test
    kb_faprotax and kb_PICRUSt2

    Write that
    """

    
    
    otu_metdat_flpth = sys.argv[1]

    pwd = os.path.dirname(os.path.realpath(__file__))
    otu_metdat_flnm = os.path.basename(otu_metdat_flpth)
    otu_metdat_reduced_flpth = os.path.join(pwd, otu_metdat_flnm[:-4] + '_reduced' + otu_metdat_flnm[-4:])

    otu_metdat_df = pd.read_csv(
        otu_metdat_flpth, sep='\t', header=0, usecols=['#OTU ID', 'FAPROTAX Traits', 'PiCrust2 Traits'])


    print(f"Writing to {otu_metdat_reduced_flpth}")

    otu_metdat_df.to_csv(otu_metdat_reduced_flpth, sep='\t', index=False)
























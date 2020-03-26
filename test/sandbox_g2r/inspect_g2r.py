#!/home/sumin/anaconda3/bin/python


import pandas as pd
import numpy as np
import sys





if __name__ == '__main__':
    """
    Hacked up script, preliminary comparison of kb_faprotax's output to PI's results (not public)
    """

    g2r_flpth = sys.argv[1]
    otu_metdat_flpth = sys.argv[2]

    print(g2r_flpth, otu_metdat_flpth)

    g2r_df = pd.read_csv(g2r_flpth, sep='\t', header=0, comment='#')
    g2r_df.drop_duplicates(inplace=True)
    g2r_df.set_index('record', inplace=True) 

    assert len(set(g2r_df.index)) == len(list(g2r_df.index))

    print(g2r_df)


    func_npArr = np.array(g2r_df.columns)
    funcsAssigned_l = []

    for record, row in g2r_df.iterrows():
        b_l = list(row)
        nonzero_ind_npArr = np.nonzero(b_l)[0]
        funcsAssigned_l.append(':'.join(list(func_npArr[nonzero_ind_npArr])))

    print(funcsAssigned_l[:20])

    g2r_df['funcsAssigned'] = funcsAssigned_l


    #---------------------------------------------------------------


    otu_metdat_df = pd.read_csv(
        otu_metdat_flpth, sep='\t', header=0, usecols=['taxonomy', 'FAPROTAX Traits']).fillna('')
    otu_metdat_df.drop_duplicates(inplace=True)
    otu_metdat_df.set_index('taxonomy', inplace=True)

    assert len(set(otu_metdat_df.index)) == len(list(otu_metdat_df.index))


    #---------------------------------------------------------------



    otu_metdat_df['funcsAssigned'] = g2r_df.loc[otu_metdat_df.index, 'funcsAssigned']
    
    diff = otu_metdat_df[otu_metdat_df['funcsAssigned'] != otu_metdat_df['FAPROTAX Traits']].copy()

    print('before dedup diff:', len(diff))
    diff.drop_duplicates(inplace=True)
    print('after dedup diff:', len(diff))
    print(diff)

    diff_outfile = sys.argv[3]

    diff.to_csv(diff_outfile, index=False)

    #---------------------------------------------------------------------


    def is_subset(row):
        return set(row['FAPROTAX Traits'].split(':')).issubset(row['funcsAssigned'].split(':'))

    assert all(list(diff.apply(lambda row: is_subset(row), axis=1)))


    def highlight_diff(row):
        less = row['FAPROTAX Traits'].split(':')
        more = row['funcsAssigned'].split(':')
        return ':'.join(['<b>' + f + '</b>' if f not in less else f for f in more])


    highlighted_l = list(diff.apply(lambda row: highlight_diff(row), axis=1))

    diff_highlighted_outfile = sys.argv[4]

    with open(diff_highlighted_outfile, 'w') as fp:
        for highlighted in highlighted_l:
            fp.write('<p>' + highlighted + '</p>\n')

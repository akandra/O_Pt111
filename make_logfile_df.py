# %%
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pandas.core.common import flatten
from IPython.display import display
import json

def make_logfile_df(logfile):
    with open(logfile,'r')  as f:
        lines = f.readlines()

    #-----------------------------------------------------
    # create a list of column names by expanding the header
    # based on the type of items in the first logfile data line.
    #-----------------------------------------------------
    col_names = []
    hdr_names = lines[0].split()
    data_line = json.loads(lines[1])

    for i_name, d in enumerate( data_line ):
        if type(d) is list:
            for i_suffix in range(len(d)):
                col_names.append(hdr_names[i_name] + '_' + str(i_suffix))
        else:
            col_names.append(hdr_names[i_name])
    col_names.append('logfile_name')

    #-----------------------------------------------------
    # expand log data by expanding lists in each line
    #-----------------------------------------------------
    log_data=[]

    for i_line in range(1,len(lines)):
        data = [];
        for d in json.loads( lines[i_line] ):
                if type(d) is list:
                    for value in d:
                        data.append(value)
                else:
                    data.append(d)
        data.append(logfile)            
        log_data.append(data)
    return pd.DataFrame(log_data, columns=flatten(col_names))

if __name__ == "__main__":
    p = Path(__file__).parent
    datapath = p
    
    # Get list of folders with data
    logfiles = list(datapath.glob('*/*.log'))
    datadirs = [f.parent for f in logfiles]
    
    df = make_logfile_df(logfiles[0])
    display(df);




# %%

import pandas as pd
import numpy as np
import os
import t3_VDJ as VDJ
import numpy as np
from multiprocessing import Pool, cpu_count
import time
import sys
from functools import partial

start_time = time.time()
try:
    receptor = str(sys.argv[1])
    print("Initiating job for {receptor}".format(receptor=receptor))
    if '@' in receptor:
        g = receptor.split('@')
        receptor = g[0]
        splitn = int(g[1])
        splitid = int(g[2])
    else:
        splitn = 'nosplit'
    dir = str(sys.argv[2]) # input data dir
    cores = cpu_count()
    print("Using {cores} cores".format(cores=cores)) 

    dbpath = str(sys.argv[3]) # vdjdb dir
except Exception:
    print("You probably put a parameter wrong")
    

def build_df(path, files, receptor=None):
    df = pd.DataFrame(columns = ["Read ID", "Read", "Chromosome", "Position", "VID", "V Match", "V Match Percent", "V Match Length", "JID", "J Match", "J Match Percent", "J Match Length", "JUNCTION", "CDR3", "Reason"])
    dfs = []
    for file in files:
        mycols = ["Read ID", "B", "Chromosome", "Position", "E", "F", "G", "H", "I", "Read", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U",  "V", "W", "X", "Y"]
        try:
            df2 = pd.read_csv(path + file, sep='\t', names=mycols)
        except:
            print("File read error {file}".format(file=file))
            pass
        df2 = df2[["Read ID", "Read", "Chromosome", "Position"]]
        df2["Filename"] = file
    
        df2 = df2.sort_values("Position")
        dfs.append(df2)
    combined = pd.concat(dfs)    
    df = df.append(combined, ignore_index=True)
    #df = df.append(df2, ignore_index=True)  
    df = df.reset_index(drop=True)
    if receptor:
        df.receptor_id=receptor
    return df
    
def run_align(i, df):
    th = 80
    read = df["Read"].iloc[i]
    if receptor == "TRA":
        a = VDJ.TRA(read, dbpath, vthreshold=th, jthreshold=th)
    elif receptor == "TRB":
        a = VDJ.TRB(read, dbpath, vthreshold=th, jthreshold=th)
    elif receptor == "TRD":
        a = VDJ.TRD(read, dbpath, vthreshold=th, jthreshold=th)
    elif receptor == "TRG":
        a = VDJ.TRG(read, dbpath, vthreshold=th, jthreshold=th)
    elif receptor == "IGH":
        a = VDJ.IGH(read, dbpath, vthreshold=th, jthreshold=th)
    elif receptor == "IGK":
        a = VDJ.IGK(read, dbpath, vthreshold=th, jthreshold=th)
    elif receptor == "IGL":
        a = VDJ.IGL(read, dbpath, vthreshold=th, jthreshold=th)
    elif receptor == "TRA_UM":
        a = VDJ.TRA(read, dbpath, vthreshold=th, jthreshold=th)    
    elif receptor == "TRB_UM":
        a = VDJ.TRB(read, dbpath, vthreshold=th, jthreshold=th)    
    elif receptor == "TRD_UM":
        a = VDJ.TRD(read, dbpath, vthreshold=th, jthreshold=th)    
    elif receptor == "TRG_UM":
        a = VDJ.TRG(read, dbpath, vthreshold=th, jthreshold=th)    
    elif receptor == "IGH_UM":
        a = VDJ.IGH(read, dbpath, vthreshold=th, jthreshold=th)    
    elif receptor == "IGK_UM":
        a = VDJ.IGK(read, dbpath, vthreshold=th, jthreshold=th)    
    elif receptor == "IGL_UM":
        a = VDJ.IGL(read, dbpath, vthreshold=th, jthreshold=th)     
    try:
        result = a.run()
    except Exception as e:
        print(f"Error --> {e}")
        raise
    #print('{receptor}{i} is finished.'.format(receptor=receptor, i=i))
    return result
    
    
if __name__  == '__main__':
    path = dir + receptor + '/'
    files = [f for f in os.listdir(path) if os.path.isfile(os.path.join(path, f))]
    print("Found {nfiles} files in {path}.".format(nfiles=len(files), path=path))
    fulldf = build_df(path, files)
    del files
    df = fulldf.drop_duplicates(["Read"], keep='first')
    del fulldf
    if splitn != 'nosplit':
        df = np.array_split(df,splitn)[splitid]
    df = df.reset_index(drop=True)
    print("{0} Starting. Aligning {1} Reads".format(receptor, len(df)))
    results = pd.DataFrame()
    
    items = range(len(df["Read"]))
    predicted_hours = (6.38e-6*max(items)**1.86)/(60*60)
    print("====Estimated hours to completion: {predicted_hours:.2f}  ".format(predicted_hours=predicted_hours))

    with Pool(cores+2) as p:
        
        source = p.map(partial(run_align, df=df), items, 1) # returns a map of series objects
    '''
    Pretty sure that the map command will split the #/core_ct and send each portion to each cpu. This means nothing returns untill the whole chunk is done.
    Maybe try setting chunksize to 1 explicitly? But may make things slower, but probably not because the data size is quite low.

    '''
    print("----------Done Aligning---------")
    
    results = pd.concat(list(source), ignore_index=True,axis=1).transpose()

    cols_to_use = df.columns.difference(results.columns)

    df2 = pd.concat([df[cols_to_use], results], axis=1)
    fulldf = df2[["Filename", "Read ID", "Read", "Chromosome", "Position", "VID", "V Match", "V Match Percent", "V Match Length", "JID", "J Match", "J Match Percent", "J Match Length", "JUNCTION", "CDR3", "Reason"]]
    '''
    try:
        results = pd.concat(source, axis=1, ignore_index=True).transpose()
    except:
        pass
    cols_to_use = df.columns.difference(results.columns)
    df2 = pd.concat([df[cols_to_use], results], axis=1)
    fulldf = fulldf[["Read ID", "Read"]]
    df2 = df2.drop("Read ID", axis = 1)
    fulldf = pd.merge(fulldf, df2, how='left', on='Read')
    fulldf = fulldf[["Filename", "Read ID", "Read", "Chromosome", "Position", "VID", "V Match", "V Match Percent", "V Match Length", "JID", "J Match", "J Match Percent", "J Match Length", "JUNCTION", "CDR3", "Reason"]]
    '''
    outputfile = dir + "{receptor}.csv".format(receptor=receptor)
    if splitn != "nosplit":
       outputfile = dir + "{receptor}${splitid}.csv".format(receptor=receptor, splitid=splitid)
    print("Writing output csv")
    fulldf.to_csv(outputfile)
    del df2
    del fulldf
    del df
    del results
    sec = (time.time() - start_time)
    hr = sec/(60*60)
    print("{receptor} Complete. Time since start: --- {sec:.1f} seconds --- {hr:.2f} hours ----".format(receptor=receptor, sec=sec, hr=hr))

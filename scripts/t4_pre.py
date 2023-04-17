# python 3 script konrad cios
# combines and converts all the csv files in the current dir into xlsx files
# 
import pandas as pd 
import os
import glob
import concurrent.futures
import sys
###############################

foldername = 'XLSX'
combinedf = 'final_csv'
datafolder ='/work/pi_gblanck/Arpan/NBL_WXS/NBL_BAMS_Results'
#####################################

if len(sys.argv) > 1:
    print(sys.argv[1])
    datafolder = sys.argv[1] + '/'

else:
    datafolder = datafolder + '/'
    foldername = datafolder + foldername 
    combinedf = datafolder + combinedf 


print('----------reading files---------------')
files = glob.glob(datafolder+'*.csv')
print(f'files = {files}')

uniqueset = []

for f in files:
    rectxt = f.split('$')[0].split('/')[-1]
    if '$' in f:
        if rectxt not in uniqueset:
            uniqueset.append(rectxt)
    else: # For files that were not large enoguh to be parititoned.
        if '.csv' in rectxt:
            rectxt = rectxt.split('.')[0]
            if rectxt not in uniqueset:
                uniqueset.append(rectxt)


print(f"Found receptors: {uniqueset}")
print("----------------Starting Combining split csv files process.----------------------")

def loader(f):
    d = pd.read_csv(f)
    ret = f.split('$')[0].split('/')[-1]
    if '.csv' in ret:
        ret = ret.split('.')[0]
    #print(ret)
    return (d,ret)
with concurrent.futures.ThreadPoolExecutor(max_workers=19) as executor:
    futures = []
    rec_groups = {}
    for rec in uniqueset:
        rec_groups[rec] = []
    #print(rec_groups)
    for f in files:
        print(f'Loading {f.split("/")[-1]} ')
        futures.append(executor.submit(loader,f=f))
    for future in concurrent.futures.as_completed(futures):
        res = future.result()
        print(f"Loaded {res[1]} part.")
        rec_groups[res[1]].append( res[0])
    executor.shutdown(wait=True)
    
print("-------------------------Done reading csvs.--------------------------")

# print(f'Making XLSX results dir. -->>>  {foldername}')
# try:
#     os.mkdir(foldername)
# except OSError as e:
#     print('folder already exists')
print(f'Making CSV results dir. -->>>  {combinedf}')
try:
    os.mkdir(combinedf)
except OSError as e:
    print('folder already exists')
print("--------------Combining csvs-------------------")
def combine(rect):
    num = len(rec_groups[rect])
    print(f"Writing {rect} ")

    dfs = rec_groups[rect]
    fulldf = pd.concat(dfs)
    fulldf.to_csv(combinedf+'/'+f"{rect}.csv")
    return f"finished writing {rect}.csv"
    
with concurrent.futures.ThreadPoolExecutor(max_workers=19) as executor:
    futures = []
    for rect in rec_groups:
        futures.append(executor.submit(combine,rect=rect))
    for future in concurrent.futures.as_completed(futures):
        print(future.result())
    executor.shutdown(wait=True)
del rec_groups
print("------------Finished combining csvs-----------")
print(f"All the final csv are now in {combinedf}")

# print('---------------reading combined CSV files ----------')
# files = glob.glob(combinedf+'/'+'*.csv')

# def processFile(file):
#     filename = file.split('.')[0].split('/')[-1]
#     if not file.endswith('.csv'): # pass non csvs
#         return 0
#     if '$' in file:
#         return 1
#     print(f'Reading {filename} for conversion.')
#     a = pd.read_csv(file)
    
#     newname =  foldername + '/' + filename + '.xlsx'
    
#     a = a.iloc[:,2:] # drop first two useless index columns
#     a.to_excel(newname,index=False)
#     print(f'Wrote {filename} to {newname}.')
#     return 2
# print("================__Converting to xlsx__================")
# with concurrent.futures.ThreadPoolExecutor(max_workers=20) as executor:
#     futures = []
#     for file in files:
#         futures.append(executor.submit(processFile, file=file))
#     for future in concurrent.futures.as_completed(futures):
#         a = future.result()


print('***************************')
print('*********done**************')
print('***************************')
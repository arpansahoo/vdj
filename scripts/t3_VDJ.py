from Bio.Seq import Seq
from Bio import pairwise2
import Bio.SeqIO as SeqIO
import pandas as pd
import re
import numpy as np
import os

def nt_matches(stringa, stringb):#returns tuple: (% match, nt/total string, True/False VJ > 20)
    u = zip(stringa, stringb) # This function might have terrible performance making it n^2, not sure though.
    y=[]
    for i,j in u: # iterating over arrays? ugh....
        if i==j:
            pass
        else: 
            y.append(j)
    try:
        percentj = float(len(stringa)-len(y))/float(len(stringb))
        percentj = percentj*100
        ntmatch = str((len(stringa)-len(y))) + "/" + str((len(stringb))) + " nt"
        vj = float(len(stringa)-len(y))
        return (percentj, ntmatch, vj)
    except:
        return (np.nan, np.nan, np.nan)
        
class VDJ(object):
    
    def score_only_align(self, subtype):
        if subtype == "V":
            ref_seqs = self.v_seqs
        elif subtype == "J":
            ref_seqs = self.j_seqs
        id_list = []
        alignments_list = []
        for ref in ref_seqs:
            alignments = pairwise2.align.localms(ref.seq, self.test_seq, 5, -10, -10, -10, score_only=True)
            alignments_list.append(alignments)
            id_list.append(ref.id)
        tmp = pd.DataFrame()
        tmp['%sID'%subtype] = id_list
        tmp['%sSCORE'%subtype] = alignments_list
        tmp = tmp.sort_values('%sSCORE'%subtype, ascending=False)
        tmp = tmp.reset_index(drop=True)
        tmp = tmp.iloc[0]
        return tmp
          
    def full_align(self, subtype):
        if subtype == "V":
            ref_seqs = self.v_seqs
        elif subtype == "J":
            ref_seqs = self.j_seqs
        id_list = []
        alignments_list = []
        for ref in ref_seqs:
            alignments = pairwise2.align.localms(ref.seq, self.test_seq, 5, -10, -10, -10, score_only=False)
            for alignment in alignments:
                alignments_list.append(alignment)
                id_list.append(ref.id)
        tmp = pd.DataFrame()
        tmp['%sID'%subtype] = id_list
        tmp['%sALIGNMENT'%subtype] = alignments_list
        tmp[['%sREFSEQ'%subtype, '%sTESTSEQ'%subtype, '%sSCORE'%subtype, '%sALIGNSTART'%subtype, '%sALIGNSTOP'%subtype]] = tmp['%sALIGNMENT'%subtype].apply(pd.Series)
        test_hit = lambda x : x['%sTESTSEQ'%subtype][x['%sALIGNSTART'%subtype]:x['%sALIGNSTOP'%subtype]]
        ref_hit = lambda x : x['%sREFSEQ'%subtype][x['%sALIGNSTART'%subtype]:x['%sALIGNSTOP'%subtype]]
        remainder_hit = lambda x : (x['%sTESTSEQ'%subtype][:x['%sALIGNSTOP'%subtype]])+(x['%sREFSEQ'%subtype][x['%sALIGNSTOP'%subtype]:])
        tmp['%sTESTALIGNSEQ'%subtype] = tmp.apply(test_hit, axis = 1)
        tmp['%sREFALIGNSEQ'%subtype] = tmp.apply(ref_hit, axis = 1)
        if subtype == "J":
            tmp['%sREMAINDERSEQ'%subtype] = tmp.apply(remainder_hit, axis = 1)
        tmp = tmp.sort_values('%sSCORE'%subtype, ascending=False)
        tmp = tmp.reset_index(drop=True)
        tmp = tmp.iloc[0]
        tmp["%s Match Percent"%subtype], tmp["%s Match"%subtype], tmp["%s Match Length"%subtype] = nt_matches(tmp['%sREFALIGNSEQ'%subtype], tmp['%sTESTALIGNSEQ'%subtype]) 
        tmp = tmp.drop(labels=["{subtype}ALIGNMENT".format(subtype=subtype)]) # try dropping this, seems to be causing crashes
        return tmp # might be returning too much data we never use, can drop it before that?
      
    def score_test(self):
        jdf = self.score_only_align("J")
        if jdf["JSCORE"] <= self.jthreshold:
            jdf['JID'] = "No Match"
            self.results = jdf
            self.results["JUNCTION"], self.results["CDR3"], self.results["REASON"] = "No Match", "No Match", "No Match"
            return False
        elif jdf["JSCORE"] >= self.jthreshold:
            vdf = self.score_only_align("V")
            if vdf["VSCORE"] <= self.vthreshold:
                vdf['VID'] = "No Match"
                self.results = pd.concat([vdf,jdf])
                self.results["JUNCTION"], self.results["CDR3"], self.results["REASON"] = "No Match", "No Match", "No Match"
                return False
            elif vdf["VSCORE"] >= self.vthreshold:
                return True
    
    def threading(self, df, CDR3, threading_db_path="/work/pi_gblanck/Arpan/UVM/scripts/threading_db"):
        def lcs(S,T):
            try:
                m = len(S)
                T = T[-20:]
                n = len(T)
                counter = [[0]*(n+1) for x in range(m+1)]
                longest = 0
                lcs_set = set()
                for i in range(m):
                    for j in range(n):
                        if S[i] == T[j]:
                            c = counter[i][j] + 1
                            counter[i+1][j+1] = c
                            if c > longest:
                                lcs_set = set()
                                longest = c
                                lcs_set.add(S[i-c+1:i+1])
                            elif c == longest:
                                lcs_set.add(S[i-c+1:i+1])

                match_list = list(lcs_set)
                match_list = [i for i in match_list if i.startswith("C")]

                return list(match_list)[0]
            except:
                return None

        def compute_vcdr3(df, CDR3, threading_db_path):
            for side in ["V"]:
                accession_id = df["%sID"%side].split("|")
                df["%s_Accession"%side] = accession_id[0]
                df["%sID"%side] = accession_id[1]
                df2 = pd.read_csv(os.path.join(threading_db_path, "%s_data.csv"%side))
                df2 = df2.drop_duplicates("%s.Name"%side)
                df2["%sID"%side] = df2["%s.Name"%side]
                df2 = df2[["%sID"%side, "%s.AA.String"%side]]
                df = pd.DataFrame([df.tolist()], columns=df.index)
                df = df.merge(df2, on="%sID"%side, how='left')

            df["V_CDR3"] = df.apply(lambda x: lcs(CDR3, x["V.AA.String"]), axis=1)
            df['V_CDR3'] = df['V_CDR3'].fillna("C")
            return df['V_CDR3'].item()

        df_copy = df.copy(deep = True)
        return compute_vcdr3(df_copy, CDR3, threading_db_path)

    def junction_find(self):
        junctions = re.search(r'((tgt|tgc)(...){4,20}(ttt|ttc|tgg)(ggt|ggc|gga|ggg)(...)(ggt|ggc|gga|ggg))', self.results['JREMAINDERSEQ'])
        if junctions is None:
            self.results["JUNCTION"], self.results["CDR3"], self.results["REASON"] = "Unproductive", "Unproductive", "Frame"
        else:
            junction = junctions.group(0)[:len(junctions.group(0))-9]
            translation = str(Seq(junction).ungap('-').translate())

            if '*' not in translation and (translation[-1] == 'W' or translation[-1] == 'F'):
                self.results["JUNCTION"], self.results["CDR3"], self.results["REASON"] = junction, translation, "NA"
            else:
                self.results["JUNCTION"], self.results["CDR3"], self.results["REASON"] = "Unproductive", "Unproductive", "Stop" 
 
    def run(self):
        score_test_pass = self.score_test()
        if score_test_pass == False:
            return self.results
        elif score_test_pass == True:
            vdf = self.full_align("V")
            jdf = self.full_align("J")
            #print("full aligned")
            self.results = pd.concat([vdf, jdf])
            self.junction_find()
            return self.results
 
 
class TRA(VDJ):
    def __init__(self, sequence, dbpath, vthreshold=65, jthreshold=65):
        self.test_seq = Seq(sequence).lower()
        self.v_seqs = list(SeqIO.parse("%sTRAV.fasta" %dbpath, "fasta"))
        self.j_seqs = list(SeqIO.parse("%sTRAJ.fasta" %dbpath, "fasta"))
        self.vthreshold = vthreshold
        self.jthreshold = jthreshold
        #print("initialized")
        
class TRB(VDJ):
    def __init__(self, sequence, dbpath, vthreshold=65, jthreshold=65):
        self.test_seq = Seq(sequence).lower()
        self.v_seqs = list(SeqIO.parse("%sTRBV.fasta" %dbpath, "fasta"))
        self.j_seqs = list(SeqIO.parse("%sTRBJ.fasta" %dbpath, "fasta"))
        self.vthreshold = vthreshold
        self.jthreshold = jthreshold

class TRG(VDJ):
   def __init__(self, sequence, dbpath, vthreshold=65, jthreshold=65):
        self.test_seq = Seq(sequence).lower()
        self.test_seq = self.test_seq.reverse_complement()
        self.v_seqs = list(SeqIO.parse("%sTRGV.fasta" %dbpath, "fasta"))
        self.j_seqs = list(SeqIO.parse("%sTRGJ.fasta" %dbpath, "fasta"))
        self.vthreshold = vthreshold
        self.jthreshold = jthreshold

class TRD(VDJ):
    def __init__(self, sequence, dbpath, vthreshold=65, jthreshold=65):
        self.test_seq = Seq(sequence).lower()
        self.v_seqs = list(SeqIO.parse("%sTRDV.fasta" %dbpath, "fasta"))
        self.j_seqs = list(SeqIO.parse("%sTRDJ.fasta" %dbpath, "fasta"))
        self.vthreshold = vthreshold
        self.jthreshold = jthreshold
        
class IGH(VDJ):
    def __init__(self, sequence, dbpath, vthreshold=65, jthreshold=65):
        self.test_seq = Seq(sequence).lower()
        self.test_seq = self.test_seq.reverse_complement()
        self.v_seqs = list(SeqIO.parse("%sIGHV.fasta" %dbpath, "fasta"))
        self.j_seqs = list(SeqIO.parse("%sIGHJ.fasta" %dbpath, "fasta"))
        self.vthreshold = vthreshold
        self.jthreshold = jthreshold
    
class IGL(VDJ):
    def __init__(self, sequence, dbpath, vthreshold=65, jthreshold=65):
        self.test_seq = Seq(sequence).lower()
        self.v_seqs = list(SeqIO.parse("%sIGLV.fasta" %dbpath, "fasta"))
        self.j_seqs = list(SeqIO.parse("%sIGLJ.fasta" %dbpath, "fasta"))
        self.vthreshold = vthreshold
        self.jthreshold = jthreshold

class IGK(VDJ):
    def __init__(self, sequence, dbpath, vthreshold=65, jthreshold=65):
        self.test_seq = Seq(sequence).lower()
        self.v_seqs = list(SeqIO.parse("%sIGKV.fasta" %dbpath, "fasta"))
        self.j_seqs = list(SeqIO.parse("%sIGKJ.fasta" %dbpath, "fasta"))
        self.vthreshold = vthreshold
        self.jthreshold = jthreshold


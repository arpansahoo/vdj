import pandas as pd
import numpy as np
from openpyxl import load_workbook
import os
import sys
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.SeqUtils import ProtParamData
from Bio import SeqIO
from localcider.sequenceParameters import SequenceParameters
from lifelines.statistics import logrank_test
from lifelines import KaplanMeierFitter
from lifelines import CoxPHFitter
from scipy import stats
import re


class VDJrecord(object):
    def __init__(
        self, cancer="TCGA-COAD", samples_path=None, vdjdb_path=None, clinical_path=None
    ):
        self.cancer = cancer
        self.samples_path = samples_path
        self.vdjdb_path = vdjdb_path
        self.clinical_path = clinical_path
        self.projectid = cancer

    def load_raw_excel_path(self, path):
        print("Loading excel files from path: %s" % path)
        self.raw = pd.DataFrame()
        for file in os.listdir(path):
            if (
                file.endswith(".xlsx")
                and not file.startswith("~")
                and not file.startswith("V")
            ):
                filepath = os.path.join(path, file)
                print("Loading: %s" % filepath)
                df = pd.read_excel(filepath)
                df["Receptor"] = file.replace(".xlsx", "")
                self.raw = pd.concat([self.raw, df])
        self.filtered = self.raw
        print("Loading excel files from path complete")
        return self.raw

    def load_raw_csv_path(self, path):
        print("Loading csv files from path %s" % path)
        self.raw = pd.DataFrame()
        for file in os.listdir(path):
            if (
                file.endswith(".csv")
                and not file.startswith("~")
                and not file.startswith("V")
            ):
                filepath = os.path.join(path, file)
                print("Loading: %s" % filepath)
                df = pd.read_csv(filepath)
                df["Receptor"] = file.replace(".csv", "")
                self.raw = pd.concat([self.raw, df])
        self.filtered = self.raw
        print("Loading csv files from path complete")
        return self.raw

    def load_raw_hdf_path(self, path):
        print("Loading hdf files from path: %s" % path)
        self.raw = pd.DataFrame()
        for file in os.listdir(path):
            if (
                file.endswith(".h5")
                and not file.startswith("~")
                and not file.startswith("v")
            ):
                filepath = os.path.join(path, file)
                print("Loading: %s" % filepath)
                df = pd.read_hdf(filepath, "raw")
                df["Receptor"] = file.replace(".h5", "")
                self.raw = pd.concat([self.raw, df])
        self.filtered = self.raw
        print("Loading hdf files from path complete")
        return self.raw

    def save_hdf(self, filepath):
        print("Saving data to hdf file: %s" % filepath)
        if hasattr(self, "raw"):
            print("Saving raw data to hdf")
            self.raw.to_hdf(filepath, "raw")
        if hasattr(self, "filtered_productive"):
            print("Saving filtered productive data to hdf")
            self.filtered_productive.to_hdf(filepath, "filtered_productive")
        if hasattr(self, "filtered_unproductive"):
            print("Saving filtered unproductive data to hdf")
            self.filtered_unproductive.to_hdf(filepath, "filtered_unproductive")
        if hasattr(self, "physicochem"):
            print("Saving physicochemical data to hdf")
            self.physicochem.to_hdf(filepath, "physicochem")
        if hasattr(self, "vdjdbmatch"):
            print("Saving VDJDB matches to hdf")
            self.vdjdbmatch.to_hdf(filepath, "vdjdbmatch")

    def save_excel(self, filepath):
        print("Saving data to excel file: %s" % filepath)
        writer = pd.ExcelWriter(filepath)
        if hasattr(self, "raw"):
            print("Saving raw data to excel")
            self.raw.to_excel(writer, "raw", index=False)
        if hasattr(self, "filtered_productive"):
            print("Saving filtered productive data to excel")
            self.filtered_productive.to_excel(
                writer, "filtered_productive", index=False
            )
        if hasattr(self, "filtered_unproductive"):
            print("Saving filtered unproductive data to excel")
            self.filtered_unproductive.to_excel(
                writer, "filtered_unproductive", index=False
            )
        if hasattr(self, "physicochem"):
            print("Saving physicochemical data to excel")
            self.physicochem.to_excel(writer, "physicochem", index=False)
        if hasattr(self, "vdjdbmatch"):
            print("Saving VDJDB matches to excel")
            self.vdjdbmatch.to_excel(writer, "vdjdbmatch", index=False)
        writer.save()

    def save_csv(self, filepath):
        print("saving data to csv %s" % filepath)
        self.filtered_productive.to_csv(filepath, index=False)

    def save_csv_pk(self, filepath):
        print("saving data to csv %s" % filepath)
        self.physicochem.to_csv(filepath, index=False)

    def load_hdf(self, filepath):
        print("Loading data from hdf file: %s" % filepath)
        try:
            print("Loading raw data from hdf file")
            self.raw = pd.read_hdf(filepath, "raw")
            self.filtered = self.raw
        except:
            raise ValueError("Hdf raw file doesn't exist")
        try:
            print("Loading filtered productive data from hdf file")
            self.filtered_productive = pd.read_hdf(filepath, "filtered_productive")
        except:
            print("Filtered productive data doesn't exist in file")
        try:
            print("Loading filtered unproductive data from hdf file")
            self.filtered_unproductive = pd.read_hdf(filepath, "filtered_unproductive")
        except:
            print("Filtered unproductive data doesn't exist in file")
        try:
            print("Loading physicochemical data from hdf file")
            self.physicochem = pd.read_hdf(filepath, "physicochem")
        except:
            print("Physicochemical data doesn't exist in file")
        try:
            print("Loading VDJDB matches from hdf file")
            self.vdjdbmatch = pd.read_hdf(filepath, "vdjdbmatch")
        except:
            print("VDJDB matches don't exist in file")

    def nt_match_length_filter(self, V_length=9, J_length=9):
        print("Filtering match lengths: V = {0}, J = {1}".format(V_length, J_length))
        self.filtered = self.filtered.loc[self.filtered["V Match Length"] >= V_length]
        self.filtered = self.filtered.loc[self.filtered["J Match Length"] >= J_length]
        return self.filtered

    def nt_match_percent_filter(self, V_percent=90, J_percent=90):
        print(
            "Filtering match percents: V = {0}%, J = {1}%".format(V_percent, J_percent)
        )
        self.filtered = self.filtered.loc[self.filtered["V Match Percent"] >= V_percent]
        self.filtered = self.filtered.loc[self.filtered["J Match Percent"] >= J_percent]
        return self.filtered

    def filename_format(self):  # not used
        if self.samples_path == None:
            raise ValueError("Must define path to samples db")
        else:
            print("Loading samples db: %s" % self.samples_path)
            self.samples = pd.read_csv(self.samples_path)
            self.samples = self.samples[self.samples["Project ID"] == self.projectid]
            filenamelist = []
            samplelist = []
            sampletypelist = []
            self.filtered["Filename"] = self.filtered["Filename"].str.replace(
                "sliced_", ""
            )
            self.filtered["Filename"] = self.filtered["Filename"].str.replace(
                ".tsv", ""
            )
            print("Matching filenames to samples db")
            for i in self.filtered["Filename"]:
                filename = self.samples.loc[
                    self.samples["File Name"] == i,
                    ["Case ID", "Sample ID", "Sample Type"],
                ]
                filenamelist.append(filename["Case ID"].iloc[0])
                samplelist.append(filename["Sample ID"].iloc[0][-3:])
                sampletypelist.append(filename["Sample Type"].iloc[0])
            self.filtered["Filename"] = filenamelist
            self.filtered.insert(1, "Sample", samplelist)
            self.filtered.insert(2, "Sample Type", sampletypelist)
            return self.filtered

    def raw_filename_format(self):
        if self.samples_path == None:
            raise ValueError("Must define path to samples db")
        else:
            print("Loading samples db: %s" % self.samples_path)
            self.samples = pd.read_csv(self.samples_path, sep="\t")
            self.samples = self.samples[self.samples["project_id"] == self.projectid]
            filenamelist = []
            samplelist = []
            sampletypelist = []
            self.raw["Filename"] = self.raw["Filename"].str.replace("sliced_", "")
            """
            self.raw["Filename"] = self.raw["Filename"].str.replace("_gdc_realn.bam.tsv","")
            print("Matching filenames to samples db")
            for i in self.raw["Filename"]:
                print("matching %s" %i)
                filename = self.samples.loc[self.samples["sample_id"] == i, ["case_id", "sample_id", "sample_type"]]
                filenamelist.append(filename["case_id"].iloc[0])
                samplelist.append(filename["sample_id"].iloc[0][-3:])
                sampletypelist.append(filename["sample_type"].iloc[0])
            self.raw["Filename"] = filenamelist
            self.raw.insert(1, 'Sample', samplelist)
            self.raw.insert(2, 'Sample Type', sampletypelist)
            """
            self.filtered = self.raw
            return self.raw

    def productive_unproductive_split(self):
        print("Splitting filtered data into productive and unproductive sets")
        self.filtered_unproductive = self.filtered[
            self.filtered["CDR3"] == "Unproductive"
        ]
        self.filtered_productive = self.filtered[
            self.filtered["CDR3"] != "Unproductive"
        ]
        return self.filtered_productive, self.filtered_unproductive

    def full_filter(self):
        self.raw_filename_format()
        # FIXME v and j match length standard
        self.nt_match_length_filter(V_length=14, J_length=14)
        self.nt_match_percent_filter(V_percent=80, J_percent=80)
        # self.filename_format()
        self.productive_unproductive_split()
        return self.filtered_productive, self.filtered_unproductive

    def vdjdb_match(self):
        if self.vdjdb_path == None:
            raise ValueError("Must define path to samples db")
        else:
            print("Matching productive CDR3s to VDJDB")
            tra = self.raw.loc[self.raw["Receptor"] == "TRA"]
            trb = self.raw.loc[self.raw["Receptor"] == "TRB"]
            vdjdb = pd.read_csv(self.vdjdb_path)
            travdjdb = vdjdb[vdjdb["Gene"] == "TRA"]
            trbvdjdb = vdjdb[vdjdb["Gene"] == "TRB"]
            tra = tra.merge(travdjdb, on="CDR3", how="inner")
            trb = trb.merge(trbvdjdb, on="CDR3", how="inner")
            df = pd.concat([tra, trb])
            # df = df.drop_duplicates("Read ID") # FIXME is this necessary?
            df = df.drop(["complex.id", "Gene", "V", "J"], axis=1)
            self.vdjdbmatch = df
            return self.vdjdbmatch

    # FIXME use Bio.Align.PairwiseAligner in conjunction with LCS to align the V and CDR3 sequences
    def threading(self, threading_db_path):
        def lcs(S, T):
            try:
                m = len(S)
                T = T[-20:]
                n = len(T)
                counter = [[0] * (n+1) for x in range(m+1)]
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
                                lcs_set.add(S[(i-c+1):(i+1)])
                            elif c == longest:
                                lcs_set.add(S[(i-c+1):(i+1)])

                match_list = list(lcs_set)
                match_list = [i for i in match_list if i.startswith("C")]

                return list(match_list)[0]
            except:
                return None

        # V_CDR3 is a prefix of the CDR3, so just chop it off the V string
        def chop_v(v_string, v_cdr3):
            cutoff_index = v_string.rfind(v_cdr3)
            if cutoff_index == -1:
                return ""
            return v_string[:cutoff_index]

        # the AAs + phenylalanine preceding the G-G doublet is a suffix of the CDR3, so just chop it off the J string
        def chop_j(j_string):
            phe_gly_doublet = re.search(r"((F|W)(G)(.)(G))", j_string)
            if phe_gly_doublet is None:
                return ""
            return j_string[(phe_gly_doublet.start() + 1):]

        def perform_threading(df, threading_db_path):
            for side in ["V", "J"]:
                df[["%s_Accession" % side, "%sID" % side, "drop"]] = df[
                    "%sID" % side
                ].str.split("|", expand=True)
                df = df.drop("drop", axis=1)
                df2 = pd.read_csv(os.path.join(threading_db_path, "%s_data.csv" % side))
                df2 = df2.drop_duplicates("%s.Name" % side)
                df2["%sID" % side] = df2["%s.Name" % side]
                df2 = df2[["%sID" % side, "%s.AA.String" % side]]
                df = df.merge(df2, on="%sID" % side, how="left")

            df["V.AA.String"] = df["V.AA.String"].fillna("")
            df["J.AA.String"] = df["J.AA.String"].fillna("")

            # "V_CDR3" is the prefix of the CDR3 derived from the V segment
            df["V_CDR3"] = df.apply(lambda x: lcs(x["CDR3"], x["V.AA.String"]), axis=1)
            df["V_CDR3"] = df["V_CDR3"].fillna("C")

            # Build the V-CDR3-J sequence. 
            df["V_CDR3_J"] = df.apply(
                lambda x: chop_v(x["V.AA.String"], x["V_CDR3"]), axis=1
            )
            df["V_CDR3_J"] += df["CDR3"]
            df["V_CDR3_J"] += df.apply(lambda x: chop_j(x["J.AA.String"]), axis=1)

            return df

        self.filtered_productive = perform_threading(
            self.filtered_productive, threading_db_path
        )

    def analyze_physicochem(self):
        print("Analyzing CDR3 physicochemical data")
        # df = self.filtered_productive[["Filename", "CDR3", "Sample", "Sample Type", "Receptor"]]
        df = self.filtered_productive[["Filename", "CDR3", "Receptor"]]
        length = []
        fraction_tiny = []
        fraction_small = []
        fraction_charged = []
        fraction_positive = []
        fraction_negative = []
        fraction_expanding = []
        fraction_aromatic = []
        molecular_weight = []
        isoelectric_point = []
        gravy = []
        aromaticity = []
        instability_index = []
        secondary_structure_helix = []
        secondary_structure_turn = []
        secondary_structure_sheet = []
        mean_hydropathy = []
        uversky_hydropathy = []
        NCPR = []
        kappa = []
        omega = []
        PPII_propensity = []
        delta = []
        fraction_disorder_promoting = []

        def count_aa(sequence, aas):
            count = 0
            for aa in aas:
                count += sequence.count(aa)
            return count

        def sequence_analysis(sequence):
            try:
                cdr3 = ProteinAnalysis(str(sequence))
                cidercdr3 = SequenceParameters(str(sequence))
                molecular_weight = cdr3.molecular_weight()
                isoelectric_point = cdr3.isoelectric_point()
                gravy = cdr3.gravy()
                aromaticity = cdr3.aromaticity()
                instability_index = cdr3.instability_index()
                secondary = cdr3.secondary_structure_fraction()
                secondary_structure_helix = secondary[0]
                secondary_structure_turn = secondary[1]
                secondary_structure_sheet = secondary[2]
                length = len(sequence)
                fraction_tiny = float(count_aa(sequence, "ABCGST")) / float(
                    len(sequence)
                )
                fraction_small = float(count_aa(sequence, "ABCDGNPSTV")) / float(
                    len(sequence)
                )
                fraction_aromatic = float(count_aa(sequence, "FHWY")) / float(
                    len(sequence)
                )
                fraction_charged = cidercdr3.get_FCR()
                fraction_positive = cidercdr3.get_fraction_positive()
                fraction_negative = cidercdr3.get_fraction_negative()
                fraction_expanding = cidercdr3.get_fraction_expanding()
                mean_hydropathy = cidercdr3.get_mean_hydropathy()
                NCPR = cidercdr3.get_NCPR()
                uversky_hydropathy = cidercdr3.get_uversky_hydropathy()
                kappa = cidercdr3.get_mean_hydropathy()
                omega = cidercdr3.get_Omega()
                PPII_propensity = cidercdr3.get_PPII_propensity()
                delta = cidercdr3.get_delta()
                fraction_disorder_promoting = (
                    cidercdr3.get_fraction_disorder_promoting()
                )
                return (
                    molecular_weight,
                    isoelectric_point,
                    gravy,
                    aromaticity,
                    instability_index,
                    secondary_structure_helix,
                    secondary_structure_turn,
                    secondary_structure_sheet,
                    length,
                    fraction_tiny,
                    fraction_small,
                    fraction_aromatic,
                    fraction_charged,
                    fraction_positive,
                    fraction_negative,
                    fraction_expanding,
                    mean_hydropathy,
                    NCPR,
                    uversky_hydropathy,
                    kappa,
                    omega,
                    PPII_propensity,
                    delta,
                    fraction_disorder_promoting,
                )
            except:
                return (
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                )

        for sequence in df["CDR3"]:
            a = sequence_analysis(sequence)
            molecular_weight.append(a[0])
            isoelectric_point.append(a[1])
            gravy.append(a[2])
            aromaticity.append(a[3])
            instability_index.append(a[4])
            secondary_structure_helix.append(a[5])
            secondary_structure_turn.append(a[6])
            secondary_structure_sheet.append(a[7])
            length.append(a[8])
            fraction_tiny.append(a[9])
            fraction_small.append(a[10])
            fraction_aromatic.append(a[11])
            fraction_charged.append(a[12])
            fraction_positive.append(a[13])
            fraction_negative.append(a[14])
            fraction_expanding.append(a[15])
            mean_hydropathy.append(a[16])
            NCPR.append(a[17])
            uversky_hydropathy.append(a[18])
            kappa.append(a[19])
            omega.append(a[20])
            PPII_propensity.append(a[21])
            delta.append(a[22])
            fraction_disorder_promoting.append(a[23])
        df["length"] = length
        df["fraction_tiny"] = fraction_tiny
        df["fraction_small"] = fraction_small
        df["fraction_aromatic"] = fraction_aromatic
        df["fraction_charged"] = fraction_charged
        df["fraction_positive"] = fraction_positive
        df["fraction_negative"] = fraction_negative
        df["fraction_expanding"] = fraction_expanding
        df["fraction_disorder_promoting"] = fraction_disorder_promoting
        df["molecular_weight"] = molecular_weight
        df["isoelectric_point"] = isoelectric_point
        df["ncpr"] = NCPR
        df["mean_hydropathy"] = mean_hydropathy
        df["uversky_hydropathy"] = uversky_hydropathy
        df["ppii_propensity"] = PPII_propensity
        df["gravy"] = gravy
        df["aromaticity"] = aromaticity
        df["kappa"] = kappa
        df["omega"] = omega
        df["delta"] = delta
        df["instability_index"] = instability_index
        df["secondary_structure_helix"] = secondary_structure_helix
        df["secondary_structure_turn"] = secondary_structure_turn
        df["secondary_structure_sheet"] = secondary_structure_sheet
        self.physicochem = df
        return self.physicochem

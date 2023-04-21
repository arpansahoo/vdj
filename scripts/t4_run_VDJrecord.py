import os
import sys

from t4_VDJrecord import VDJrecord


path = sys.argv[1]
vdjrecord = VDJrecord(cancer="TARGET-NBL", samples_path=sys.argv[2])
vdjrecord.load_raw_csv_path(path=os.path.join(path))
vdjrecord.full_filter()
vdjrecord.threading(threading_db_path=sys.argv[3])
vdjrecord.analyze_physicochem()
# vdjrecord.save_hdf(os.path.join(path,'vdj_recoveries.h5'))
vdjrecord.save_csv(os.path.join(path, "VDJ_Recoveries.csv"))
vdjrecord.save_csv_pk(os.path.join(path, "VDJ_Recoveries_pk.csv"))

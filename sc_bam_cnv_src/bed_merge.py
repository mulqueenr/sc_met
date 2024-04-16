import pandas as pd
import os, glob
import scipy.io, scipy.sparse
import argparse, pathlib

parser = argparse.ArgumentParser(description='Generate merged data frame from bed files.')
parser.add_argument('--feat_name', dest='feat_name',
                    help='Outfile feature names. To be appended to the outfile names.')
parser.add_argument('--working_dir', dest='working_dir',
                    help='Working directory for all bed.gz files.')


args = parser.parse_args()

#python /src/bed_merge.py \
#--feat_name met_atlas \
#--working_dir "."

""" FUNCTIONS """
def read_bed_files(bed_file):
	cell_name=bed_file.split(".")[3]
	dat = pd.read_table(bed_file,delimiter="\t",names=['feat', 'cg_count', 'mcg_count', 'rate'],na_values=".")
	dat = pd.DataFrame(dat, columns=['rate'])
	dat = dat.rename(columns={"rate": cell_name})
	return dat

"""				"""

""" TOTAL COVERAGE MATRIX 
All CGs measured that overlap per feature """
#adding filter to remove Y chr
feat_name=args.feat_name
bed_merge = [read_bed_files(bed) for bed in glob.glob(os.path.join(args.working_dir,"*bed.gz"))]
win_feats = pd.read_table(glob.glob(os.path.join(args.working_dir,"*bed.gz"))[0],delimiter="\t",names=['feat', 'cg_count', 'mcg_count', 'rate'],na_values=".")
win_feats = pd.DataFrame(win_feats, columns=['feat'])
bed_merge_out = pd.concat(bed_merge ,axis=1).transpose()
bed_merge_out=bed_merge_out.rename(win_feats['feat'],axis=1)

#output with cells as rows and consistent column names, such that everything can be concatenated together
bed_merge_out.to_csv(feat_name+".read_met.tsv.gz",sep="\t",header=True,index=True,compression="gzip")
"""					"""

#bsub -Is -W 36:00 -q long -n 10 -M 100 -R rusage[mem=100] /bin/bash
#cd met_work/9e/0a49b75fb2fb8b514f92af1e1feb2c/
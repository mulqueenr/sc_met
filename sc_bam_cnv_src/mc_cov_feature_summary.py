from pybedtools import BedTool #
import multiprocessing
import pandas as pd
import numpy as np
import os, glob
import scipy.io, scipy.sparse
import numpy as np
import argparse, pathlib

parser = argparse.ArgumentParser(description='Generate methylation summary statistics per feature.')
parser.add_argument('--features', dest='feat',
                    help='Bed or bed.gz formatted input list of features, must contain a unique feature name as column 4.')
parser.add_argument('--feat_name', dest='feat_name',
                    help='Outfile feature names. To be appended to the outfile names.')
parser.add_argument('--cov_folder', dest='cov_folder',
                    help='Folder output of ScaleMethyl Pipeline containing cov.gz files per chromosome.')
args = parser.parse_args()

""" FUNCTIONS """
def cov_per_feat(dat,win,cellid_i,mc_cov_only=False):
	# calculate methylation coverage for cell at window
	tmp=dat[(dat['cellid']==cellid_i)]
	if mc_cov_only:
		print("Calculating "+dat['cg_chrom'][0]+" mC counts for "+cellid_i)
		result = [tmp[(tmp['win_name']==win_j) & (tmp['cg_val']=="Z")].shape[0] for win_j in win.win_name.unique()]
	else:
		print("Calculating "+dat['cg_chrom'][0]+" total coverage for "+cellid_i)
		result = [tmp[(tmp['win_name']==win_j)].shape[0] for win_j in win.win_name.unique()]
	return result

def cov_per_chrom(cov,feat,mc_cov_only=False):
	dat_cov = BedTool(cov)
	bed = BedTool(feat)  
	overlap = dat_cov.intersect(bed,wa=True,wb=True)
	dat = pd.read_table(overlap.fn,delimiter="\t",names=['cg_chrom', 'cg_start', 'cg_stop', 'cellid', 'cg_val','win_chrom','win_start','win_end','win_name'])
	win = pd.read_table(feat, delimiter="\t",names=['win_chrom','win_start','win_end','win_name'])
	win = win[(win['win_chrom']==dat['cg_chrom'][0])]
	args=[(dat,win,x,mc_cov_only) for x in dat.cellid.unique().tolist()]
	if __name__ == '__main__':
		# generate met coverage pandas df
		with multiprocessing.Pool(processes=96) as pool:
			out=pool.starmap(cov_per_feat,args)
	df_cov = pd.DataFrame(out,columns=win.win_name.unique(),index=dat.cellid.unique())
	df_cov = df_cov.replace({None: np.nan})  # Replacing None with NaN for missing values
	return df_cov 


def posterior_mcrate_estimate(cov_out,mc_out,cellid):
	if sum(cov_out.index == mc_out.index)==len(cov_out.index):
		#from https://github.com/lhqing/ALLCools/blob/master/ALLCools/mcds/utilities.py
		print("Calculating posterior estimate of methylation for "+cellid)
		raw_frac = mc_out[cellid]/cov_out[cellid]
		cell_rate_mean = np.nanmean(raw_frac)
		cell_rate_var = np.nanvar(raw_frac)
		# based on beta distribution mean, var
		# a / (a + b) = cell_rate_mean
		# a * b / ((a + b) ^ 2 * (a + b + 1)) = cell_rate_var
		# calculate alpha beta value for each cell
		cell_a = (1 - cell_rate_mean) * (cell_rate_mean**2) / cell_rate_var - cell_rate_mean
		cell_b = cell_a * (1 / cell_rate_mean - 1)
		post_frac = (mc_out[cellid] + cell_a) / (cov_out[cellid] + cell_a + cell_b)
		prior_mean = cell_a / (cell_a + cell_b)
		post_frac = post_frac / prior_mean
		return(post_frac)
	else:
		print("The rows aren't aligned right.")

"""				"""

""" TOTAL COVERAGE MATRIX 
All CGs measured that overlap per feature """
cov_chr = [cov_per_chrom(cov,args.feat,mc_cov_only=False) for cov in glob.glob(os.path.join(args.cov_folder,"chr*bed.gz"))]
cov_out = pd.concat(cov_chr,axis=1).transpose().sort_index()
sample_name = args.cov_folder.split("/")[-1].split(".")[0]+"_"+args.cov_folder.split("/")[-1].split(".")[1]+"_"
cov_out.columns =[sample_name+x for x in cov_out.columns]
cov_out = cov_out.replace({np.nan : 0})  # NaN with 0 since for coverage that is the same thing
cov_out.to_csv(sample_name+args.feat_name+".total_count.tsv.gz",sep="\t",header=True,index=True,compression="gzip")
"""					"""

"""	METHYLATION CG COVERAGE MATRIX
Z/(z+Z) that overlap per  feature, NaN if < min_count_for_rate """
mc_chr = [cov_per_chrom(cov,args.feat,mc_cov_only=True) for cov in glob.glob(os.path.join(args.cov_folder,"chr*bed.gz"))]
mc_out = pd.concat(mc_chr,axis=1).transpose().sort_index()
sample_name = args.cov_folder.split("/")[-1].split(".")[0]+"_"+args.cov_folder.split("/")[-1].split(".")[1]+"_"
mc_out.columns =[sample_name+x for x in mc_out.columns]
mc_out.to_csv(sample_name+args.feat_name+".mc_count.tsv.gz",sep="\t",header=True,index=True,compression="gzip")
"""				"""


""" NEG BINOM RATE ESTIMATE
#modified from https://github.com/lhqing/ALLCools/blob/master/ALLCools/mcds/utilities.py """
args=[(cov_out,mc_out,x) for x in cov_out.columns]
if __name__ == '__main__':
	# generate met coverage pandas df
	with multiprocessing.Pool(processes=96) as pool:
		out=pool.starmap(posterior_mcrate_estimate,args)
df_posterior = pd.DataFrame(out).transpose()
df_posterior.to_csv(sample_name+args.feat_name+".mc_posteriorest.tsv.gz",sep="\t",header=True,index=True,compression="gzip")


#cd /volumes/seq/projects/metACT/240205_RMMM_scalebiotest2/cg_sort_cov
#for i in $(ls -d ./*/*); do python mc_cov_feature_summary.py --features genome_windows.100kb.bed --cov_folder $i ; done
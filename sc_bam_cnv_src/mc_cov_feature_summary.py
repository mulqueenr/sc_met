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
parser.add_argument('--min_cg', dest='min_cg',
                    help='Minimum count of CG sites for feature reporting (otherwise NA)')

#parser.add_argument('--cpus', dest='task_cpus',
                    #help='Cpus to use for process pool.')

args = parser.parse_args()
#python /src/mc_cov_feature_summary.py \
#--features GRCh38_transcripts.longest.bed \
#--feat_name genebody \
#--cov_folder MCF10A.1A11.CG.chroms.sort #--cpus 5

#feat="GRCh38_transcripts.longest.bed"
#cov_folder="MCF10A.1A11.CG.chroms.sort"
#feat_name="genebody"
#args.feat=feat
#args.cov_folder=cov_folder
#args.feat_name=feat_name

""" FUNCTIONS """
def cov_per_feat(dat,win,cellid_i,mc_cov_only=False):
	# calculate methylation coverage for cell at window
	tmp=dat[(dat['cellid']==cellid_i)]
	if mc_cov_only:
		print("Calculating "+str(dat['cg_chrom'][0])+" mC counts for "+cellid_i)
		result = [tmp[(tmp['win_name']==win_j) & (tmp['cg_val']=="Z")].shape[0] for win_j in win.win_name.unique()]
	else:
		print("Calculating "+str(dat['cg_chrom'][0])+" total coverage for "+cellid_i)
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
		with multiprocessing.Pool(processes=1) as pool:
			out=pool.starmap(cov_per_feat,args)
	df_cov = pd.DataFrame(out,columns=win.win_name.unique(),index=dat.cellid.unique())
	df_cov = df_cov.replace({None: np.nan})  # Replacing None with NaN for missing values
	return df_cov 


def mcrate_simple(cov_out,mc_out,cellid,cutoff):
	if sum(cov_out.index == mc_out.index)==len(cov_out.index):
		print("Calculating methylation rates for "+cellid+" putting NA for coverage less than "+str(cutoff))
		raw_frac = mc_out[cellid]/cov_out[cellid]
		raw_frac = raw_frac.mask(cov_out[cellid] < cutoff)
		return(raw_frac)
	else:
		print("The rows aren't aligned right.")

def posterior_mcrate_estimate(cov_out,mc_out,cellid,cutoff):
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
		post_frac = post_frac.mask(cov_out[cellid] < cutoff)
		return(post_frac)
	else:
		print("The rows aren't aligned right.")

"""				"""

""" TOTAL COVERAGE MATRIX 
All CGs measured that overlap per feature """
#adding filter to remove Y chr
feat_name=args.feat_name
cov_chr = [cov_per_chrom(cov,args.feat,mc_cov_only=False) for cov in glob.glob(os.path.join("./"+args.cov_folder,"chr*bed.gz")) if "chrY" not in cov ]
cov_out = pd.concat(cov_chr,axis=1).transpose().sort_index()
sample_name = args.cov_folder.split("/")[-1].split(".")[0]+"_"+args.cov_folder.split("/")[-1].split(".")[1]
cov_out.columns =[sample_name+x for x in cov_out.columns]
cov_out = cov_out.replace({np.nan : 0})  # NaN with 0 since for coverage that is the same thing
#output with cells as rows and consistent column names, such that everything can be concatenated together
cov_out = pd.DataFrame(cov_out).transpose()
cov_out.to_csv(sample_name+"."+feat_name+".total_count.tsv.gz",sep="\t",header=True,index=True,compression="gzip")
"""					"""

"""	METHYLATION CG COVERAGE MATRIX
Z/(z+Z) that overlap per  feature"""
mc_chr = [cov_per_chrom(cov,args.feat,mc_cov_only=True) for cov in glob.glob(os.path.join(args.cov_folder,"chr*bed.gz")) if "chrY" not in cov ]
mc_out = pd.concat(mc_chr,axis=1).transpose().sort_index()
sample_name = args.cov_folder.split("/")[-1].split(".")[0]+"_"+args.cov_folder.split("/")[-1].split(".")[1]
mc_out.columns =[sample_name+x for x in mc_out.columns]
mc_out = pd.DataFrame(mc_out).transpose()
mc_out.to_csv(sample_name+"."+feat_name+".mc_count.tsv.gz",sep="\t",header=True,index=True,compression="gzip")
"""				"""


""" SIMPLE RATE ESTIMATE
met/all_c with a coverage cutoff for reporting (default is 10)"""
args=[(cov_out.T,mc_out.T,x,min_cg) for x in cov_out.index]
if __name__ == '__main__':
	# generate met coverage pandas df
	with multiprocessing.Pool(processes=1) as pool:
		out=pool.starmap(mcrate_simple,args)

df_rate = pd.DataFrame(out)
df_rate.to_csv(sample_name+"."+feat_name+".mc_simplerate.tsv.gz",sep="\t",header=True,index=True,compression="gzip")


""" NEG BINOM RATE ESTIMATE
#modified from https://github.com/lhqing/ALLCools/blob/master/ALLCools/mcds/utilities.py """
args=[(cov_out.T,mc_out.T,x,min_cg) for x in cov_out.index]
if __name__ == '__main__':
	# generate met coverage pandas df
	with multiprocessing.Pool(processes=1) as pool:
		out=pool.starmap(posterior_mcrate_estimate,args)

df_posterior = pd.DataFrame(out)
df_posterior.to_csv(sample_name+"."+feat_name+".mc_posteriorest.tsv.gz",sep="\t",header=True,index=True,compression="gzip")


#cd /volumes/seq/projects/metACT/240205_RMMM_scalebiotest2/cg_sort_cov
#for i in $(ls -d ./*/*); do python mc_cov_feature_summary.py --features genome_windows.100kb.bed --cov_folder $i ; done
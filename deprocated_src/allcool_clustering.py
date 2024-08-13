import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import anndata
import scanpy as sc
import dask
from ALLCools.clustering import \
    tsne, \
    significant_pc_test, \
    filter_regions, \
    remove_black_list_region, \
    lsi, \
    binarize_matrix
from ALLCools.plot import *
from ALLCools.mcds import MCDS
from pybedtools import BedTool
from ALLCools.mcds import RegionDS
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--metadata')
parser.add_argument('--mcds')
parser.add_argument('--allcool_path')
parser.add_argument('--outname')
parser.add_argument('--var_dim',default='chrom5k')
args = parser.parse_args()

metadata_path=args.metadata
mcds_path=args.mcds
allcool_path=args.allcool_path
var_dim=args.var_dim
outname=args.outname
outdir="."

#Set up metadata
metadata = pd.read_csv(metadata_path)
#correct BC names
repl_dict = {"\+":'_'}
metadata.replace({'BC':repl_dict}, inplace=True, regex=True)
metadata['cellid']=metadata['sampleName']+'.'+metadata['tgmt_well']+'.'+metadata['BC']
metadata.set_index('cellid',drop=True,inplace=True)
mcds = MCDS.open(mcds_path, var_dim=var_dim)
mcds.add_cell_metadata(metadata)


########Filter chromosomes##########
#https://lhqing.github.io/ALLCools/_modules/ALLCools/mcds/mcds.html#MCDS.remove_chromosome
exclude_chromosome = [ 'M','Y']
judge = mcds.coords[f"{var_dim}_chrom"].isin(exclude_chromosome)
print(f"{int(judge.sum())} {var_dim} features in {exclude_chromosome} removed.")
mcds = mcds.sel({f"{var_dim}": ~judge.to_numpy(),
    f"{var_dim}_chrom": ~judge.to_numpy(),
    f"{var_dim}_start": ~judge.to_numpy(),
    f"{var_dim}_end": ~judge.to_numpy()})

########Filter blacklist##########
# black_list_path='/ref/ENCFF356LFX.bed.gz'
# black_list_bed = BedTool(black_list_path)
# black_feature = feature_bed.intersect(black_list_bed, f=f, wa=True)
# test=mcds.coords.to_dataset()[[f"{var_dim}_chrom",f"{var_dim}_start",f"{var_dim}_end"]]

# bed_df = pd.DataFrame(
#     [chrom[chrom.columns[1]].to_series(),start[start.columns[1]],end[end.columns[1]]],
#     index=["chrom", "start", "end"],columns=mcds.get_index(var_dim)).T


# def check_dups(index, n):
#     s = set(index)
#     if len(s) != len(index):
#         print(f'df_list[{n}]:', index, end='\n\n')

# for n, df in enumerate(df_list):
#     check_dups(mcds.coords.columns, n)

# feature_bed_df = mcds.get_feature_bed(var_dim=var_dim)
# feature_bed = BedTool.from_dataframe(feature_bed_df)


######Clustering####    
#https://lhqing.github.io/ALLCools/api/ALLCools/mcds/mcds/index.html    
mcad = mcds.get_score_adata(mc_type='CGN', quant_type='hypo-score',sparse=False)
pc_cutoff = 0.1 # PC cutoff
knn = -1  # -1 means auto determine
resolution = 0.5 # Leiden

binarize_matrix(mcad, cutoff=0.95)
filter_regions(mcad)
lsi(mcad, algorithm='arpack', obsm='X_pca')
significant_pc_test(mcad, p_cutoff=pc_cutoff, update=True)
knn = max(15, int(np.log2(mcad.shape[0])*2))

sc.pp.neighbors(mcad, n_neighbors=knn)
sc.tl.leiden(mcad, resolution=resolution)

tsne(mcad,
     obsm='X_pca',
     metric='euclidean',
     exaggeration=-1,  # auto determined
     perplexity=30,
     n_jobs=-1)

sc.tl.umap(mcad)

#add some additional metadata
fig, axes = plt.subplots(3, 2, figsize=(6, 4), dpi=300)
ax=axes[0][0]
_ = categorical_scatter(data=mcad,ax=ax,
                        coord_base='tsne',
                        hue='leiden',
                        text_anno='leiden',
                        show_legend=False)
ax.set_title('TSNE')
ax.set_ylabel('Leiden')
ax=axes[0][1]
_ = categorical_scatter(data=mcad,ax=ax,
                        coord_base='umap',
                        hue='leiden',
                        text_anno='leiden',
                        show_legend=True)
ax.set_title('UMAP')
ax=axes[1][0]
_ = categorical_scatter(data=mcad,ax=ax,
                        coord_base='umap',
                        hue=mcds.cell_sampleName,
                        show_legend=True)
ax.set_title('SampleName')
ax=axes[1][1]
_ = categorical_scatter(data=mcad,ax=ax,
                        coord_base='umap',
                        hue=mcds.cell_tgmt_well,
                        show_legend=False)
ax.set_title('Tagmentation Well')

ax=axes[2][0]
_ = continuous_scatter(data=mcad,ax=ax,
                        coord_base='umap',
                        hue=mcds.cell_uniq)

ax.set_title('Unique Reads')
ax=axes[2][1]
_ = continuous_scatter(data=mcad,ax=ax,
                        coord_base='umap',
                        hue=mcds.cell_MitoReads)
ax.set_title('Mito Reads (Perc)')


outfilename=outname+'.mCG'+'_'+var_dim
plt.savefig(outdir+'/'+outfilename+'_clustering.png')


mcad.write_h5ad(outdir+'/'+outfilename+'-clustering.h5ad')
mcad.obs['allc_file']=[allcool_path+'/'+x+'/'+x+'.gz' for x in mcad.obs.index.to_list()]
mcad.obs.to_csv(outdir+'/'+outfilename+'.ClusteringResults.csv.gz',sep="\t")




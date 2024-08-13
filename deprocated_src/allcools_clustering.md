```bash
bsub -Is -W 36:00 -q long -n 10 -M 100 -R rusage[mem=100] /bin/bash

sif="/rsrch5/home/genetics/NAVIN_LAB/Ryan/projects/metact/src/allcool.sif"
module load singularity
singularity shell \
--bind /rsrch5/home/genetics/NAVIN_LAB/Ryan/projects/metact/ref:/ref/ \
--bind /rsrch5/home/genetics/NAVIN_LAB/Ryan/projects/metact/240205_RMMM_scalebiotest2 \
--bind /rsrch4/scratch/genetics/rmulqueen \
$sif 

projDir="/rsrch5/home/genetics/NAVIN_LAB/Ryan/projects/metact/240205_RMMM_scalebiotest2"
cd ${projDir}/postprocessing/bismark/sc_allcool


cat ${projDir}/report/HBCA-16R/csv/HBCA-16R.allCells.csv | grep "pass" > ${projDir}/cell_metadata.csv
cat ${projDir}/report/HBCA-83L/csv/HBCA-83L.allCells.csv | grep "pass" | tail -n +1 >> ${projDir}/cell_metadata.csv
```
#FILTER BLACKLIST REGIONS
#DEFINE GENE DMRS

# Cluster on 5k regions
```python
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
from ALLCools.dmr import call_dms, call_dmr
from ALLCools.clustering import one_vs_rest_dmg
from ALLCools.clustering import one_vs_rest_dmg

metadata_path="/rsrch5/home/genetics/NAVIN_LAB/Ryan/projects/metact/240205_RMMM_scalebiotest2/cell_metadata.csv"
mcds_path="/rsrch5/home/genetics/NAVIN_LAB/Ryan/projects/metact/240205_RMMM_scalebiotest2/postprocessing/bismark/sc_allcool/data.mcds"
var_dim='chrom5k'


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

outdir="/rsrch5/home/genetics/NAVIN_LAB/Ryan/projects/metact/240205_RMMM_scalebiotest2/postprocessing/bismark/sc_allcool"
plt.savefig(outdir+'/HBCA.mCG-5K-clustering.png')


mcad.write_h5ad('HBCA.mCG-5K-clustering.h5ad')
mcad.obs.to_csv('HBCA.ClusteringResults.csv.gz')

mcad.obs['allc_file']=['./'+x+'/'+x+'.gz' for x in mcad.obs.index.to_list()]
mcad.obs.to_csv('HBCA.ClusteringResults.csv.gz')

for i in mcad.obs.leiden.unique():
    out=mcad.obs[mcad.obs.leiden.isin([i])]['allc_file']
    out=pd.Series(out.to_list()) 
    out.to_csv('HBCA.Cluster_'+str(i)+'.allcool.tsv',header=False,index=False)
    merge_allc_files(allc_paths=out.to_list(),
                        output_path='HBCA.Cluster_'+str(i),
                        chrom_size_path='/ref/grch38/genome.txt',
                        cpu=5)

from ALLCools._merge_allc import *
import glob 
# Using '?' pattern 
for name in glob.glob('./HBCA.Cluster_*.allcool.tsv'): 
    print(name) 

cluster_col = 'leiden'
obs_dim = 'cell'
var_dim = 'transcripts'
mc_type = 'CGN'
coord_base = 'umap'
min_cov = 5

#Set up ensembl to gene id
gene_meta = pd.read_csv(f'/ref/grch38/filteredGTF/GRCh38_transcripts.longest.bed', sep="\t",names=['chrom','start','end','gene_id','gene_name']) #accidentally had gene length space separated on end
gene_meta['gene_name']=[x.split()[0] for x in gene_meta['gene_name']]
gene_name_to_gene_id = dict(zip(gene_meta['gene_name'],gene_meta['gene_id']))
gene_meta.index.name = 'gene_id'

def get_gene_values_by_name(gene_name):
    data = gene_frac_da.sel(transcripts=gene_name_to_gene_id[gene_name]).to_pandas()
    data.name = gene_name
    return data


cell_meta = mcad.obs.copy()

gene_mcds = MCDS.open(mcds_path, var_dim=var_dim)
gene_mcds.add_cell_metadata(cell_meta)
gene_mcds.add_cell_metadata(metadata)
gene_mcds.add_mc_rate(var_dim=var_dim,normalize_per_cell=True,clip_norm_value=10)
gene_mcds.add_mc_frac(var_dim=var_dim,normalize_per_cell=True, clip_norm_value=10)
gene_mcds.add_feature_cov_mean(var_dim=var_dim,)

genes_to_skip = set()
use_features = gene_meta.index
feature_cov_mean = gene_mcds.coords[f'{var_dim}_cov_mean'].to_pandas()
use_features = feature_cov_mean[feature_cov_mean > min_cov].index
print(f'{use_features.size} features remained')

gene_frac_da = gene_mcds[[f'transcripts_da_frac']]

#write out filtered gene list
gene_frac_da.write_dataset('geneslop2k_frac.mcds', var_dims=['geneslop2k'])
use_gene_meta = gene_meta.sel(use_features)
use_gene_meta.to_csv('GeneMetadata.csv.gz')


top_n = 100
auroc_cutoff = 0.8
adj_p_cutoff = 0.001
fc_cutoff = 0.5
max_cluster_cells = 2000
max_other_fold = 5
cpu = 10

dmg_table = one_vs_rest_dmg(cell_meta,
                            group=cluster_col,
                            mcds=gene_frac_da,
                            obs_dim=obs_dim,
                            var_dim=var_dim,
                            mc_type=mc_type,
                            top_n=top_n,
                            adj_p_cutoff=adj_p_cutoff,
                            fc_cutoff=fc_cutoff,
                            auroc_cutoff=auroc_cutoff,
                            max_cluster_cells=max_cluster_cells,
                            max_other_fold=max_other_fold,
                            cpu=cpu)
dmg_table.to_hdf(f'{cluster_col}.OneVsRestDMG.hdf', key='data')

#plot cluster DMGs
use_cells = gene_mcad.obs_names
cluster_dmg_path = f'{cluster_col}.OneVsRestDMG.hdf'

def plot_cluster_and_genes(cluster, cell_meta, cluster_col, genes_data,coord_base='umap', ncols=5, axes_size=3, dpi=150, hue_norm=(0.67, 1.5)):
    ncols = max(2, ncols)
    nrows = 1 + (genes_data.shape[1] - 1) // ncols + 1

    # figure
    fig = plt.figure(figsize=(ncols * axes_size, nrows * axes_size), dpi=dpi)
    gs = fig.add_gridspec(nrows=nrows, ncols=ncols)

    # cluster axes
    ax = fig.add_subplot(gs[0, 0])
    categorical_scatter(data=cell_meta,
                        ax=ax,
                        coord_base=coord_base,
                        axis_format=None,
                        hue=cluster_col,
                        palette='tab20')
    ax.set_title('All Clusters')
    ax = fig.add_subplot(gs[0, 1])
    categorical_scatter(data=cell_meta,
                        ax=ax,
                        coord_base=coord_base,
                        hue=cell_meta[cluster_col] == cluster,
                        axis_format=None,
                        palette={True: 'red',False: 'lightgray'})
    ax.set_title('This Cluster')

    # gene axes
    for i, (gene, data) in enumerate(genes_data.iteritems()):
        col = i % ncols
        row = i // ncols + 1
        ax = fig.add_subplot(gs[row, col])

        if ax.get_subplotspec().is_first_col() and ax.get_subplotspec().is_last_row():
            axis = 'tiny'
        else:
            axis = None

        continuous_scatter(ax=ax,
                           data=cell_meta,
                           hue=data,
                           axis_format=axis,
                           hue_norm=hue_norm,
                           coord_base=coord_base)
        ax.set_title(f'{data.name}')
    fig.suptitle(f'Cluster {cluster} Top Markers')
    return fig

cluster = 'c1'
genes = cluster_dmgs[cluster_dmgs['cluster'] == cluster].sort_values(
    'AUROC', ascending=False)[:10]
genes_data = gene_frac_da.sel(geneslop2k=genes.index).to_pandas()
genes_data.columns = genes_data.columns.map(gene_meta['gene_name'])

fig = plot_cluster_and_genes(cluster=cluster,
                             cell_meta=adata.obs,
                             cluster_col=cluster_col,
                             genes_data=genes_data,
                             coord_base='tsne',
                             ncols=5,
                             axes_size=3,
                             dpi=250,
                             hue_norm=(0.67, 1.5))

fig.savefig(f'{cluster}.TopMarker.png', bbox_inches='tight')

```

#collate clusters to higher coverage
```bash
for i in HBCA.Cluster*.allcool.tsv; do
    allcools merge-allc \
    --cpu 10 \
    --allc_paths $i \
    --output_path ${i::-12} \
    --chrom_size_path /ref/grch38/genome.txt
done

```
```python
from ALLCools.clustering import PairwiseDMG
from ALLCools.plot import *
import pandas as pd
import anndata
import matplotlib.pyplot as plt
from ALLCools.plot import *
from ALLCools.mcds import MCDS
import seaborn as sns
from scipy.stats import zscore

metadata_path="/rsrch5/home/genetics/NAVIN_LAB/Ryan/projects/metact/240205_RMMM_scalebiotest2/cell_metadata.csv"
mcds_path="/rsrch5/home/genetics/NAVIN_LAB/Ryan/projects/metact/240205_RMMM_scalebiotest2/postprocessing/bismark/sc_allcool/data.mcds"
cluster_col = 'leiden'
var_dim='transcripts'
mc_type='GCN'

mcad = anndata.read_h5ad('HBCA.mCG-5K-clustering.h5ad') #clustering info
cell_meta = mcad.obs.copy()

gene_mcds = MCDS.open(mcds_path, var_dim=var_dim)
gene_mcds.add_cell_metadata(cell_meta)
gene_mcds.add_mc_rate(var_dim=var_dim,normalize_per_cell=True,clip_norm_value=10)
gene_frac_da = gene_mcds[f'transcripts_da_frac']
gene_frac_da = gene_frac_da.sel().load()
gene_mcad = gene_mcds.get_adata()

#snRNA markers
markers={
"lumhr":["ANKRD30A","AFF3","ERBB4","TTC6","MYBPC1","NEK10","THSD4"],
"lumsec":["AC011247.1","COBL","GABRP","ELF5","CCL28","KRT15","KIT"],
"basal":["AC044810.2","CARMN","LINC01060","ACTA2","KLHL29","DST","IL1RAPL2"],
"fibro":["LAMA2","DCLK1","NEGR1","LINC02511","ANK2","KAZN","SLIT2"],
"lymphatic":["AL357507.1","PKHD1L1","KLHL4","LINC02147","RHOJ","ST6GALNAC3","MMRN1"],
"vascular":["MECOM","BTNL9","MCTP1","PTPRB","VWF","ADGRL4","LDB2"],
"perivasc":["RGS6","KCNAB1","COL25A1","ADGRL3","PRKG1","NR2F2-AS1","AC012409.2"],
"myeloid":["F13A1","MRC1","RBPJ","TBXAS1","FRMD4B","CD163","RAB31"],
"tcells":["SKAP1","ARHGAP15","PTPRC","THEMIS","IKZF1","PARP8","CD247"],
"mast":["NTM","IL18R1","SYTL3","SLC24A3","HPGD","TPSB2","HDC"],
"adipo":["PDE3B","ACACB","WDPCP","PCDH9","CLSTN2","ADIPOQ","TRHDE"]
}

#Set up ensembl to gene id
gene_meta = pd.read_csv(f'/ref/grch38/filteredGTF/GRCh38_transcripts.longest.bed', sep="\t",names=['chrom','start','end','gene_id','gene_name']) #accidentally had gene length space separated on end
gene_meta['gene_name']=[x.split()[0] for x in gene_meta['gene_name']]
gene_name_to_gene_id = dict(zip(gene_meta['gene_name'],gene_meta['gene_id']))
gene_meta.index.name = 'gene_id'

def get_gene_values_by_name(gene_name):
    data = gene_frac_da.sel(transcripts=gene_name_to_gene_id[gene_name]).to_pandas()
    data.name = gene_name
    return data


for key in markers.keys():
    celltype=key
    marker_genes=markers[key] #get marker genes
    marker_genes=[x for x in marker_genes if x in gene_meta['gene_name'].to_list()] #filter to those in our list
    if len(marker_genes)>0:
        nrows=len(marker_genes)+1
        ncols=2
        axes_size=5
        dpi=300
        fig = plt.figure(figsize=(ncols*axes_size,nrows*axes_size), dpi=dpi)
        gs = fig.add_gridspec(nrows=nrows, ncols=ncols)
        ax = fig.add_subplot(gs[0, 0])
        _ = categorical_scatter(data=mcad,ax=ax,
                        coord_base='umap',
                        hue='leiden',palette='tab20b',
                        text_anno='leiden',
                        show_legend=True)
        for i in range(1,nrows):
            gene=marker_genes[i-1]
            row=i
            col=0
            ax = fig.add_subplot(gs[row, col])
            continuous_scatter(data=mcad,
                            ax=ax,
                            coord_base='umap',
                            axis_format=None,
                            hue=get_gene_values_by_name(gene))
            ax.set_ylabel(f'{gene}')
            gene_dat=pd.merge(get_gene_values_by_name(gene),cell_meta.leiden,left_index=True,right_index=True)
            gene_dat[gene]=zscore(gene_dat[gene])
            col=1
            ax = fig.add_subplot(gs[row, col])
            sns.boxenplot(data=gene_dat, 
                            x=cluster_col,hue=cluster_col,legend=False, y=gene, 
                            ax=ax,palette='tab20b')
    fig.savefig(f'{celltype}.Markers.png', bbox_inches='tight')


```

Scan for TFmotifs
```python
from ALLCools.mcds import RegionDS
from ALLCools.motif import MotifSet, get_default_motif_set
from ALLCools.motif.parse_meme import _parse_meme_database as parse_meme_database
from ALLCools import __path__
import pandas as pd
# check out the default motif set
#https://lhqing.github.io/ALLCools/cluster_level/RegionDS/04.motif_scan.html
#default_motif_set = get_default_motif_set() #broken due to removal of pandas squeeze

PACKAGE_DIR = __path__[0]
DEFAULT_MOTIF_DIR = f"{PACKAGE_DIR}/motif/default_motif_set/"
def get_default_motif_set(database="three_databases"):
    if database == "three_databases":
        motif_set = parse_meme_database(
            f"{DEFAULT_MOTIF_DIR}/JASPAR2018HOCOMOCOv11Jolma2013.meme",
            f"{DEFAULT_MOTIF_DIR}/JASPAR2018HOCOMOCOv11Jolma2013.metadata.csv",
        )        # default thresholds, fnr/fpr = 1000
        motif_set.thresholds = pd.read_csv(
            f"{DEFAULT_MOTIF_DIR}/JASPAR2018HOCOMOCOv11Jolma2013.thresholds.csv",
            header=None,
            index_col=0,
        ).squeeze('columns').to_dict()
        return motif_set
    else:
        # TODO: allow user create motif set by providing meme and metadata (optional)
        raise NotImplementedError

default_motif_set = get_default_motif_set() 
default_motif_set.n_motifs

dmr_ds = RegionDS.open('test_HIP', select_dir=['dmr'])

```











































# Cluster on CGI
```python
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


metadata_path="/rsrch5/home/genetics/NAVIN_LAB/Ryan/projects/metact/240205_RMMM_scalebiotest2/cell_metadata.csv"
mcds_path="/rsrch5/home/genetics/NAVIN_LAB/Ryan/projects/metact/240205_RMMM_scalebiotest2/postprocessing/bismark/sc_allcool/data.mcds"
var_dim='CGI'


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


######Clustering####    
#https://lhqing.github.io/ALLCools/api/ALLCools/mcds/mcds/index.html    
mcad = mcds.get_score_adata(mc_type='CGN', quant_type='hypo-score',sparse=False)
pc_cutoff = 0.1 # PC cutoff
knn = -1  # -1 means auto determine
resolution = 1 # Leiden

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
#mcad.obs['cell_tgmt_well'] = mcds.variables['cell_tgmt_well'].to_dataframe()[1].to_list()
#mcad.obs['cell_uniq'] = mcds.variables['cell_uniq'].to_dataframe()[1].to_list()

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




def get_gene_values_by_name(gene_name):
    data = gene_frac_da.sel(transcripts=gene_name_to_gene_id[gene_name]).to_pandas()
    data.name = gene_name
    return data


gene = 'KRT15'
hue_norm = (0.67, 1.5)
coord_base = 'umap'

_ = continuous_scatter(ax=ax,
                       data=mcad,
                       hue=get_gene_values_by_name(gene),
                       coord_base='tsne',
                       max_points=None,
                       s=4)

ax.set_title('UMAP')
outdir="/rsrch5/home/genetics/NAVIN_LAB/Ryan/projects/metact/240205_RMMM_scalebiotest2/postprocessing/bismark/sc_allcool"
plt.savefig(outdir+'/HBCA.mCG-CGI-clustering.png')


# axes = categorical_scatter(data=mcad,ax=axes[1][1],
#                         coord_base='umap',
#                         hue='cell_tgmt_well',
#                         text_anno='cell_tgmt_well',
#                         show_legend=False)

mcad.write_h5ad('HBCA.mCG-5K-clustering.h5ad') #adata = anndata.read_h5ad(adata_path)

mcad.obs.to_csv('HBCA.ClusteringResults.csv.gz')

```

#### Differentially Methylated Genes ####
Using 5k regions
```python
from ALLCools.clustering import PairwiseDMG
from ALLCools.plot import *
import pandas as pd
import anndata
import matplotlib.pyplot as plt
from ALLCools.plot import *
from ALLCools.mcds import MCDS


metadata_path="/rsrch5/home/genetics/NAVIN_LAB/Ryan/projects/metact/240205_RMMM_scalebiotest2/cell_metadata.csv"
mcds_path="/rsrch5/home/genetics/NAVIN_LAB/Ryan/projects/metact/240205_RMMM_scalebiotest2/postprocessing/bismark/sc_allcool/data.mcds"
var_dim='CGI'

cluster_col = 'leiden'
var_dim='transcripts'
mc_type='GCN'
cell_meta = mcad.obs.copy()

mcad.write_h5ad('HBCA.mCG-5K-clustering.h5ad')

gene_mcds = MCDS.open(mcds_path, var_dim=var_dim)
gene_mcds.add_cell_metadata(cell_meta)
gene_mcds.add_mc_rate(var_dim=var_dim,normalize_per_cell=True,clip_norm_value=10)
gene_frac_da = gene_mcds[f'transcripts_da_frac']
gene_frac_da = gene_frac_da.sel().load()
gene_mcad = gene_mcds.get_adata()

# mc_type = 'CGN'
# top_n = 1000
# adj_p_cutoff = 1e-3
# delta_rate_cutoff = 0.3
# auroc_cutoff = 0.9
# random_state = 0
# n_jobs = 10

# pwdmg = PairwiseDMG(max_cell_per_group=1000,
#                     top_n=top_n,
#                     adj_p_cutoff=adj_p_cutoff,
#                     delta_rate_cutoff=delta_rate_cutoff,
#                     auroc_cutoff=auroc_cutoff,
#                     random_state=random_state,
#                     n_jobs=n_jobs)

# pwdmg.fit_predict(x=gene_mcds[f'{var_dim}_da_frac'], 
#                   var_dim=var_dim,groups=cell_meta[cluster_col])

gene_meta = pd.read_csv(f'/ref/grch38/filteredGTF/GRCh38_transcripts.longest.bed', sep="\t",names=['chrom','start','end','gene_id','gene_name']) #accidentally had gene length space separated on end
gene_meta['gene_name']=[x.split()[0] for x in gene_meta['gene_name']]
gene_name_to_gene_id = dict(zip(gene_meta['gene_name'],gene_meta['gene_id']))
gene_meta.index.name = 'gene_id'

def get_gene_values_by_name(gene_name):
    data = gene_frac_da.sel(transcripts=gene_name_to_gene_id[gene_name]).to_pandas()
    data.name = gene_name
    return data


gene = 'KRT15'
hue_norm = (0.67, 1.5)
coord_base = 'umap'

fig, ax = plt.subplots(figsize=(4, 4), dpi=300)
_ = continuous_scatter(ax=ax,
                       data=mcad,
                       hue=get_gene_values_by_name(gene),
                       coord_base=coord_base,
                       max_points=None,
                       s=4)


plt.savefig(outdir+'/HBCA.mCG-5K-clustering.'+gene+'.png')

```
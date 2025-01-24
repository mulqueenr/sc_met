
```bash
singularity shell \
--bind ~/projects \
--bind /volumes/seq/projects/metACT \
--bind /volumes/seq/projects/gccACT \
~/singularity/copykit.sif
source ~/.bashrc
```

# Reprocessing of all scalebio data generated thus far. 
```R

setwd("/volumes/USR2/Ryan/projects/metact/cnv_processing") #set wd
source("~/projects/metact/src/copykit_custom_functions.R") #to load in custom functions and libraries
cpu_count=300
register(MulticoreParam(progressbar = T, workers = cpu_count), default = T)

act_dir<- "/volumes/seq/projects/gccACT/mdamb231_ACTseq/cells"
cellline_dir<- "/volumes/USR2/Ryan/projects/metact/240205_RMMM_scalebiotest2/postprocessing/sc_bam"
dcis1_dir<-"/volumes/USR2/Ryan/projects/metact/240526_RMMM_scalebio_dcis/sc_bams"
dcis_79t_dir <- "/volumes/USR2/Ryan/projects/metact/241007_RM_scalebio_dcis2/alignments/dedup/DCIS-79T_split_bam"
dcis_92t_dir <- "/volumes/USR2/Ryan/projects/metact/241007_RM_scalebio_dcis2/alignments/dedup/DCIS-92T_split_bam"
hbca_17_dir <- "/volumes/USR2/Ryan/projects/metact/241007_RM_scalebio_dcis2/alignments/dedup/HBCA-17T_split_bam"
hbca_19_dir <- "/volumes/USR2/Ryan/projects/metact/241007_RM_scalebio_dcis2/alignments/dedup/HBCA-19T_split_bam"
idc_79t_dir <- "/volumes/USR2/Ryan/projects/metact/241007_RM_scalebio_dcis2/alignments/dedup/IDC-79T_split_bam"
dcis_66t_dir<-"/volumes/USR2/Ryan/projects/metact/241007_RM_scalebio_dcis2/alignments/dedup/DCIS-66T_split_bam"

#cell line processing and comparison
    act<-generate_copykit_object(indir=act_dir,info='act',experiment_name="gccact",sample="mda-mb-231")#act mdamb231
    cellline<-generate_copykit_object(indir=cellline_dir, experiment_name="cellline_scalebio",sample="cellline")#cellline
    colData(cellline)$sample_name<-tolower(unlist(lapply(strsplit(colnames(cellline),"[.]"),"[",1))) #rename to be cell line specific after generation
    merged_cellline<-cbind(act,cellline[,colData(cellline)$sample_name %in% c("mcf10a","mcf7","mda-mb-231")])
    merged_cellline<-process_copykit_object(dat=merged_cellline,prefix="cellline")

#dcis set 1
    dcis1 <- generate_copykit_object(indir=dcis1_dir,experiment_name="dcis1_scalebio",sample="dcis")
    colData(dcis1)$sample_name<-tolower(unlist(lapply(strsplit(colnames(dcis1),"[.]"),"[",1)))
    dcis_41t<-dcis1[,colData(dcis1)$sample_name %in% c("dcis-41t")]
    #split and run each separately
    for (i in unique(dcis1$sample_name)){
        if(sum(colData(dcis1)$sample_name == i)>50){
        dcis1_sample<-dcis1[,colData(dcis1)$sample_name == i ]
        cellline<-process_copykit_object(dat=dcis1_sample,prefix=i)
        }
    }

#dcis set 2
    #dcis_66t
    dcis_66t <- generate_copykit_object(indir=dcis_66t_dir,experiment_name="dcis2_scalebio",sample="dcis_66t")
    process_copykit_object(dat=dcis_66t,prefix="scalebio_dcis66t",min_readcount=250000)

    #79t
    dcis_79t <- generate_copykit_object(indir=dcis_79t_dir,experiment_name="dcis2_scalebio",sample="dcis_79t")
    idc_79t <- generate_copykit_object(indir=idc_79t_dir,experiment_name="dcis2_scalebio",sample="idc_79t")
    dcis_idc_79t<-cbind(dcis_79t,idc_79t)
    process_copykit_object(dat=dcis_idc_79t,prefix="scalebio_idcdcis79t",min_readcount=250000)

    #92t
    dcis_92t<- generate_copykit_object(indir=dcis_92t_dir,experiment_name="dcis2_scalebio",sample="dcis_92t")
    process_copykit_object(dat=dcisc_92t,prefix="scalebio_dcis92t",min_readcount=250000)

    #hbca 17
    hbca_17<- generate_copykit_object(indir=hbca_17_dir,experiment_name="dcis2_scalebio",sample="hbca_17")
    process_copykit_object(dat=hbca_17,prefix="scalebio_hbca_17",min_readcount=250000)

    #hbca 19
    hbca_19<- generate_copykit_object(indir=hbca_19_dir,experiment_name="dcis2_scalebio",sample="hbca_19")
    process_copykit_object(dat=hbca_19,prefix="scalebio_hbca_19",min_readcount=250000)

#all hbca
    
    hbca_cnv<-cellline[,colData(cellline)$sample_name %in% c("hbca-83l")]
    process_copykit_object(dat=hbca_cnv,prefix="scalebio_hbca_83l",min_readcount=250000)


    hbca_cnv<-cbind(hbca_17,hbca_19,
            cellline[,colData(cellline)$sample_name %in% c("hbca-16r","hbca-83l")])
    process_copykit_object(dat=hbca_cnv,prefix="scalebio_hbca_all",min_readcount=250000)
```

DMR Analysis
```
DCIS-66T.amethyst.rds
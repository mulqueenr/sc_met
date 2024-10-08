# You'll need a docker login and a AWS login.

Follow this to make a EC2 instance.
https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/get-set-up-for-amazon-ec2.html

### Set these settings for the EC2 instance.
- Make with ubuntu
- Free tier 64bit (x86)
- t3.xlarge
- Generated pem
- 100gb gp2 storage

Change permissions for pem key-pair
```bash
chmod 400 ~/Downloads/newkey2.pem
```

ssh -i <path to pem> ubuntu@<public ip4>
```bash
ssh -i ~/Downloads/newkey2.pem ubuntu@54.187.193.117
```


## Install singularity, go, and docker. Use docker to pull the ubuntu image
https://docs.sylabs.io/guides/latest/user-guide/quick_start.html
```bash
#install basics ##
# Ensure repositories are up-to-date
sudo apt-get update
# Install debian packages for dependencies
sudo apt-get install -y \
   autoconf \
   automake \
   cryptsetup \
   git \
   libfuse-dev \
   libglib2.0-dev \
   libseccomp-dev \
   libtool \
   pkg-config \
   runc \
   squashfs-tools \
   squashfs-tools-ng \
   uidmap \
   wget \
   zlib1g-dev \
   make

## Install go ##
export VERSION=1.21.0 OS=linux ARCH=amd64 && \
  wget https://dl.google.com/go/go$VERSION.$OS-$ARCH.tar.gz && \
  sudo tar -C /usr/local -xzvf go$VERSION.$OS-$ARCH.tar.gz && \
  rm go$VERSION.$OS-$ARCH.tar.gz

echo 'export PATH=/usr/local/go/bin:$PATH' >> ~/.bashrc && \
  source ~/.bashrc

## Install singularity ##
export VERSION=4.1.0 && \
    wget https://github.com/sylabs/singularity/releases/download/v${VERSION}/singularity-ce-${VERSION}.tar.gz && \
    tar -xzf singularity-ce-${VERSION}.tar.gz && \
    cd singularity-ce-${VERSION}

./mconfig && \
    make -C builddir && \
    sudo make -C builddir install

## Install docker and use to pull the ubuntu image. ##
#Convert ubuntu image to sif file for building our own SIFs
cd 
curl -fsSL https://get.docker.com -o get-docker.sh
sudo sh ./get-docker.sh
sudo chmod 666 /var/run/docker.sock

#Pull ubuntu image to act as our base ##
#docker pull ubuntu:latest
#docker create --name ubuntu -p 80:80 ubuntu:latest
#docker images #take IMAGE ID 
#sudo docker save fd1d8f58e8ae -o ubuntu.tar #save image as tar
#sudo singularity build ubuntu.sif docker-archive://ubuntu.tar #convert tar to sif to act as local bootstrap

```

```bash
sudo singularity build --sandbox copykit/ docker://ubuntu:latest
sudo singularity shell --writable copykit/
```

## Define the copykit image to be created

copykit.def
```bash
Bootstrap: docker
From: ubuntu:latest

%environment
	# set up all essential environment variables
	export LC_ALL=C
	export PATH=/opt/miniconda3/bin:$PATH
	export PYTHONPATH=/opt/miniconda3/lib/python3.9/:$PYTHONPATH
	export LC_ALL=C.UTF-8
%post
	# update and install essential dependencies
	apt-get -y update
	apt-get update && apt-get install -y automake \
	build-essential \
	bzip2 \
	wget \
	git \
	default-jre \
	unzip \
	zlib1g-dev \
	parallel \
	libglpk40 \
	gfortran

	# download, install, and update miniconda3
	wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
	bash Miniconda3-latest-Linux-x86_64.sh -b -f -p /opt/miniconda3/
	rm Miniconda3-latest-Linux-x86_64.sh

	# install dependencies via conda
	export PATH="/opt/miniconda3/bin:$PATH"
	conda install -y -c conda-forge mamba 
	mamba install -y -f bioconda::samtools #
	mamba install -y -f bioconda::bedtools #
	mamba install -y -f conda-forge::parallel #

	#install R packages
	mamba install -y -f conda-forge::r-base #=4.2
	mamba install -y -f conda-forge::r-devtools
	#mamba install -y -f bioconda::bioconductor-ensdb.hsapiens.v86

	R --slave -e 'devtools::install_github("navinlabcode/copykit")'
	conda install -y -f --no-deps conda-forge::r-igraph
	conda install -y -f --no-deps bioconda::bioconductor-bluster
	conda install -y -f --no-deps bioconda::bioconductor-copynumber
	conda install -y -f --no-deps bioconda::bioconductor-ggtree
	wget https://github.com/navinlabcode/copykit/releases/download/v.0.1.2/copykit_0.1.2.tar.gz
	R --slave -e 'install.packages("copykit_0.1.2.tar.gz", repos = NULL)' # the install_github is broken so pulling from archive

	R --slave -e 'install.packages("optparse", repos="http://cran.us.r-project.org")' #


%labels
	Author Ryan Mulqueen
	Version v0.1
	MyLabel Copykit 

```

```bash
sudo singularity build archr_signac.sif archr_signac.def

sudo singularity build --sandbox archr_signac/ docker://ubuntu:latest
sudo singularity shell --writable archr_signac/

```

## Define the archR/Signac image to be created

archr_signac.def
```bash
Bootstrap: docker
From: ubuntu:latest

%environment
	# set up all essential environment variables
	export LC_ALL=C
	export PATH=/opt/miniconda3/bin:$PATH
	export PYTHONPATH=/opt/miniconda3/lib/python3.9/:$PYTHONPATH
	export LC_ALL=C.UTF-8
%post
	# update and install essential dependencies
	apt-get -y update
	apt-get update && apt-get install -y automake \
	build-essential \
	bzip2 \
	wget \
	git \
	default-jre \
	unzip \
	zlib1g-dev \
	parallel \
	libglpk40 \
	gfortran

	# download, install, and update miniconda3
	wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
	bash Miniconda3-latest-Linux-x86_64.sh -b -f -p /opt/miniconda3/
	rm Miniconda3-latest-Linux-x86_64.sh

	# install dependencies via conda
	export PATH="/opt/miniconda3/bin:$PATH"
	conda install -y -c conda-forge mamba 

	#install R packages
	mamba install -y -f conda-forge::r-base=4.1
	mamba install -y -f conda-forge::r-devtools
	mamba install -y -f conda-forge::r-essentials
	#mamba install -y -f bioconda::bioconductor-bsgenome.hsapiens.ucsc.hg38 #
	#mamba install -y -f bioconda::bioconductor-ensdb.hsapiens.v86 #

	R --slave -e 'install.packages("BiocManager", repos="http://cran.us.r-project.org")'
	R --slave -e 'BiocManager::install("biovizBase")'
	R --slave -e 'devtools::install_github("GreenleafLab/ArchR", ref="master", repos = BiocManager::repositories())'
	R --slave -e 'install.packages("Signac", repos="http://cran.us.r-project.org")'
	R --slave -e 'install.packages("Seurat", repos="http://cran.us.r-project.org")'
	R --slave -e 'install.packages("stringr", repos="http://cran.us.r-project.org")'
	R --slave -e 'devtools::install_github("swaruplabUCI/ArchRtoSignac", repos = BiocManager::repositories())' #
	R --slave -e 'BiocManager::install("EnsDb.Hsapiens.v86")'

	#library(ArchRtoSignac)
	#packages <- c("ArchR","Seurat", "Signac","stringr","EnsDb.Hsapiens.v86") # required packages
	#loadinglibrary(packages)


%labels
	Author Ryan Mulqueen
	Version v0.1
	MyLabel ArchR_Signac

```

```bash
sudo singularity build copykit.sif copykit.def

rsync -alPvz \
/home/ubuntu/copykit.sif \
mulqueen@qcprpgeo.mdanderson.edu

#sudo singularity build copykit.sif copykit.def 

#sudo singularity build copykit.sif copykit.def
#sudo singularity shell copykit.sif
```

## Generat scmet analysis sif

```bash
sudo singularity build --sandbox scmetR/ docker://ubuntu:latest
sudo singularity shell --writable scmetR/
#test r installs and stuff (just go line by line of def below)

```

scmetR.def
```bash
Bootstrap: docker
From: ubuntu:latest

%environment
# set up all essential environment variables
export LC_ALL=C
export PATH=/opt/miniconda3/bin:$PATH
export PYTHONPATH=/opt/miniconda3/lib/python3.9/:$PYTHONPATH

%post
# update and install essential dependencies
apt-get -y update
apt-get update && apt-get install -y automake \
build-essential \
bzip2 \
wget \
git \
default-jre \
unzip \
zlib1g-dev \
parallel

# download, install, and update miniconda3
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b -f -p /opt/miniconda3/
rm Miniconda3-latest-Linux-x86_64.sh

# install dependencies via conda
export PATH="/opt/miniconda3/bin:$PATH"
conda config --add channels bioconda
conda config --add channels conda-forge
conda install -y -c conda-forge mamba # general dependencies

#mamba installs
mamba install -y -f pip #
pip install pybedtools
pip install pandas
pip install scipy
mamba install -y -f numpy #
mamba install -y -f bioconda::bwa #
mamba install -y -f bioconda::samtools #
mamba install -y -f bioconda::bedtools #
mamba install -y -f conda-forge::parallel #

#install R packages
mamba install -y -f conda-forge::r-base #
mamba install -y -f conda-forge::r-devtools #
mamba install -y -f conda-forge::icu
mamba install -y -f bioconda::bioconductor-biocgenerics 
mamba install -y -f bioconda::bioconductor-sparsearray
mamba install -y -f bioconda::bioconductor-s4arrays
mamba install -y -f conda-forge::r-biocmanager #

#Funner stuff!
R --slave -e 'install.packages("Seurat", repos="http://cran.us.r-project.org")' #
R --slave -e 'devtools::install_github("stuart-lab/signac", "develop")' #
R --slave -e 'remotes::install_github("satijalab/seurat-wrappers")' #
R --slave -e 'install.packages(c("DescTools", "reshape2", "ggridges", "mice"), repos="http://cran.us.r-project.org")' #

mamba install -y -f conda-forge::r-rlang #
mamba install -y -f conda-forge::r-ggplot2 #
mamba install -y -f bioconda::bioconductor-dirichletmultinomial #
mamba install -y -f conda-forge::r-igraph #
mamba install -y -f conda-forge::r-rjags #
mamba install -y -f conda-forge::r-leiden #
mamba install -y -f conda-forge::r-hdf5r #
mamba install -y -f conda-forge::r-rmpfr #
mamba install -y -f conda-forge::r-ggraph #
mamba install -y -f conda-forge::r-nloptr #
mamba install -y -f conda-forge::r-jomo #

#R utility libraries
R --slave -e 'install.packages("remotes", repos="http://cran.us.r-project.org")' #
R --slave -e 'install.packages("circlize", repos="http://cran.us.r-project.org")' #
R --slave -e 'install.packages("optparse", repos="http://cran.us.r-project.org")' #
R --slave -e 'install.packages("patchwork", repos="http://cran.us.r-project.org")' #
R --slave -e 'install.packages("plyr", repos="http://cran.us.r-project.org")' #
R --slave -e 'install.packages("stringr", repos="http://cran.us.r-project.org")' #
R --slave -e 'install.packages("tidyverse", repos="http://cran.us.r-project.org")' #
R --slave -e 'install.packages("RColorBrewer", repos="http://cran.us.r-project.org")' #

#Bioconductor packages through conda
mamba install -y -f bioconda::bioconductor-biocparallel #
mamba install -y -f bioconda::bioconductor-bsgenome.hsapiens.ucsc.hg38 #
mamba install -y -f bioconda::bioconductor-ensdb.hsapiens.v86 #
mamba install -y -f bioconda::bioconductor-genomicranges
mamba install -y -f bioconda::bioconductor-jaspar2020 #
mamba install -y -f bioconda::bioconductor-org.hs.eg.db #
mamba install -y -f bioconda::bioconductor-tfbstools #
mamba install -y -f bioconda::bioconductor-txdb.hsapiens.ucsc.hg38.knowngene #
mamba install -y -f bioconda::bioconductor-universalmotif #
mamba install -y -f bioconda::bioconductor-chromvar #
mamba install -y -f bioconda::bioconductor-motifmatchr #
mamba install -y -f bioconda::bioconductor-scran #
mamba install -y -f bioconda::bioconductor-complexheatmap #
mamba install -y -f bioconda::bioconductor-biovizbase #

#correct cistopic install and matrix install
cd #
wget https://github.com/aertslab/cisTopic/archive/refs/tags/v2.1.0.tar.gz 
R --slave -e 'install.packages("v2.1.0.tar.gz", repos = NULL)' # the install_github is broken so pulling from archive

#Correct matrix version
mamba install -y -f conda-forge::r-matrix=1.6_1 # set this version
mamba install -y -f r::r-irlba # set this version
R --slave -e 'install.packages("Matrix", type = "source",repos="http://cran.us.r-project.org")' #reinstall from source
R --slave -e 'install.packages("irlba", type = "source",repos="http://cran.us.r-project.org")' #reinstall from source
R --slave -e 'install.packages("SeuratObject", type = "source",repos="http://cran.us.r-project.org")' #reinstall from source
R --slave -e 'oo <- options(repos = "https://cran.r-project.org/");
	tools::package_dependencies("Matrix", which = "LinkingTo", reverse = TRUE)[[1L]];
	install.packages("lme4", type = "source",repos = "http://cran.us.r-project.org");
	options(oo)'

%labels
    Author Ryan Mulqueen
    Version v0.0
    MyLabel Singlecell Methylation

```

```bash
sudo singularity build scmetR.sif scmetR.def
sudo singularity shell scmetR.sif
```

allcool.def
```bash
Bootstrap: docker
From: ubuntu:latest

%environment
# set up all essential environment variables
export LC_ALL=C
export PATH=/opt/miniconda3/bin:$PATH
export PYTHONPATH=/opt/miniconda3/lib/python3.9/:$PYTHONPATH
export PYTHONPATH=/opt/miniconda3/lib/python3.12/:$PYTHONPATH
source /opt/etc/bashrc

%post
# update and install essential dependencies
apt-get -y update
apt-get update && apt-get install -y automake \
build-essential \
bzip2 \
wget \
git \
default-jre \
unzip \
zlib1g-dev \
parallel \
nano

# download, install, and update miniconda3
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b -f -p /opt/miniconda3/
rm Miniconda3-latest-Linux-x86_64.sh

# install dependencies via conda
export PATH="/opt/miniconda3/bin:$PATH"
conda config --add channels bioconda
conda config --add channels conda-forge
conda install -y -c conda-forge mamba # general dependencies

#mamba installs
mamba install -y -f pip #
pip install pybedtools
pip install pandas
pip install scipy
mamba install -y -f numpy #
mamba install -y -f bioconda::bwa #
mamba install -y -f bioconda::samtools #
mamba install -y -f bioconda::bedtools #
mamba install -y -f conda-forge::parallel #

#from yap env yaml
mamba create -y -n bismark_env \
python=3.7 \
cutadapt=2.10 \
bismark \
picard \
subread=2.0 \
bowtie2=2.3 \
bowtie=1.3 \
matplotlib

#from allcool env yaml
mamba create -y -n allcool_env \
python=3.8 \
pip \
anndata \
biopython \
dask \
numba \
htslib>=1.9 \
jupyter \
jupyter_contrib_nbextensions \
leidenalg \
natsort \
netCDF4 \
networkx \
opentsne \
plotly \
pybedtools \
pyBigWig \
pynndescent \
pysam \
pytables \
scanpy \
scikit-learn=1.2.2 \
seaborn \
statsmodels \
xarray \
yaml \
zarr \
&& echo -e "#! /bin/bash\n\n# script to activate the conda environment" > ~/.bashrc \
&& conda init bash \
&& echo -e "\nsource activate allcool_env" >> ~/.bashrc \
&& conda clean -y -a \
&& mkdir -p /opt/etc/bashrc \
&& cp ~/.bashrc /opt/etc/bashrc
conda install numba
pip install anndata scanpy ALLCools xarray pybedtools dask leidenalg plotly imblearn openTSNE
pip install papermill imblearn 
pip install zarr
pip install dask
pip install numba
pip install matplotlib==3.6
conda install numba
pip install tables
pip install ctxcore
pip install psutil
pip install biopython
pip install scikit-learn==1.2.2


%labels
    Author Ryan Mulqueen
    Version v0.0
    MyLabel AllCools SIF

```

```bash
singularity build --fakeroot allcool.sif allcool.def #trying on geo
sudo singularity shell allcool.sif
```

## Pull Images
Use sftp to get images off cluster. Move to local computer > geo > seadragon2

```bash
sftp -i ~/Downloads/newkey2.pem ubuntu@54.187.193.117
get copykit.sif

sftp mulqueen@qcprpgeo.mdanderson.edu
cd ./singularity

put copykit.sif
```

Now move to seadragon
```bash
ssh Rmulqueen@seadragon2
bsub -Is -W 6:00 -q transfer -n 1 -M 16 -R rusage[mem=16] /bin/bash #get interactive transfer node this has internet access for environment set up
lcd /rsrch4/home/genetics/rmulqueen/singularity
get scmetR.sif

```

```bash
#singularity shell --bind /rsrch4/home/genetics/rmulqueen/singularity/scmetR.sif
```



## Define the amethyst image to be created

```bash
Bootstrap: docker
From: ubuntu:latest

%environment
	# set up all essential environment variables
	export LC_ALL=C
	export PATH=/opt/miniconda3/bin:$PATH
	export PYTHONPATH=/opt/miniconda3/lib/python3.9/:$PYTHONPATH
	export PYTHONPATH=/opt/miniconda3/lib/python3.12/:$PYTHONPATH
	source /opt/etc/bashrc

%post
	# update and install essential dependencies
	apt-get -y update
	apt-get update && apt-get install -y automake \
	build-essential \
	bzip2 \
	wget \
	git \
	default-jre \
	unzip \
	parallel \
	nano

	# download, install, and update miniconda3
	wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
	bash Miniconda3-latest-Linux-x86_64.sh -b -f -p /opt/miniconda3/
	rm Miniconda3-latest-Linux-x86_64.sh

	# install dependencies via conda
	export PATH="/opt/miniconda3/bin:$PATH"
	conda install -y -c conda-forge mamba 
	mamba install -y -c conda-forge r

	#mamba installs
	mamba install -y -f pip #
	pip install pybedtools
	pip install pandas
	pip install scipy
	conda install -y -f numpy bioconda::bwa bioconda::samtools bioconda::bedtools conda-forge::parallel #

	#install python libraries
	pip install h5py os numpy argparse sys

	#install R packages
	R --slave -e 'install.packages("Seurat",repos="http://cran.us.r-project.org")'
	R --slave -e 'install.packages("Signac",repos="http://cran.us.r-project.org")'
	R --slave -e 'install.packages("BiocManager", repos="http://cran.us.r-project.org")'
	R --slave -e 'BiocManager::install(c("caret", "devtools", "data.table", "dplyr", "furrr", "future", "future.apply",
	"ggplot2", "grDevices", "gridExtra", "igraph", "irlba", "janitor", "methods", 
	"plotly", "plyr", "purrr", "randomForest", "rhdf5", "rtracklayer", "scales", "stats", "stringr", 
	"tibble", "tidyr", "umap", "utils"))'
	R --slave -e 'devtools::install_github("JinmiaoChenLab/Rphenograph")'
	R --slave -e 'devtools::install_github("KrishnaswamyLab/MAGIC/Rmagic")'
	R --slave -e 'devtools::install_github("lrylaarsdam/amethyst", ref = "dev")'
	wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.annotation.gtf.gz

```


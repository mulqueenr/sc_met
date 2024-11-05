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
chmod 400 <path to pem>
```

ssh -i <path to pem> ubuntu@<public ip4>

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

# Build amethyst sandbox
```bash
sudo singularity build --sandbox amethyst2/ docker://ubuntu:latest
sudo singularity shell --writable amethyst2/
# %environment
# set up all essential environment variables
export LC_ALL=C
export PATH=/opt/miniconda3/bin:$PATH
export PYTHONPATH=/opt/miniconda3/lib/python3.9/:$PYTHONPATH
export PYTHONPATH=/opt/miniconda3/lib/python3.12/:$PYTHONPATH

# %post
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
nano \
libfontconfig1-dev \
libtiff-dev

# download, install, and update miniconda3
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b -f -p /opt/miniconda3/
rm Miniconda3-latest-Linux-x86_64.sh

# install dependencies via conda
export PATH="/opt/miniconda3/bin:$PATH"
conda config --add channels bioconda
conda config --add channels conda-forge
conda install -y -c conda-forge mamba # general dependencies

#mamba installs pip
mamba install -y -f pip

#pip installs for premethyst
pip install h5py numpy argparse pybedtools pandas scipy
mamba install -y -f -c bioconda bwa samtools bedtools
mamba install -y -f -c conda-forge parallel
mamba install -y -f -c conda-forge r r-devtools

#install R packages
R --slave -e 'install.packages("Seurat",repos="http://cran.us.r-project.org")'
R --slave -e 'install.packages("BiocManager", repos="http://cran.us.r-project.org")'
R --slave -e 'BiocManager::install(c("caret", "data.table", "dplyr", "furrr", "future", "future.apply",
"ggplot2", "grDevices", "gridExtra", "igraph", "irlba", "janitor", "methods", 
"plotly", "plyr", "purrr", "randomForest", "rhdf5", "rtracklayer", "scales", "stats", "stringr", 
"tibble", "tidyr", "umap", "utils"))'
R --slave -e 'install.packages("Signac",repos="http://cran.us.r-project.org")'

R --slave -e 'devtools::install_github("JinmiaoChenLab/Rphenograph")'
R --slave -e 'devtools::install_github("KrishnaswamyLab/MAGIC/Rmagic")'
R --slave -e 'devtools::install_github("lrylaarsdam/amethyst")'
R --slave -e 'install.packages("pheatmap",repos="http://cran.us.r-project.org")'
R --slave -e 'install.packages("plyr",repos="http://cran.us.r-project.org")'

#install chromvar for motif enrichment in DMRs
conda install conda-forge::gsl
R --slave -e 'BiocManager::install(c("DirichletMultinomial"))'
R --slave -e 'BiocManager::install(c("GO.db","motifmatchr","JASPAR2020","rGREAT"))'
R --slave -e 'devtools::install_github("GreenleafLab/chromVARmotifs")'
R --slave -e 'BiocManager::install("chromVAR")'

#predownload the hg38 ref data
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.annotation.gtf.gz

#Add premethyst dependencies
pip install h5py numpy pandas argparse
```

```bash
sudo singularity build amethyst.sif amethyst2/ 
```

# then build with proper env variables preloaded

amethyst.def
```bash
Bootstrap: localimage
From: /home/ubuntu/amethyst.sif

%environment
# set up all essential environment variables
export LC_ALL=C
export PATH=/opt/miniconda3/bin:$PATH
export PYTHONPATH=/opt/miniconda3/lib/python3.9/:$PYTHONPATH
export PYTHONPATH=/opt/miniconda3/lib/python3.12/:$PYTHONPATH
source /opt/etc/bashrc

%files
  gencode.v43.annotation.gtf.gz /container_ref/gencode.v43.annotation.gtf.gz
  premethyst_bed2h5.py /container_src/premethyst_bed2h5.py

%labels
    Author Ryan Mulqueen
    Version v0.2
    MyLabel Amethyst Container
```

```bash
sudo singularity build amethyst.sif amethyst.def
```

# Copykit image

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
sudo singularity build copykit.sif copykit.def

```
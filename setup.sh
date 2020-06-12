sudo apt install -y apt-transport-https software-properties-common
sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
sudo add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/'
sudo apt update
sudo apt install -y r-base
sudo apt install build-essential
sudo chmod -R 777 /usr/local/lib/R/site-library

sudo apt-get -y install aptitude libcurl4-openssl-dev libxml2-dev


R -e "install.packages('BiocManager'); BiocManager::install('SummarizedExperiment'); BiocManager::install('SingleCellExperiment'); BiocManager::install('MAST')"

conda create -n seqwell python=3.7.3 jupyterlab==2.0.1 numpy==1.18.1 pandas==1.0.3 matplotlib==3.2.1 scikit-learn==0.22.2.post1 scipy==1.4.1 palettable==3.3.0 statsmodels==0.11.1 fastcluster==1.1.26 scanpy==1.4.4.post1 leidenalg==0.7.0 louvain==0.6.1 python-igraph==0.8.0 umap-learn==0.3.10 joblib==0.14.1 seaborn==0.10.0 gcsfs==0.6.1 pyyaml==5.3.1 xlrd==1.2.0 mygene==3.1.0 dsub parallel

conda activate seqwell

pip install harmony-pytorch==0.1.2
pip install --user git+git://github.com/KrishnaswamyLab/MAGIC.git#subdirectory=python
pip install scrublet
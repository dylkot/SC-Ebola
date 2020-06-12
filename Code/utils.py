import os
import pandas as pd
import numpy as np
import scanpy as sc


###########################################################################################
### Functions for loading files from google cloud storage, and H5AD files in particular ###

from tempfile import NamedTemporaryFile
import subprocess

def save_adata(adata, filepath, ext='.h5ad', gcs=False):
    if gcs:
        temp = NamedTemporaryFile(suffix=ext, delete=False)
        temp.close()
        sc.write(temp.name, adata)
        subprocess.call('gsutil -m cp %s %s' % (temp.name, filepath), shell=True)
        subprocess.call('rm %s' % temp.name, shell=True)
    else:
        sc.write(filepath, adata)
    
    
    
def read_adata(filepath, ext='.h5ad', gcs=False):
    if gcs:
        temp = NamedTemporaryFile(suffix=ext, delete=False)
        temp.close()
        subprocess.call('gsutil -m cp %s %s' % (filepath, temp.name), shell=True)
        adata = sc.read(temp.name)
        subprocess.call('rm %s' % temp.name, shell=True)
    else:
        adata = sc.read(filepath)
    return(adata)

###########################################################################################
### Functions for running an individual subclustering iteration on scRNA-Seq data       ###

from harmony import harmonize

def subcluster_iteration(adata_in, min_cells=10, nhvgs=2000, npcs=20, n_neighbors=50,
                        min_dist=1.0, spread=2.0,resolution=1.,
                        umap_genestoplot=['CD14'], pc_genestoplot=['CD14'],
                        other_plot=['DPIc', 'louvain'], random_state=14,
                        harmony=False, harmony_key='frz_status',
                        regress_out_keys=None, n_jobs_regress=1,
                        harmony_theta=2, scale=True):
    ''' Assumes input data is already log TP10K normalized'''

    _adata = adata_in.copy()
    sc.pp.filter_genes(_adata, min_cells=min_cells)
    sc.pp.highly_variable_genes(_adata, n_top_genes=nhvgs)
    _adata = _adata[:, _adata.var['highly_variable']]
    
    if regress_out_keys is not None:
        _adata = _adata.copy()
        sc.pp.regress_out(_adata, regress_out_keys, n_jobs=n_jobs_regress, copy=False)
    
    if scale:
        sc.pp.scale(_adata, max_value=10)

    sc.tl.pca(_adata, svd_solver='arpack', random_state=14)

    sc.pl.pca(_adata, components=['1,2', '3,4', '5,6', '7,8'], color=pc_genestoplot,
          ncols=4, use_raw=True)

    sc.pl.pca_loadings(_adata, components=[1,2,3,4,5])
    sc.pl.pca_variance_ratio(_adata, log=True)

    if harmony:
        Z = harmonize(_adata.obsm['X_pca'], _adata.obs, batch_key = harmony_key,
                      random_state=random_state, theta=harmony_theta)
        _adata.obsm['X_harmony'] = Z
        sc.pp.neighbors(_adata, n_neighbors=n_neighbors, n_pcs=npcs, random_state=random_state,
                        use_rep='X_harmony')
    else:
        sc.pp.neighbors(_adata, n_neighbors=n_neighbors, n_pcs=npcs, random_state=random_state)
        
    sc.tl.umap(_adata, min_dist=min_dist, spread=spread, random_state=random_state)
    
    np.random.seed(random_state)
    sc.tl.leiden(_adata, resolution=resolution, random_state=random_state)
    
    fig = sc.pl.umap(_adata, color=umap_genestoplot, use_raw=True)
    fig = sc.pl.umap(_adata, color=other_plot)

    sc.tl.rank_genes_groups(_adata, 'leiden', method='wilcoxon')
    display(pd.DataFrame(_adata.uns['rank_genes_groups']['names']).head(20))
    return(_adata)



###########################################################################################
### Helper function for CyTOF                                                           ###

from sklearn.decomposition import PCA

def pca_cytof(X, random_state=14):
    PCAmod = PCA(n_components=20, random_state=random_state)
    pcs = PCAmod.fit_transform(X.X)
    X.obsm['X_pca'] = pcs
    X.varm['PCs'] = PCAmod.components_.T
    X.uns['pca'] = {} 
    X.uns['pca']['variance'] = PCAmod.explained_variance_
    X.uns['pca']['variance_ratio'] = PCAmod.explained_variance_ratio_
    

###########################################################################################
### Draw boxes on a scatter plot and capture indeces for the corresponding points       ###

def draw_box(xvars, yvars, ax, linestyle='--', color='k',
             linewidth=.5, skipline=[0,0,0,0]):
    cords = [0, 1, 2, 3, 0]
    for i in range(4):
        p1 = cords[i]
        p2 = cords[i+1]
        if skipline[i]!=1:
            ax.plot([xvars[p1], xvars[p2]], [yvars[p1], yvars[p2]],
                linestyle=linestyle, color=color, linewidth=linewidth)
        
        
def get_box_ind(xbound, ybound, xvals, yvals):
    cords = [0, 1, 2, 3, 0]
    ind = np.array([True]*len(xvals))
    for i in range(4):
        p1 = cords[i]
        p2 = cords[i+1]
        if (xbound[p2] - xbound[p1]) != 0:
            m = (ybound[p2] - ybound[p1]) / (xbound[p2] - xbound[p1])
            if i in [0, 3]:
                ind = ind & (yvals > (m*xvals-m*xbound[p1]+ybound[p1]))
            else:
                ind = ind & (yvals < (m*xvals-m*xbound[p1]+ybound[p1]))
        else:
            if i in [0, 3]:
                ind = ind & (xvals > xbound[p2])
            else:
                ind = ind & (xvals < xbound[p2])
            
    return(ind)
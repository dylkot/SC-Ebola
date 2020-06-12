import numpy as np
import pandas as pd
import os, errno
import datetime
import uuid
import itertools
import yaml
import subprocess
import scipy.sparse as sp


from scipy.spatial.distance import squareform
from sklearn.decomposition.nmf import non_negative_factorization
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from sklearn.utils import sparsefuncs


from fastcluster import linkage
from scipy.cluster.hierarchy import leaves_list

import matplotlib.pyplot as plt

import scanpy as sc

from cnmf import *

class cNMF_nonorm(cNMF):

    def get_norm_counts(self, counts, tpm,
                         high_variance_genes_filter = None,
                         num_highvar_genes = None
                         ):
        """
        Parameters
        ----------

        counts : anndata.AnnData
            Scanpy AnnData object (cells x genes) containing raw counts. Filtered such that
            no genes or cells with 0 counts
        
        tpm : anndata.AnnData
            Scanpy AnnData object (cells x genes) containing tpm normalized data matching
            counts

        high_variance_genes_filter : np.array, optional (default=None)
            A pre-specified list of genes considered to be high-variance.
            Only these genes will be used during factorization of the counts matrix.
            Must match the .var index of counts and tpm.
            If set to None, high-variance genes will be automatically computed, using the
            parameters below.

        num_highvar_genes : int, optional (default=None)
            Instead of providing an array of high-variance genes, identify this many most overdispersed genes
            for filtering

        Returns
        -------

        normcounts : anndata.AnnData, shape (cells, num_highvar_genes)
            A counts matrix containing only the high variance genes. No normalization.

        """

        if high_variance_genes_filter is None:
            ## Get list of high-var genes if one wasn't provided
            if sp.issparse(tpm.X):
                (gene_counts_stats, gene_fano_params) = get_highvar_genes_sparse(tpm.X, numgenes=num_highvar_genes)  
            else:
                (gene_counts_stats, gene_fano_params) = get_highvar_genes(np.array(tpm.X), numgenes=num_highvar_genes)
                
            high_variance_genes_filter = list(tpm.var.index[gene_counts_stats.high_var.values])
                
        ## Subset out high-variance genes
        norm_counts = counts[:, high_variance_genes_filter]                 
                    
        ## Save a \n-delimited list of the high-variance genes used for factorization
        open(self.paths['nmf_genes_list'], 'w').write('\n'.join(high_variance_genes_filter))

        ## Check for any cells that have 0 counts of the overdispersed genes
        zerocells = norm_counts.X.sum(axis=1)==0
        if zerocells.sum()>0:
            examples = norm_counts.obs.index[zerocells]
            print('Warning: %d cells have zero counts of overdispersed genes. E.g. %s' % (zerocells.sum(), examples[0]))
            print('Consensus step may not run when this is the case')
        
        return(norm_counts)
    
    
    
    def consensus(self, k, density_threshold_str='0.5', local_neighborhood_size = 0.30,show_clustering = False,
                  skip_density_and_return_after_stats = False, close_clustergram_fig=True):
        merged_spectra = load_df_from_npz(self.paths['merged_spectra']%k)
        norm_counts = sc.read(self.paths['normalized_counts'])

        if skip_density_and_return_after_stats:
            density_threshold_str = '2'
        density_threshold_repl = density_threshold_str.replace('.', '_')
        density_threshold = float(density_threshold_str)
        n_neighbors = int(local_neighborhood_size * merged_spectra.shape[0]/k)

        # Rescale topics such to length of 1.
        l2_spectra = (merged_spectra.T/np.sqrt((merged_spectra**2).sum(axis=1))).T


        if not skip_density_and_return_after_stats:
            # Compute the local density matrix (if not previously cached)
            topics_dist = None
            if os.path.isfile(self.paths['local_density_cache'] % k):
                local_density = load_df_from_npz(self.paths['local_density_cache'] % k)
            else:
                #   first find the full distance matrix
                topics_dist = squareform(fast_euclidean(l2_spectra.values))
                #   partition based on the first n neighbors
                partitioning_order  = np.argpartition(topics_dist, n_neighbors+1)[:, :n_neighbors+1]
                #   find the mean over those n_neighbors (excluding self, which has a distance of 0)
                distance_to_nearest_neighbors = topics_dist[np.arange(topics_dist.shape[0])[:, None], partitioning_order]
                local_density = pd.DataFrame(distance_to_nearest_neighbors.sum(1)/(n_neighbors),
                                             columns=['local_density'],
                                             index=l2_spectra.index)
                save_df_to_npz(local_density, self.paths['local_density_cache'] % k)
                del(partitioning_order)
                del(distance_to_nearest_neighbors)

            density_filter = local_density.iloc[:, 0] < density_threshold
            l2_spectra = l2_spectra.loc[density_filter, :]

        kmeans_model = KMeans(n_clusters=k, n_init=10, random_state=1)
        kmeans_model.fit(l2_spectra)
        kmeans_cluster_labels = pd.Series(kmeans_model.labels_+1, index=l2_spectra.index)

        # Find median usage for each gene across cluster
        median_spectra = l2_spectra.groupby(kmeans_cluster_labels).median()

        # Normalize median spectra to probability distributions.
        median_spectra = (median_spectra.T/median_spectra.sum(1)).T

        # Compute the silhouette score
        stability = silhouette_score(l2_spectra.values, kmeans_cluster_labels, metric='euclidean')

        # Obtain the reconstructed count matrix by re-fitting the usage matrix and computing the dot product: usage.dot(spectra)
        refit_nmf_kwargs = yaml.load(open(self.paths['nmf_run_parameters']), Loader=yaml.FullLoader)
        refit_nmf_kwargs.update(dict(
                                    n_components = k,
                                    H = median_spectra.values,
                                    update_H = False
                                    ))
        
        _, rf_usages = self._nmf(norm_counts.X,
                                          nmf_kwargs=refit_nmf_kwargs)
        rf_usages = pd.DataFrame(rf_usages, index=norm_counts.obs.index, columns=median_spectra.index)
        rf_pred_norm_counts = rf_usages.dot(median_spectra)

        # Compute prediction error as a frobenius norm
        if sp.issparse(norm_counts.X):
            prediction_error = ((norm_counts.X.todense() - rf_pred_norm_counts)**2).sum().sum()
        else:
            prediction_error = ((norm_counts.X - rf_pred_norm_counts)**2).sum().sum()
        
        consensus_stats = pd.DataFrame([k, density_threshold, stability, prediction_error],
                    index = ['k', 'local_density_threshold', 'stability', 'prediction_error'],
                    columns = ['stats'])

        if skip_density_and_return_after_stats:
            return consensus_stats
        
        save_df_to_npz(median_spectra, self.paths['consensus_spectra']%(k, density_threshold_repl))
        save_df_to_npz(rf_usages, self.paths['consensus_usages']%(k, density_threshold_repl))
        save_df_to_npz(consensus_stats, self.paths['consensus_stats']%(k, density_threshold_repl))
        save_df_to_text(median_spectra, self.paths['consensus_spectra__txt']%(k, density_threshold_repl))
        save_df_to_text(rf_usages, self.paths['consensus_usages__txt']%(k, density_threshold_repl))

        # Compute gene-scores for each GEP by regressing usage on Z-scores of TPM
        tpm = sc.read(self.paths['tpm'])
        tpm_stats = load_df_from_npz(self.paths['tpm_stats'])
        
        if sp.issparse(tpm.X):
            norm_tpm = (np.array(tpm.X.todense()) - tpm_stats['__mean'].values) / tpm_stats['__std'].values
        else:
            norm_tpm = (tpm.X - tpm_stats['__mean'].values) / tpm_stats['__std'].values
            
        usage_coef = fast_ols_all_cols(rf_usages.values, norm_tpm)
        usage_coef = pd.DataFrame(usage_coef, index=rf_usages.columns, columns=tpm.var.index)

        save_df_to_npz(usage_coef, self.paths['gene_spectra_score']%(k, density_threshold_repl))
        save_df_to_text(usage_coef, self.paths['gene_spectra_score__txt']%(k, density_threshold_repl))

        # Convert spectra to TPM units, and obtain results for all genes by running last step of NMF
        # with usages fixed and TPM as the input matrix
        norm_usages = rf_usages.div(rf_usages.sum(axis=1), axis=0)
        refit_nmf_kwargs.update(dict(
                                    H = norm_usages.T.values,
                                ))
        
        _, spectra_tpm = self._nmf(tpm.X.T, nmf_kwargs=refit_nmf_kwargs)
        spectra_tpm = pd.DataFrame(spectra_tpm.T, index=rf_usages.columns, columns=tpm.var.index)
        save_df_to_npz(spectra_tpm, self.paths['gene_spectra_tpm']%(k, density_threshold_repl))
        save_df_to_text(spectra_tpm, self.paths['gene_spectra_tpm__txt']%(k, density_threshold_repl))

        if show_clustering:
            if topics_dist is None:
                topics_dist = squareform(fast_euclidean(l2_spectra.values))
                # (l2_spectra was already filtered using the density filter)
            else:
                # (but the previously computed topics_dist was not!)
                topics_dist = topics_dist[density_filter.values, :][:, density_filter.values]


            spectra_order = []
            for cl in sorted(set(kmeans_cluster_labels)):

                cl_filter = kmeans_cluster_labels==cl

                if cl_filter.sum() > 1:
                    cl_dist = squareform(topics_dist[cl_filter, :][:, cl_filter])
                    cl_dist[cl_dist < 0] = 0 #Rarely get floating point arithmetic issues
                    cl_link = linkage(cl_dist, 'average')
                    cl_leaves_order = leaves_list(cl_link)

                    spectra_order += list(np.where(cl_filter)[0][cl_leaves_order])
                else:
                    ## Corner case where a component only has one element
                    spectra_order += list(np.where(cl_filter)[0])


            from matplotlib import gridspec
            import matplotlib.pyplot as plt

            width_ratios = [0.5, 9, 0.5, 4, 1]
            height_ratios = [0.5, 9]
            fig = plt.figure(figsize=(sum(width_ratios), sum(height_ratios)))
            gs = gridspec.GridSpec(len(height_ratios), len(width_ratios), fig,
                                    0.01, 0.01, 0.98, 0.98,
                                   height_ratios=height_ratios,
                                   width_ratios=width_ratios,
                                   wspace=0, hspace=0)

            dist_ax = fig.add_subplot(gs[1,1], xscale='linear', yscale='linear',
                                      xticks=[], yticks=[],xlabel='', ylabel='',
                                      frameon=True)

            D = topics_dist[spectra_order, :][:, spectra_order]
            dist_im = dist_ax.imshow(D, interpolation='none', cmap='viridis', aspect='auto',
                                rasterized=True)

            left_ax = fig.add_subplot(gs[1,0], xscale='linear', yscale='linear', xticks=[], yticks=[],
                xlabel='', ylabel='', frameon=True)
            left_ax.imshow(kmeans_cluster_labels.values[spectra_order].reshape(-1, 1),
                            interpolation='none', cmap='Spectral', aspect='auto',
                            rasterized=True)


            top_ax = fig.add_subplot(gs[0,1], xscale='linear', yscale='linear', xticks=[], yticks=[],
                xlabel='', ylabel='', frameon=True)
            top_ax.imshow(kmeans_cluster_labels.values[spectra_order].reshape(1, -1),
                              interpolation='none', cmap='Spectral', aspect='auto',
                                rasterized=True)


            hist_gs = gridspec.GridSpecFromSubplotSpec(3, 1, subplot_spec=gs[1, 3],
                                   wspace=0, hspace=0)

            hist_ax = fig.add_subplot(hist_gs[0,0], xscale='linear', yscale='linear',
                xlabel='', ylabel='', frameon=True, title='Local density histogram')
            hist_ax.hist(local_density.values, bins=np.linspace(0, 1, 50))
            hist_ax.yaxis.tick_right()

            xlim = hist_ax.get_xlim()
            ylim = hist_ax.get_ylim()
            if density_threshold < xlim[1]:
                hist_ax.axvline(density_threshold, linestyle='--', color='k')
                hist_ax.text(density_threshold  + 0.02, ylim[1] * 0.95, 'filtering\nthreshold\n\n', va='top')
            hist_ax.set_xlim(xlim)
            hist_ax.set_xlabel('Mean distance to k nearest neighbors\n\n%d/%d (%.0f%%) spectra above threshold\nwere removed prior to clustering'%(sum(~density_filter), len(density_filter), 100*(~density_filter).mean()))

            fig.savefig(self.paths['clustering_plot']%(k, density_threshold_repl), dpi=250)
            if close_clustergram_fig:
                plt.close(fig)
                
                


if __name__=="__main__":
    """
    Example commands for now:

        output_dir="/Users/averes/Projects/Melton/Notebooks/2018/07-2018/cnmf_test/"


        python cnmf.py prepare --output-dir $output_dir \
           --name test --counts /Users/averes/Projects/Melton/Notebooks/2018/07-2018/cnmf_test/test_data.df.npz \
           -k 6 7 8 9 --n-iter 5

        python cnmf.py factorize  --name test --output-dir $output_dir

        THis can be parallelized as such:

        python cnmf.py factorize  --name test --output-dir $output_dir --total-workers 2 --worker-index WORKER_INDEX (where worker_index starts with 0)

        python cnmf.py combine  --name test --output-dir $output_dir

        python cnmf.py consensus  --name test --output-dir $output_dir

    """

    import sys, argparse
    parser = argparse.ArgumentParser()

    parser.add_argument('command', type=str, choices=['prepare', 'factorize', 'combine', 'consensus', 'k_selection_plot'])
    parser.add_argument('--name', type=str, help='[all] Name for analysis. All output will be placed in [output-dir]/[name]/...', nargs='?', default='cNMF')
    parser.add_argument('--output-dir', type=str, help='[all] Output directory. All output will be placed in [output-dir]/[name]/...', nargs='?', default='.')

    parser.add_argument('-c', '--counts', type=str, help='[prepare] Input (cell x gene) counts matrix as df.npz or tab delimited text file')
    parser.add_argument('-k', '--components', type=int, help='[prepare] Numper of components (k) for matrix factorization. Several can be specified with "-k 8 9 10"', nargs='+')
    parser.add_argument('-n', '--n-iter', type=int, help='[prepare] Numper of factorization replicates', default=100)
    parser.add_argument('--total-workers', type=int, help='[all] Total number of workers to distribute jobs to', default=1)
    parser.add_argument('--seed', type=int, help='[prepare] Seed for pseudorandom number generation', default=None)
    parser.add_argument('--genes-file', type=str, help='[prepare] File containing a list of genes to include, one gene per line. Must match column labels of counts matrix.', default=None)
    parser.add_argument('--numgenes', type=int, help='[prepare] Number of high variance genes to use for matrix factorization.', default=2000)
    parser.add_argument('--tpm', type=str, help='[prepare] Pre-computed (cell x gene) TPM values as df.npz or tab separated txt file. If not provided TPM will be calculated automatically', default=None)
    parser.add_argument('--beta-loss', type=str, choices=['frobenius', 'kullback-leibler', 'itakura-saito'], help='[prepare] Loss function for NMF.', default='frobenius')
    parser.add_argument('--densify', dest='densify', help='[prepare] Treat the input data as non-sparse', action='store_true', default=False)

    
    parser.add_argument('--worker-index', type=int, help='[factorize] Index of current worker (the first worker should have index 0)', default=0)
    
    parser.add_argument('--local-density-threshold', type=str, help='[consensus] Threshold for the local density filtering. This string must convert to a float >0 and <=2', default='0.5')
    parser.add_argument('--local-neighborhood-size', type=float, help='[consensus] Fraction of the number of replicates to use as nearest neighbors for local density filtering', default=0.30)
    parser.add_argument('--show-clustering', dest='show_clustering', help='[consensus] Produce a clustergram figure summarizing the spectra clustering', action='store_true')

    args = parser.parse_args()

    cnmf_obj = cNMF(output_dir=args.output_dir, name=args.name)
    cnmf_obj._initialize_dirs()
    
    if args.command == 'prepare':

        if args.counts.endswith('.h5ad'):
            input_counts = sc.read(args.counts)
        else:
            ## Load txt or compressed dataframe and convert to scanpy object
            if args.counts.endswith('.npz'):
                input_counts = load_df_from_npz(args.counts)
            else:
                input_counts = pd.read_csv(args.counts, sep='\t', index_col=0)
                
            if args.densify:
                input_counts = sc.AnnData(X=input_counts.values,
                                       obs=pd.DataFrame(index=input_counts.index),
                                       var=pd.DataFrame(index=input_counts.columns))
            else:
                input_counts = sc.AnnData(X=sp.csr_matrix(input_counts.values),
                                       obs=pd.DataFrame(index=input_counts.index),
                                       var=pd.DataFrame(index=input_counts.columns))

                
        if sp.issparse(input_counts.X) & args.densify:
            input_counts.X = np.array(input_counts.X.todense())
 
        if args.tpm is None:
            tpm = compute_tpm(input_counts)
            sc.write(cnmf_obj.paths['tpm'], tpm)
        elif args.tpm.endswith('.h5ad'):
            subprocess.call('cp %s %s' % (args.tpm, cnmf_obj.paths['tpm']), shell=True)
            tpm = sc.read(cnmf_obj.paths['tpm'])
        else:
            if args.tpm.endswith('.npz'):
                tpm = load_df_from_npz(args.tpm)
            else:
                tpm = pd.read_csv(args.tpm, sep='\t', index_col=0)
            
            if args.densify:
                tpm = sc.AnnData(X=tpm.values,
                            obs=pd.DataFrame(index=tpm.index),
                            var=pd.DataFrame(index=tpm.columns)) 
            else:
                tpm = sc.AnnData(X=sp.csr_matrix(tpm.values),
                            obs=pd.DataFrame(index=tpm.index),
                            var=pd.DataFrame(index=tpm.columns)) 

            sc.write(cnmf_obj.paths['tpm'], tpm)
        
        if sp.issparse(tpm.X):
            gene_tpm_mean = np.array(tpm.X.mean(axis=0)).reshape(-1)
            gene_tpm_stddev = var_sparse_matrix(tpm.X)**.5
        else:
            gene_tpm_mean = np.array(tpm.X.mean(axis=0)).reshape(-1)
            gene_tpm_stddev = np.array(tpm.X.std(axis=0, ddof=0)).reshape(-1)
            
            
        input_tpm_stats = pd.DataFrame([gene_tpm_mean, gene_tpm_stddev],
             index = ['__mean', '__std']).T
        save_df_to_npz(input_tpm_stats, cnmf_obj.paths['tpm_stats'])
        
        if args.genes_file is not None:
            highvargenes = open(args.genes_file).read().rstrip().split('\n')
        else:
            highvargenes = None

        norm_counts = cnmf_obj.get_norm_counts(input_counts, tpm, num_highvar_genes=args.numgenes,
                                               high_variance_genes_filter=highvargenes)
        cnmf_obj.save_norm_counts(norm_counts)
        (replicate_params, run_params) = cnmf_obj.get_nmf_iter_params(ks=args.components, n_iter=args.n_iter, random_state_seed=args.seed, beta_loss=args.beta_loss)
        cnmf_obj.save_nmf_iter_params(replicate_params, run_params)


    elif args.command == 'factorize':
        cnmf_obj.run_nmf(worker_i=args.worker_index, total_workers=args.total_workers)

    elif args.command == 'combine':
        run_params = load_df_from_npz(cnmf_obj.paths['nmf_replicate_parameters'])

        if type(args.components) is int:
            ks = [args.components]
        elif args.components is None:
            ks = sorted(set(run_params.n_components))
        else:
            ks = args.components

        for k in ks:
            cnmf_obj.combine_nmf(k)

    elif args.command == 'consensus':
        run_params = load_df_from_npz(cnmf_obj.paths['nmf_replicate_parameters'])

        if type(args.components) is int:
            ks = [args.components]
        elif args.components is None:
            ks = sorted(set(run_params.n_components))
        else:
            ks = args.components

        for k in ks:
            merged_spectra = load_df_from_npz(cnmf_obj.paths['merged_spectra']%k)
            cnmf_obj.consensus(k, args.local_density_threshold, args.local_neighborhood_size, args.show_clustering)

    elif args.command == 'k_selection_plot':
        cnmf_obj.k_selection_plot()

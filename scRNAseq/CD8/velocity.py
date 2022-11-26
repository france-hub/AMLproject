import os
import scvelo as scv
import pandas as pd
os.getcwd()
os.chdir("/home/francesco/Documents/AML_project/scRNA_AMLproj/python_scRNAseq/")
adata = scv.read("CD8_veloscv.h5ad")
pal_ident = ["#F0E442","#E69F00", "#009E73", "#56B4E9" ]
adata.uns['clusters']= pal_ident
scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
#Plot velocity onto UMAP with clusters
scv.tl.recover_dynamics(adata)
scv.tl.velocity(adata, mode = 'dynamical')
#save results
adata.__dict__['_raw'].__dict__['_var'] = adata.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'})
#adata.write('adata.h5ad')
adata = scv.read('adata.h5ad')
scv.tl.velocity_graph(adata)
scv.pl.velocity_embedding_stream(adata, basis="umap", title = '', color="clusters2", palette=pal_ident)#, save = "velocity_dyn")
#Define Naive as root
df = pd.DataFrame(index=adata.obs_names).reset_index()
adata.uns['root_key'] = df.index[df['index'].isin(adata.obs_names[adata.obs['clusters'] == "Naive"])][0]
scv.tl.latent_time(adata, root_key= "root_key")
scv.pl.scatter(adata, color='latent_time', cmap='viridis')
#latent time
scv.pl.scatter(adata, color='velocity_pseudotime', title = '', cmap='viridis', alpha = 0.4, save = "latent_time.png")

#Run Paga clusters
scv.tl.paga(adata, use_time_prior= "latent_time", groups='clusters')
scv.pl.paga(adata,basis = "umap", title = "", threshold = 0.05,arrowsize = 10,edge_width_scale = 0.5, palette = pal_ident, save = "paga.png")

#Run Paga clusters2
pal_ident2 = ["#F0E442", "#E69F00", "#009E73", "#0072B2", "#56B4E9"]
adata.uns['clusters2']= pal_ident2
scv.tl.paga(adata, use_time_prior= "latent_time", groups='clusters2')
scv.pl.paga(adata,basis = "umap", title = "", threshold = 0.05,arrowsize = 10,edge_width_scale = 0.5, palette = pal_ident2, save = "paga.png")
scv.set_figure_params(figsize="5, 5")
scv.pl.paga_compare(adata, threshold = 0.1, title = '', arrowsize = 10, edge_width_scale = 0.5, transitions = "transitions_confidence", dashed_edges = "connectivities", save = "paga.compare.png")


top_genes = adata.var['fit_likelihood'].sort_values(ascending=False).index[:100]
genes = ["CCR7","LEF1","BACH2","SELL","IL7R","CD27","GZMK","PDCD1","CX3CR1","GZMA","GZMB","PRF1","FGFBP2","ADGRG1","KLRD1","NKG7","GNLY","ZNF683"]
scv.pl.heatmap(adata, var_names=genes, sort = False,  sortby='latent_time', color_map = "coolwarm", col_color='clusters',  yticklabels=True, figsize=(4, 4), save='.png')
l_time = (

    adata.obs.groupby(["clusters2"])
    .apply(
        lambda x: x["latent_time"]
        .value_counts(dropna=False)
        .rename_axis("latent_time")
    )
)
l_time.to_csv("l_time.csv")

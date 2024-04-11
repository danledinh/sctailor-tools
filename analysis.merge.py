import scanpy as sc
import anndata as ad
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import ShuffleSplit

def ingest_counts(path, name_conversion_dict):
    """
    Ingest Isoquant cell x gene count matrix.
    path = str, path to TSV
    name_conversion_dict = dict, covert between TSV feature names and desired naming system
    """
    mat = pd.read_csv(value, index_col=0, sep="\t")
    mat.index = mat.index.map(name_conversion_dict)
    mat = mat.reset_index().dropna().groupby("#feature_id").sum()
    mat = mat.astype(int)
    return mat

# transform to cell-norm + log1p
imm_mat = ingest_counts(imm_path).loc[common_imm_tx,:].loc[:,common_cb]
imm_ad = make_adata(imm_mat)
sc.pp.filter_cells(imm_ad, min_counts=10)
sc.pp.normalize_total(imm_ad, target_sum=1e6)
sc.pp.log1p(imm_ad)

can_mat = ingest_counts(can_path).loc[common_can_tx,:].loc[:,common_cb]
can_ad = make_adata(can_mat)
sc.pp.filter_cells(can_ad, min_counts=10)
sc.pp.normalize_total(can_ad, target_sum=1e6)
sc.pp.log1p(can_ad)

common_cb = set(imm_ad.obs_names) & set(can_ad.obs_names)
imm_mat = pd.DataFrame(imm_ad.X.toarray(),columns=imm_ad.var_names, index=imm_ad.obs_names).T.loc[:,common_cb]
can_mat = pd.DataFrame(can_ad.X.toarray(),columns=can_ad.var_names, index=can_ad.obs_names).T.loc[:,common_cb]

# compute lin reg + cross val
common_intersect = set(imm_mat.index) & set(can_mat.index) & set(intersect_l)
X = can_mat.loc[common_intersect,:].mean(1).values.reshape((-1,1))
y = imm_mat.loc[common_intersect,:].mean(1).values.reshape((-1,1))
ss = ShuffleSplit(n_splits=10, test_size=0.9, random_state=0)
slope_arr = np.empty(0)
for train_index, test_index in ss.split(X):

    reg = LinearRegression(fit_intercept=False).fit(X[train_index], y[train_index])
    slope_arr = np.append(slope_arr,reg.coef_)

mu_slope = np.mean(slope_arr)
reg = LinearRegression().fit(X, y)
r2 = reg.score(X,y)**2

# adjust values and merge
can_mat = can_mat * mu_slope
merge_mat = imm_mat.loc[set(common_imm_tx)-set(common_intersect),:].append(can_mat)

# predict doublets
dubs_ad = ad.AnnData(np.expm1(merge_mat.T).astype(int))
sc.external.pp.scrublet(dubs_ad)

# create adata
n_neighbors=50
n_pcs=50
merge_gene_ad = make_adata(merge_mat)
sc.pp.calculate_qc_metrics(merge_gene_ad, inplace=True)
merge_gene_ad.raw = merge_gene_ad
sc.pp.scale(merge_gene_ad)
sc.tl.pca(merge_gene_ad, svd_solver='arpack')
sc.pp.neighbors(merge_gene_ad, n_neighbors=n_neighbors, n_pcs=n_pcs)
sc.tl.umap(merge_gene_ad)

merge_gene_ad.obs["doublet"] = merge_gene_ad.obs.index.map(dubs_ad.obs["predicted_doublet"].to_dict())
merge_gene_ad.obs["doublet"] = merge_gene_ad.obs["doublet"].astype(int)

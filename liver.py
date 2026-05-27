pip install scanpy
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import scipy.sparse as sp
import gseapy as gp
import statsmodels.api as sm
from scipy.stats import mannwhitneyu
import seaborn as sns

# load dataset 
ad = sc.read_h5ad("liver_aging_subsampled.h5ad")

sc.pp.normalize_total(ad)
sc.pp.log1p(ad)
sc.pp.highly_variable_genes(ad, n_top_genes=2000)

# principal component analysis 
sc.tl.pca(ad)

sc.pp.neighbors(ad)
sc.tl.umap(ad)
sc.pl.umap(ad, color=["AgeBracket","cell_type","sex"])

print(ad.obs.columns)

print(ad.obs[["cell_type", "donor_age", "AgeGroup"]])

# counts total cells per cell type 
counts = ad.obs.groupby(["AgeGroup", "cell_type"]).size().reset_index(name="n_cells")
counts.head()

counts = ad.obs.groupby(["AgeGroup", "cell_type"]).size().reset_index(name="n_cells")

abundance = counts.pivot(
    index="cell_type",
    columns="AgeGroup",
    values="n_cells"
).fillna(0)

abundance.plot(kind="barh", figsize=(8,6))

plt.xlabel("Number of cells")
plt.ylabel("Cell type")
plt.title("Cell-type composition by AgeGroup")
plt.legend(title="AgeGroup")

plt.tight_layout()
plt.show()

# get hepatocyte subset for downstream analysis 
ad_hep = ad[ad.obs["cell_type"] == "hepatocyte", :].copy()
print(ad_hep.shape)

X = (ad_hep.layers["counts"]
     if "counts" in ad_hep.layers
     else (ad_hep.raw.X if ad_hep.raw is not None else ad_hep.X))

if sp.issparse(X):
    X = X.tocsr()

obs = ad_hep.obs.loc[:, ~ad_hep.obs.columns.duplicated()].copy()

# create a pseudo-bulk ID per donor
obs["pb_id"] = obs["donor_id"].astype(str) + " || " + obs["cell_type"].astype(str)

pb_ids, inv = np.unique(obs["pb_id"].values, return_inverse=True)
G = np.zeros((len(pb_ids), X.shape[1]), dtype=np.float64)

for i in range(len(pb_ids)):
    rows = np.where(inv == i)[0]
    s = X[rows].sum(axis=0)
    if sp.issparse(s):
        s = s.A1
    G[i, :] = s

counts = pd.DataFrame(G, index=pb_ids, columns=ad_hep.var_names).round().astype(int)

# sample metadata (per donor)
meta = obs.groupby("pb_id")[["donor_id","AgeGroup","sex","institute","assay"]].first()
meta = meta.dropna(subset=["AgeGroup"])
counts = counts.loc[meta.index]

# filter low-count genes
keep = (counts >= 10).sum(axis=0) >= 2
counts = counts.loc[:, keep]
print(counts.shape)
meta.head()

log_counts = np.log1p(counts)

# numeric covariate for AgeGroup
meta["AgeGroup_code"] = meta["AgeGroup"].map({"Ped": 0, "Adult": 1})

# differential expression 
results = []
for gene in log_counts.columns:
    y = log_counts[gene].values
    X = sm.add_constant(meta["AgeGroup_code"].values)
    model = sm.OLS(y, X).fit()
    coef = model.params[1]
    pval = model.pvalues[1]
    results.append({"gene": gene, "log2FC": coef, "pval": pval})

de = pd.DataFrame(results)
de["FDR"] = sm.stats.multipletests(de["pval"], method="fdr_bh")[1]
de = de.sort_values("pval")
de.head()

de.to_csv(("output.csv"))
print(de[5:])
print(de[:-5])

de.sort_values("log2FC")

sig_de = de[de["pval"] < 0.05]

# genes highly expressed in Adults
sig_de = sig_de.sort_values("log2FC", ascending=False) 

sig_de.head()

gene_list = [
    "PLA2G2A",  
    "XIST",      
    "NEAT1",     # ENSG00000115421
    "IGKC",      # ENSG00000242498
    "LGALS4" # ENSG00000114786
]

# gene set enrichment analysis using the GSEApy library
enr = gp.enrichr(
    gene_list=gene_list,
    gene_sets=['GO_Biological_Process_2023'],
    organism='Human',  
    outdir=None,        
)

# sort biological pathways by significance in adults
results = enr.results.sort_values("Adjusted P-value")
results[["Term", "Overlap", "Adjusted P-value", "Combined Score"]].head(10)

gp.barplot(
    enr.results,
    column="P-value",
    title="GO Biological Process Enrichment"
)

# sort by genes highly expressed in Children
sig_de = sig_de.sort_values("log2FC", ascending=True) 

sig_de.head()

gene_list = [
    "EIF3J-DT",  
    "LATS1",      
    "PAPOLG",     # ENSG00000115421
    "ARPIN",      # ENSG00000242498
    "ABHD14A-ACY1" # ENSG00000114786
]

# sort biological pathways by significance in children
enr = gp.enrichr(
    gene_list=gene_list,
    gene_sets=['GO_Biological_Process_2023'],
    organism='Human', 
    outdir=None,      
)

results = enr.results.sort_values("Adjusted P-value")
results[["Term", "Overlap", "Adjusted P-value", "Combined Score"]].head(10)

gp.barplot(
    enr.results,
    column="P-value",
    title="GO Biological Process Enrichment"
)

# run preranked GSEA to identify biological pathways enriched in Adult vs Pediatric hepatocytes
gene_map = ad_hep.var.copy()
gene_map["ensembl_id"] = gene_map.index
print(gene_map.columns)

# add gene symbols to DE table
de2 = de.merge(
    gene_map[["ensembl_id", "feature_name"]],
    left_on="gene",
    right_on="ensembl_id",
    how="left"
)

# prepare ranks using gene name, not ENSG IDs
ranks = de2[["feature_name", "log2FC"]].dropna()
ranks = ranks.rename(columns={"feature_name": "gene"})

# remove duplicated gene symbols
ranks = ranks.groupby("gene", as_index=False)["log2FC"].mean()

# sort strongest Adult-up genes first
ranks = ranks.sort_values("log2FC", ascending=False)

enr = gp.prerank(
    rnk=ranks,
    gene_sets=["Reactome_2022", "KEGG_2021_Human", "MSigDB_Hallmark_2020"],
    permutation_num=1000,
    outdir=None,
    seed=42,
    min_size=5,
    max_size=1000
)

enr.res2d.sort_values("FDR q-val").head(20)

# cell-type abundance analysis between peds and adults donors

obs = ad.obs.copy()  

# count cells per donor per cell_type
counts = (
    obs.groupby(["donor_id", "AgeGroup", "cell_type"])
    .size()
    .reset_index(name="n_cells")
)

# total cells per donor
totals = counts.groupby("donor_id")["n_cells"].sum().rename("total_cells")

# merge and compute proportions
counts = counts.merge(totals, on="donor_id")
counts["prop"] = counts["n_cells"] / counts["total_cells"]

display(counts.head())

results = []
for ct in counts["cell_type"].unique():
    df = counts[counts["cell_type"] == ct]
    ped = df.loc[df["AgeGroup"] == "Ped", "prop"]
    adult = df.loc[df["AgeGroup"] == "Adult", "prop"]
    if len(ped) >= 2 and len(adult) >= 2:
        stat, p = mannwhitneyu(ped, adult, alternative="two-sided")
        results.append({
            "cell_type": ct,
            "Ped_mean": ped.mean(),
            "Adult_mean": adult.mean(),
            "pvalue": p
        })

res_df = pd.DataFrame(results).sort_values("pvalue")
display(res_df)

top_cts = res_df.head(5)["cell_type"]

sns.boxplot(
    data=counts[counts["cell_type"].isin(top_cts)],
    x="AgeGroup",
    y="prop",
    hue="AgeGroup"
)

plt.title("Top differential cell-type abundances by AgeGroup")
plt.show()

keep = [
    "hepatocyte",
    "cholangiocyte",
    "Kupffer cell",
    "endothelial cell of hepatic sinusoid",
    "stellate cell"
]

filtered = res_df[res_df["cell_type"].isin(keep)]
print(filtered)


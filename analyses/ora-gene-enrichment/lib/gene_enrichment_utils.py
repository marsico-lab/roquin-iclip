import pandas as pd

import math
import ast
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.axes import Axes

# set matplotlib font to arial
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 10

import scipy.stats
from statsmodels.stats.multitest import multipletests

## IO functions

def read_txt_file(file_path):
    """
    Reads a text file and returns its contents as a list of stripped lines.

    Args:
        file_path (str): The path to the text file to be read.

    Returns:
        list: A list of strings, each representing a line from the file with leading and trailing whitespace removed.
    """
    with open(file_path, "r") as file:
        lines = file.readlines()
    return [line.strip() for line in lines]

## Data wrangling functions

def convert_to_array(s):
    # Use ast.literal_eval to safely evaluate the string as a Python literal (list)
    return np.array(ast.literal_eval(s))

def convert_ensembl_id(ensembl_ids, ensembl_map):
    """
    Convert a list of Ensembl gene IDs to their corresponding gene symbols or entrez ids, depending on the dict passed.

    This function strips the version number from Ensembl gene IDs (i.e., the substring after the dot)
    and then maps the stripped Ensembl IDs to their corresponding gene symbols or entrez IDs using the provided mapping.

    Args:
        ensembl_ids (list of str): A list of Ensembl gene IDs.
        ensembl_symbol_map (dict): A dictionary mapping Ensembl gene IDs (without version numbers) to gene symbols.

    Returns:
        list of str: A list of gene symbols corresponding to the provided Ensembl gene IDs.
    """

    # Strip version number from ensembl gene ids i.e. the substring after the dot
    ensembl_ids = [gene.split(".")[0] for gene in ensembl_ids]
    # Do conversion
    ensembl_symbols = [ensembl_map[ens_id] for ens_id in ensembl_ids if ens_id in ensembl_map]
    return ensembl_symbols

## Pathway enrichment analysis (ORA)

def perform_ora(study_genes, gene_sets, background_genes):
    """
    Performs Overrepresentation Analysis (ORA) using a hypergeometric test (Fisher's exact test).
    
    Parameters:
    - study_genes: a collection (list or set) of genes of interest.
    - gene_sets: a dictionary where keys are term names (e.g. pathway names)
                 and values are collections (list or set) of genes associated with that term.
    - background_genes: a collection (list or set) representing the background gene set.
    
    Returns:
    A list of dictionaries. Each dictionary contains:
        - term: the term name
        - k: number of study genes in the term
        - K: number of background genes in the term
        - n: number of study genes in the background
        - N: total number of background genes
        - p_value: the enrichment p-value from the hypergeometric test
        - FDR: FDR-adjusted p-value (Benjamini–Hochberg correction)

    Example usage:

    # Define a background gene set
    background = {
        "gene1", "gene2", "gene3", "gene4", "gene5", 
        "gene6", "gene7", "gene8", "gene9", "gene10"
    }

    # Define study genes (e.g., from an experiment)
    study = {"gene1", "gene3", "gene5", "gene7"}

    # Define gene sets (e.g., pathways or terms)
    gene_sets = {
        "Pathway_A": {"gene1", "gene2", "gene3"},
        "Pathway_B": {"gene3", "gene4", "gene6"},
        "Pathway_C": {"gene7", "gene8", "gene9", "gene10"},
        "Pathway_D": {"gene2", "gene5", "gene8"}
    }

    # Perform ORA
    ora_results = perform_ora(study, gene_sets, background)

    """
    results = []
    
    # Convert inputs to sets (if not already) for efficiency
    study_genes = set(study_genes)
    background_genes = set(background_genes)
    
    # Only consider study genes that are in the background
    n = len(study_genes & background_genes)
    N = len(background_genes)
    
    for term, genes in gene_sets.items():
        term_genes = set(genes) & background_genes
        K = len(term_genes)
        k = len(study_genes & term_genes)
        
        # If there are no genes in this term in the background, skip or assign a p-value of 1.0
        if K == 0:
            p_value = 1.0
        else:
            # Calculate p-value using the survival function:
            # P(X >= k) = hypergeom.sf(k-1, N, K, n)
            p_value = scipy.stats.hypergeom.sf(k - 1, N, K, n)
            # We could use also Fisher's exact test, they are equivalent
            # contingency_table = [
            # [k, n - k],
            # [K - k, N - K - (n - k)]
            # ]
            # _, p_value = scipy.stats.fisher_exact(contingency_table, alternative='greater')
        # study_genes = study_genes & term_genes
        results.append({
            "term": term,
            "k": k,
            "K": K,
            "n": n,
            "N": N,
            "p_value": p_value,
            "study_genes": list(study_genes & term_genes)
        })
    
    results_df = pd.DataFrame(results)
    results_df['p_value'] = results_df['p_value'].astype(float)
    # Perform multiple testing correction with Benjamini-Hochberg
    reject, fdrs, _, _ = multipletests(results_df['p_value'], method='fdr_bh')
        
    results_df["FDR"] = fdrs
    results_df = results_df.sort_values(by='FDR')
    # Return only significant results
    results_df = results_df[results_df['FDR'] <= 0.05]

    # Make column names pretty
    # results_df['ratio_in_study'] = results_df[['k', 'n']].apply(lambda x: f"{str(x[0])}/{str(x[1])}", axis=1)
    # results_df['ratio_in_background'] = results_df[['K', 'N']].apply(lambda x: f"{str(x[0])}/{str(x[1])}", axis=1)
    results_df['study_count'] = results_df['k']
    results_df['background_count'] = results_df['K']

    results_df = results_df[['term', 'study_count', 'background_count', 'p_value', 'FDR', 'study_genes']]
    
    return results_df

## Plotting functions

# Modified from # https://github.com/harryhaller001/keggtools/blob/main/keggtools/analysis.py line 288
def plot_enrichment_scatter(
    result_df: pd.DataFrame = None,
    ax: Axes | None = None,
    x_column: str = 'study_count',
    y_column: str = 'name',
    pvalue_column: str = 'p_fdr_bh', 
    bg_count_column: str = 'background_count',
    sc_column: str = 'study_count',
    figsize: tuple[int, int] = (7, 7),
    cmap: str = "inferno",
    min_study_count: int = 1,
    max_pval: float | None = None,
    use_percent_study_count: bool = True,
    top_n=None, 
    outpath=None) -> Axes:
    """
    Plots a scatter plot for gene enrichment analysis results.

    Parameters:
    result_df (pd.DataFrame): DataFrame containing enrichment terms data.
    x_column (str): Column name for x-axis values (default is 'study_count').
    y_column (str): Column name for y-axis values (default is 'name').
    pvalue_column (str): Column name for p-values to be used for coloring the points (default is 'p_fdr_bh').
    figsize (tuple): Size of the figure (default is (7, 7)).
    top_n (int, optional): Number of top terms to plot based on x_column values. If None, plot all terms (default is None).
    outpath (str, optional): Path to save the plot. If None, the plot is not saved (default is None).

    Returns:
    None
    """

    data = result_df.copy()

    # Filter out pathways with no genes found
    data = data[data[x_column] >= min_study_count]

    if max_pval is not None:
        data = data[data[pvalue_column] <= max_pval]
    
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    # Sort by ratio of genes in study and select the 'top_n'
    data['ratio_in_study'] = data[x_column] / data[bg_count_column]
    data = data.sort_values('ratio_in_study', ascending=True)

    if top_n:
        data = data.head(top_n)

    # NOTE: needed only for long GO terms as y-labels
    # import textwrap
    # # For clarity, you might want to separate long GO term names into two lines
    # data['short_name'] = data[y_column].apply(lambda x: '\n'.join(textwrap.wrap(x, width=50)))

    x_values = data[x_column]
    x_label = "genes in study"

    if use_percent_study_count is True:
        x_values = data['ratio_in_study'] * 100
        x_label = "genes in study [%]"

    scatter = ax.scatter(
        x=x_values,
        y=data[y_column],
        c=[-math.log10(x) for x in data[pvalue_column]],
        cmap=cmap,
    )

    cbar = plt.colorbar(scatter)
    cbar.set_label('- log10(p value)')
    ax.set_xlabel(x_label)
    ax.grid(visible=None)

    plt.tight_layout()

    if outpath:
        plt.savefig(outpath, bbox_inches='tight')
        plt.savefig(outpath.replace(".pdf", ".png"), bbox_inches='tight')

    plt.show()
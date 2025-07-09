# figure-1
Figure 1: Sequence and structure determinants of Roquin-1 binding to the CD4+ T cell transcriptome.

- Panel b): Proportions of binding sites across protein coding biotypes, before and after iDKO filter, per batch

```
TODO
```

- Panel c): Correlation of cross-link counts at the final binding sites between the two batches
```
Rscript panel_c_scatter_batches_reproducibility.R \
    --bs-bed-path <path_to_bs_bed_file> \
    --batch1-bw-path <path_to_batch1_bw_file> \
    --batch2-bw-path <path_to_batch2_bw_file> \
    --out-filename <output_filename>
```

- Panel d): Number of binding sites and their biotype annotations with relative proportions

```
python panel_d_bs_4tx_pie_chart.py \
    --bs-bed-path <path_to_bs_bed_file> \
    --palette mouse \ # or 'human' for human palette
    --out-filename <output_filename>
```

- Panel f): Co-occurrence analysis with relative positions between structure and primary sequence motifs for four of the sequence motifs and the top scoring structural motifs

```
TODO
```
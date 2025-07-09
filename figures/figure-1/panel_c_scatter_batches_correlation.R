suppressMessages(library(SummarizedExperiment))
suppressMessages(library(EnrichedHeatmap))
suppressMessages(library(circlize))
suppressMessages(library(ggpubr))
suppressMessages(library(ggplot2))
suppressMessages(library(argparse))

parser = ArgumentParser(description="Correlation plot of crosslink signal between batches")
parser$add_argument('--bs-bed-path', help='Path to the .BED file (BED6 format) containing binding sites')
parser$add_argument('--batch1-bw-path', help='Path to the .bigwig file with crosslink events for the first batch')
parser$add_argument('--batch2-bw-path', help='Path to the .bigwig file with crosslink events for the second batch')
parser$add_argument('--out-filename', default="", help="Filename of the output figure")

args = parser$parse_args()

bs_path = args$bs_bed_path
roq_batch1_bw_path = args$batch1_bw_path
roq_batch2_bw_path = args$batch2_bw_path
out_filename = args$out_filename

# Load either batch1 bs or batch2 bs or the intersection bs
merged_bs = rtracklayer::import(con = bs_path, format = "BED")

# Load .bigwig of merged (crosslinks)
roq_batch1_bw = rtracklayer::import(con = roq_batch1_bw_path, format = "bw")
# Load .bigwig of batch2merged (crosslinks)
roq_batch2_bw = rtracklayer::import(con = roq_batch2_bw_path, format = "bw")

# Get signal of batch1merged at bs positions
b1_mat = normalizeToMatrix(roq_batch1_bw, merged_bs, value_column = "score", 
                           extend = c(0, 0), mean_mode = "w0", w = 1,
                           background = 0, smooth = FALSE)

# dim(b1_mat)

# Get signal of batch2merged at bs positions
b2_mat = normalizeToMatrix(roq_batch2_bw, merged_bs, value_column = "score", 
                           extend = c(0, 0), mean_mode = "w0", w = 1,
                           background = 0, smooth = FALSE)

# Compute sum of xlinks for each row
b1_sums = rowSums(b1_mat)
b2_sums = rowSums(b2_mat)

# Combine signals in a df
df = data.frame(b1_xlinks = b1_sums, b2_xlinks = b2_sums)

# Use quantiles to set the upper limit for the plot
# quantile(roq_batch1_bw$score, c(0, 0.25, 0.5, 0.75,  0.999))
# quantile(roq_batch2_bw$score, c(0, 0.25, 0.5, 0.75,  0.999))

upper_limit = 15
# sum((log2(df$b1_xlinks) > upper_limit) | (log2(df$b2_xlinks) > upper_limit)) # Check how many points "outliers" we're cutting from the plot

# Plot scatter using batch1merged signal as x-axis and batch2merged signal as y-axis
# For contour plots
# https://r-charts.com/correlation/contour-plot-ggplot2/?utm_content=cmp-true
plot = ggplot(log2(df + 1), aes(x = b1_xlinks, y = b2_xlinks)) +
  geom_point(shape = 21, colour = "black", fill = "white", size = 1, stroke = 1) +
  # geom_density_2d(aes(color = ..level..)) +
  # geom_density_2d_filled() +
  geom_abline(slope = 1, intercept = 0, colour = "lightgrey", linetype = "dashed", linewidth = 0.8) +
  stat_cor(method = "pearson") +
  labs(x = "Batch 1 crosslink counts (log2)", 
       y = "Batch 2 crosslink counts (log2)", 
       title = "Signal correlation between batches\n at unified binding sites") +
  xlim(0, upper_limit) +
  ylim(0, upper_limit) +
  # scale_y_continuous(expand = c(0, 0)) + # Removes space between bars and x-axis - this removes the previously set y scale
  theme_classic(base_size = 14) +
  theme(
    # text = element_text(family = "Arial"),
    plot.title = element_text(size = 14),  # Set font size for the title
    legend.text = element_text(size = 14, color = 'black'),  # Adjusts size of legend text
    legend.title = element_text(size = 14, color = 'black'),  # Adjusts size of legend title
    axis.title.x = element_text(size = 14, color = 'black'),  # X-axis label font size
    axis.title.y = element_text(size = 14, color = 'black'),  # Y-axis label font size
    axis.text.x = element_text(size = 14, color = 'black'),   # X-axis tick labels font size
    axis.text.y = element_text(size = 14, color = 'black'),   # Y-axis tick labels font size
    axis.line.x = element_line(linewidth = 1.5),  # Thicker x-axis line
    axis.line.y = element_line(linewidth = 1.5),   # Thicker y-axis line
    # legend.position = 'top'
    )

# plot # Show plot

# Save the plot in both PNG and PDF formats

if (out_filename == "") {
  png_filename = paste0("./figures/scatter_plot_batches_crosslinks.png")
  pdf_filename = paste0("./figures/scatter_plot_batches_crosslinks.pdf")
} else {
  png_filename = paste0("./figures/", out_filename, ".png")
  pdf_filename = paste0("./figures/", out_filename, ".pdf")
}
ggsave(png_filename, plot, width = 5, height = 5)
ggsave(pdf_filename, plot, width = 5, height = 5)


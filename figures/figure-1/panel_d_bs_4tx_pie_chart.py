""" Read a .bed file of binding sites and plot a pie chart divided by biotype """

import os, argparse
import pandas as pd 

import matplotlib
matplotlib.use('Agg') # don't use XWindows - this has to go before importing pyplot

import matplotlib.pyplot as plt

plt.rcParams['font.sans-serif'] = "Arial"
plt.rcParams['font.family'] = "sans-serif"

# base_color_mouse = "#D41159"
# base_color_human = "#1A85FF"

# colors_dict = {
#     "indigo": "#332288",
#     "cyan": "#88CCEE",
#     "teal": "#44AA99",
#     "green": "#117733",
#     "olive": "#999933",
#     "sand": "#DDCC77",
#     "rose": "#CC6677",
#     "wine": "#882255",
# }

colors_mouse = ["#D41159", "#B12E7D", "#82448D", "#524D89", "#314E74"]
colors_human = ["#1A85FF", "#00AEFF", "#00CDF1", "#00E4BB", "#8AF387"]

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot a pie chart of binding sites by transcript region")

    parser.add_argument('--input-bed', type=str, help="Path to annotated binding sites in .bed format")
    parser.add_argument('--outpath', type=str, help="Path to save pie chart")
    parser.add_argument('--palette', type=str, default="mouse", help="Select between two possible color palettes, [mouse|human]")

    args = parser.parse_args()
    print("############## Arguments ##############")
    print(args)
    print("#######################################")

    # bs = pd.read_csv(args.input_bed, sep="\t", names=["chr", "start", "end", "name"])
    bs = pd.read_csv(args.input_bed, sep="\t", names=["chr", "start", "end", "name", "score", "strand"])

    bs['tx_region'] = bs['name'].apply(lambda x: x.split(";")[0])
    bs = bs.groupby("tx_region")['chr'].count()
    print(bs.head(n=10))
    
    # We rm exon, stop_codon, start_codon cause we're not interested in these annotations
    bs = bs.drop(["exon", "start_codon", "stop_codon"])
    print(bs.head(n=10))

    renamed_columns = ["CDS", "5'UTR", "Intron", "3'UTR"]
    bs = bs.rename(index=dict(zip(bs.index, renamed_columns)))
    # Labels are dataframe groups
    #labels = bs.index
    labels = renamed_columns
    # Add absolute counts for each category to the legend
    labels = [f"{label} (n={bs[label]})" for label in labels]
    
    # Set color palette
    colors = colors_mouse
    # colors = [colors_dict[color] for color in colors]

    if args.palette == "human":
        colors = colors_human
    
    # Set figure size
    # plt.figure(figsize = (10, 10)) - if legend inside
    plt.figure(figsize = (20, 20)) # if legend outside, lower right 

    plt.pie(bs, 
        colors=colors, 
        autopct='%.0f%%', # percent format
        pctdistance=1.15, # distance of the percentage from pie center
        textprops={'fontsize': 45})
    
    # Add legend
    # lgd = plt.legend(labels, loc="lower right", fontsize=50) # place legend inside, lower right
    #lgd = plt.legend(labels, loc="right", fontsize=50) # place legend inside, right
    lgd = plt.legend(labels, loc="right", bbox_to_anchor=(1.2, 1), fontsize=35) # place legend outside plot, setting distance manually with bbox_to_anchor


    #plt.savefig(args.outpath)
    plt.savefig(args.outpath, dpi=300, bbox_extra_artists=(lgd,), bbox_inches='tight') # avoid pyplot from cropping legend, this fixes legend automatically

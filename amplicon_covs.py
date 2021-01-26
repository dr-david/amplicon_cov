import argparse
import numpy as np
import pandas as pd
import os
import re
import gzip
import subprocess
import matplotlib.pyplot as plt
import seaborn as sns

def parse_args():
    '''Parsing of command line args'''
    parser = argparse.ArgumentParser(
            description="Script to calculate primer rebalancings according to november 2020 version 5 of the ARTIC V3 protocol for sars-cov-2 sequencing.",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument("-r", required=True, default=None, metavar='BED',
                        dest='bedfile_addr', type=str,
                        help="Bedfile of the articV3 primers, eg. from: \
                        https://raw.githubusercontent.com/artic-network/artic-ncov2019/master/primer_schemes/nCoV-2019/V3/nCoV-2019.bed")
    parser.add_argument("-s", required=False, metavar='TSV',
                        dest='samp_file', help="tsv file like samples.tsv.", 
                        default='/cluster/project/pangolin/working/samples.tsv')
    parser.add_argument("-f", required=False, metavar='PATH',
                        dest='samp_path', help="main path to samples", 
                        default='/cluster/project/pangolin/working/samples')
    parser.add_argument("-o", required=False, default=os.getcwd(),
                        metavar='PATH', dest='outdir',
                        help="Output directory")
    parser.add_argument("-p", dest='makeplots', help="Output plots.", action='store_true')
    parser.add_argument("-v", help="Verbose", action='store_true')


    return parser.parse_args()

def get_samples_paths(main_samples_path='/cluster/project/pangolin/working/samples', samplestsv='/cluster/project/pangolin/working/samples.tsv'):
    '''make list of sample paths by combining main path and samples.tsv'''
#    sam_names_list = []
    sam_paths_list = []
    with open(samplestsv, 'r') as f:
        for line in f:
            tmp = line.rstrip("\n").split("\t")
#             sam_names_list.append((tmp[0], tmp[1]))
            sam_paths_list.append(main_samples_path+"/"+tmp[0]+"/"+tmp[1]+"/alignments/coverage.tsv.gz")
    return sam_paths_list

def load_bedfile(bed="articV3primers.bed"):
    '''function to load a bed file of primers'''
    bedfile = pd.read_table(bed, header=None)
    bedfile["sense"] = [re.search("(LEFT|RIGHT)",i).group(1) for i in bedfile[3]]
    bedfile["primer_num"] = [int(re.search("_([0-9]+)_",i).group(1)) for i in bedfile[3]]
    bedfile["pool"] = [int(re.search("_([1-2])",i).group(1)) for i in bedfile[4]]
    bedfile = bedfile[[re.search("alt", i) is None for i in bedfile[3]]]
    # bedfile["alt"] = [re.search("(_alt[0-9]+)", i).group(1) if re.search("(_alt[0-9]+)", i) is not None else " " for i in bedfile[3]]
    # bedfile["primer_code"] = bedfile["primer_num"].astype(str) + bedfile["alt"]
    return bedfile



def make_amplicons_df(bedfile):
    '''function to collapse loaded bedfile into a list of amplicons with start and stop positions of primers, sequences and query'''
    amplicons = []
    for i in np.unique(bedfile["primer_num"]):
        pr_num = i
        seq_start = bedfile[(bedfile["primer_num"] == pr_num) & (bedfile["sense"] == "LEFT")][2].values[0]
        primer_start = bedfile[(bedfile["primer_num"] == pr_num) & (bedfile["sense"] == "LEFT")][1].values[0]
        seq_end = bedfile[(bedfile["primer_num"] == pr_num) & (bedfile["sense"] == "RIGHT")][1].values[0]
        primer_end = bedfile[(bedfile["primer_num"] == pr_num) & (bedfile["sense"] == "RIGHT")][2].values[0]
        pool = bedfile[bedfile["primer_num"] == pr_num]["pool"].values[1]
        
        amplicons.append([pool, pr_num, primer_start, seq_start, seq_end, primer_end])

    amplicons_df = pd.DataFrame(np.array(amplicons),
                                columns=["pool", "primer_num", "primer_start", "seq_start", "seq_end", "primer_end"])
    
    # make query_start and query_stop
    q_starts = []
    q_stops = []
    for i in range(amplicons_df.shape[0]):
        if i>0:
            query_start = amplicons_df.iloc[i-1]["primer_end"] + 5
        else:
            query_start = amplicons_df.iloc[i]["primer_start"]
        
        if i < amplicons_df.shape[0]-1:
            query_stop = amplicons_df.iloc[i+1]["primer_start"] - 5
        else:
            query_stop = amplicons_df.iloc[i]["seq_end"]
            
        q_starts.append(query_start)
        q_stops.append(query_stop)

    amplicons_df["query_start"] = q_starts
    amplicons_df["query_end"] = q_stops

    return amplicons_df


def get_amplicon_cov(cov_df, start, stop, length=20):
    '''function to compute the median coverage in a start:stop positions slice of a cov_df'''
    amplicon_slice = cov_df.iloc[np.r_[start:length, (stop-length):stop],[2]]
    return np.median(amplicon_slice)


def get_frac_reads(cov_df, amplicons_df):
    '''function to return estimated fraction of the reads in cov_df aligned in each query window of the amplicon df'''
    cov = amplicons_df.apply(lambda x: get_amplicon_cov(cov_df, x["query_start"], x["query_end"]), axis=1)
    frac_reads = cov / np.sum(cov)

    return frac_reads


def make_cov_heatmap(cov_df, output=None):

    plt.figure(figsize=(15,8*2.5))

    split_at = round(cov_df.shape[0]/2)

    plt.subplot(1,2,1)
    ax = sns.heatmap(cov_df.iloc[0:split_at,1:], cmap='Reds', vmin=0, square=True,
                     cbar_kws={"shrink": .2, "anchor": (0.0, 0.8)})
    sns.heatmap(cov_df.iloc[0:split_at,1:],
                cmap=plt.get_cmap('binary'), vmin=0, vmax=2, mask=cov_df.iloc[0:split_at,1:] > 0, cbar=False, ax=ax)
    plt.xlabel("amplicon")
    plt.ylabel("sample")
    plt.title("Samples 0:{}".format(split_at))

    plt.subplot(1,2,2)
    ax = sns.heatmap(cov_df.iloc[split_at:,1:], cmap='Reds', vmin=0, square=True,
                     cbar_kws={"shrink": .2, "anchor": (0.0, 0.8)})
    sns.heatmap(cov_df.iloc[split_at:,1:],
                cmap=plt.get_cmap('binary'), vmin=0, vmax=2, mask=cov_df.iloc[split_at:,1:] > 0, cbar=False, ax=ax)
    plt.xlabel("amplicon")
    plt.ylabel("sample")
    plt.title("Samples {}:{}".format(split_at, cov_df.shape[0]-1))
    
    if output is not None:
        plt.savefig(output)

def make_median_cov_hist(cov_df, output=None):
    median = np.nanmedian(cov_df.iloc[:,1:].values, axis=0)
    
    plt.figure(figsize=(12,6))
    sns.histplot(y=median, binwidth=0.002, stat="density")
    plt.title("Median coverage histogram")
    plt.ylabel("median fraction of reads aligned on amplicon")
    plt.xlabel("density")
#     plt.ylim((-0.005,0.1))
#     plt.xlim((0,175))
    plt.axhline(1/98, linestyle="--", color="black")

    if output is not None:
        plt.savefig(output)

def make_median_coverage_barplot(cov_df, output=None):
    
    cov_df_long = pd.melt(cov_df.iloc[:,1:])
    cov_df_long["pool"] = cov_df_long["variable"].astype("int").mod(2) + 1
    
    plt.figure(figsize=(22, 9))
    sns.barplot(x="variable", y="value", hue="pool", data=cov_df_long, estimator=np.median)
    plt.axhline(1/98, linestyle="--", color="black")
    # plt.ylim((0, 0.1))
    plt.xlabel("amplicon")
    plt.ylabel("median fraction of reads")
    plt.title("Median coverage barplot")

    if output is not None:
        plt.savefig(output)


def main():
    # parse arguments
    args = parse_args()
    samp_file = args.samp_file
    samp_path = args.samp_path
    bedfile_addr = args.bedfile_addr
    outdir = args.outdir
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    # make amplicons df
    if args.v: print("Loading primers bedfile.")   
    amplicons_df = make_amplicons_df(load_bedfile(bedfile_addr))
    
    # read list of samples
    if args.v: print("Reading list of coverage files.")
    sam_list = get_samples_paths(samp_path, samp_file)

    # iterate through list of samples
    if args.v: print("Loading and parsing coverage files.")
    all_covs = []
    indexes = []
    i = 1
    for sam in sam_list:
        if args.v: print("Parsing coverage file {}/{}".format(i, len(sam_list)), end="\r")
        try:
            temp_cov_df = pd.read_csv(sam, sep="\t", compression="gzip")
            temp_frac_read_df = pd.DataFrame(get_frac_reads(temp_cov_df, amplicons_df)).T
            indexes.append(sam.split("/")[-4])
            all_covs.append(temp_frac_read_df)
        except FileNotFoundError:
            if args.v: print("WARNING: file {} not found.".format(sam))
    #         all_covs.append([])
        i += 1
    all_covs = pd.concat(all_covs, axis=0)
    all_covs = pd.concat([pd.DataFrame({"sample":indexes}), all_covs.reset_index(drop=True)], axis=1, ignore_index=False)
#    all_covs.set_index(pd.Index(indexes))
    
    # output DF
    if args.v: print("\nOutputting csv.")
    all_covs.to_csv(outdir + "/amplicons_coverages.csv", index = False)

    # make plots
    if args.makeplots: 
        if args.v: print("\nOutputting plots.")

        make_cov_heatmap(all_covs, outdir + "/cov_heatmap.pdf")
        make_median_cov_hist(all_covs, outdir + "/median_cov_hist.pdf")
        make_median_coverage_barplot(all_covs, outdir + "/median_coverage_barplot.pdf")



if __name__ == '__main__':
    main()


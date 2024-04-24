
from pathlib import Path
import numpy as np
import pandas as pd
import os
from tqdm import tqdm
from optparse import OptionParser
from collections import defaultdict

def read_diamond(infile, mode = "read_sum", translate_dict = None, min_identity=50, min_read_bases=140):
    cols = "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore".split(" ")
    if os.stat(infile).st_size == 0:
        return {}
    df = pd.read_csv(infile, sep = "\t", header=None)
    df.columns = cols
    df["read_bases"] = np.absolute(df["qend"] - df["qstart"])
    df = df.query('(pident >= {}) & (read_bases >= {})'.format(min_identity, min_read_bases))
    if translate_dict is not None:
        df["sseqid"] = [translate_dict[x] if x in translate_dict else x for x in df["sseqid"]]

    # Option #1, just sum the the unique number of reads
    if mode == "read_sum":
        return {"read_sum": len(df["qseqid"].unique())}
    elif mode == "best_hit":
        # Option #2, best hit
        results = defaultdict(int)
        for (read, dat) in df.groupby("qseqid"):
            dat = dat.sort_values("pident", ascending=False)
            #TODO: Discount cases where there is another close hit
            results[dat["sseqid"].values[0]] += 1
        return results

def load_dict(infile, col_oi = "Species"):
    if infile is None:
        return None
    df = pd.read_csv(infile)
    name_dict = {}
    for (i,dat) in df.iterrows():
        name_dict[dat["Protein"]] = dat[col_oi]
    return name_dict


if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("-m", "--mode", default="best_hit")
    parser.add_option("--min_identity", default=90)
    parser.add_option("--min_bases", default=90)
    parser.add_option("--tax_info", default="Species") #Only used if you include --infile_dict
    parser.add_option("--infile_dict") #something like Devlin_TnaA_referencespecies.csv 
    parser.add_option("-o", "--outfile")
    (options, args) = parser.parse_args()
    
    name_dict = load_dict(options.infile_dict, options.tax_info)

    master_dict = {}
    for infile in tqdm(args):
        infile = Path(infile)
        master_dict[infile.name] = read_diamond(infile, options.mode, name_dict, options.min_identity, options.min_bases)
    
    df_out = pd.DataFrame(master_dict).fillna(0)
    df_out.to_csv(options.outfile, index=True)


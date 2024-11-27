#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
'''
@File    :   card_trim.py
@Time    :   2024/11/18 17:07:42
@Author  :   Wanjin.Hu 
@Version :   1.0
@Contact :   wanjin.hu@outlook.com
@Description :
'''

import argparse
import pandas as pd
import os

def card_trim(card_raw,card_trim):
    "card results trim"
    with open(card_raw, 'r') as f, open(card_trim, 'w') as t:
        for lines in f:
            line = lines.strip().split("\t")
            orf_line = line[0].split(" ")
            orf_id = orf_line[0]
            other_info = line[5:]
            t.write(orf_id + "\t" + "\t".join(other_info) + "\n")
        t.close()

def card_merge(card_trim,merge_table,card_merge):
    "card results merge"
    df_card_trim = pd.read_csv(card_trim, header=0, sep="\t")
    df_merge = pd.read_csv(merge_table, header=0, sep="\t")
    df_card_merge = pd.merge(df_card_trim, df_merge, 
                             left_on="ORF_ID", right_on="gene", how="inner")
    df_card_merge.to_csv(card_merge, header=True, index=False, sep="\t")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog="card",
                                     usage="python card_trim.py -i card.txt -m merged_file.txt -o result",
                                     description="CARD results merge.")
    parser.add_argument('-i', '--input', dest="inCARD",
                        type=str, required=True, help="Sample's card information, such as card.txt")
    parser.add_argument('-m', '--merge_table', dest="mergeTable",
                        type=str, required=True, help="Sample's merged gene count table, such as merged_file.txt")
    parser.add_argument('-o', '--output', dest="outDir",
                        type=str, required=True, help="Output directory")
    args = parser.parse_args()
    if not os.path.exists(args.outDir):
        os.makedirs(args.outDir)
    card_trim(card_raw=args.inCARD,
              card_trim=os.path.join(args.outDir,"card_trim.txt"))
    card_merge(card_trim=os.path.join(args.outDir,"card_trim.txt"),
               merge_table=args.mergeTable,
               card_merge=os.path.join(args.outDir,"card_merge.txt"))

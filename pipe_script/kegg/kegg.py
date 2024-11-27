#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
'''
@File    :   kegg.py
@Time    :   2023/03/27 19:25:33
@Author  :   Wanjin.Hu 
@Version :   1.0
@Contact :   wanjin.hu@diprobio.com
@Description : 对使用 eggnog 方式获取的 kegg 信息进行处理，生成后续可以分析的 KO 基因表 (KO_samples.csv) 和 pathway 代谢通路表 (pathway_samples.csv)
'''

import pandas as pd
import os
import argparse

class KeggCount(object):
    def __init__(self):
        super(KeggCount,self).__init__()
        self.kegg_db = "/root/pipe_script/metagenome/pipe_script/kegg/KEGG_KO_ENZYME_PATHWAY.txt"    # 原始的 KO 对应的 KEGG 信息表【功能：统计 KO 基因】
        self.kegg_pathway = "/root/pipe_script/metagenome/pipe_script/kegg/ko_pathway_uniq.txt"  # 处理过的 ko 代谢通路三个层级的信息表【功能：统计 pathway 代谢通路】
        # self.KO = r"F:\Dipro-Project\metagenome\KEGG_KO.txt"    # 宏基因组分析得到的非冗余基因集对应的 KO 基因信息表【示例】
        # self.pathway = r"F:\Dipro-Project\metagenome\KEGG_PATHWAY.txt"  # 宏基因组分析得到的非冗余基因集对应的 ko 的通路信息表【示例】
        # self.sample_merge = r"F:\Dipro-Project\metagenome\merged_file.txt"  # 宏基因组分析得到的非冗余基因集在样本中的分布表【示例】
        if not os.path.exists(args.tmp):
            os.mkdir(args.tmp)
        
    def KO_count(self,df_gene,df_ko):
        """
        @功能：生成 KO 基因组成表
        @df_gene: 宏基因组分析得到的非冗余基因集在样本中的分布表
        @df_ko: 宏基因组分析得到的非冗余基因集对应的 KO 基因信息表
        """
        df_gene_ko = pd.merge(df_gene, df_ko, on="gene", how="inner")
        if os.path.exists("df_gene_ko.txt"):
            os.remove("df_gene_ko.txt")
            df_gene_ko.to_csv("df_gene_ko.txt", header=True, index=False, sep="\t")
        if not os.path.exists("df_gene_ko.txt"):
            df_gene_ko.to_csv("df_gene_ko.txt", header=True, index=False, sep="\t")
        with open("df_gene_ko.txt", 'r') as f1, open(self.kegg_db, 'r') as f2, open('tmp_KO.txt', 'w') as f3:
            sample_line = f1.readline()
            db_line = f2.readline()
            sample_line = sample_line.strip().split("\t")
            a_dict = {}
            for line in f1:
                line = line.strip().split('\t')
                a_dict[line[0]] = line[1:]
            b_dict = {}
            for line in f2:
                line = line.strip().split('\t')
                b_dict[line[0]] = line[1:3]
            # KO_name KO_des KO  [sample1,sample2...]
            sam_list = sample_line[1:-1]
            # f3.write("meta_gene\tKO_name\tKO_des\tKO"+"\t"+"\t".join(sam_list) + "\n")
            f3.write("KO_name\tKO_des\tKO"+"\t"+"\t".join(sam_list) + "\n")
            for k,v in a_dict.items():
                KO_ids = v[-1]
                if "," in KO_ids:
                    KO_ids = KO_ids.split(",")
                    for KO_id in KO_ids:
                        if KO_id in b_dict:
                            res_line = b_dict[KO_id][0]+"\t"+b_dict[KO_id][1]+"\t"+KO_id+"\t"+"\t".join(v[0:-1])
                            f3.write(res_line+"\n")
                        else:
                            continue
                            # print("{} not in kegg_db".format(KO_id))
                else:
                    if KO_ids in b_dict:
                        res_line = b_dict[KO_ids][0]+"\t"+b_dict[KO_ids][1]+"\t"+KO_ids+"\t"+"\t".join(v[0:-1])
                        f3.write(res_line+"\n")
                    else:
                        continue
                        # print("{} not in kegg_db".format(KO_ids))
            f3.close()
        # tmp_KO.txt 包含了所有冗余的 KO 信息，这一步需要对相同的 KO 求和合并
        KO_tmp = pd.read_csv("tmp_KO.txt", header=0, sep="\t")
        cols = KO_tmp.columns[2:]
        KO_tmp_1 = KO_tmp.groupby("KO")[cols].sum(numeric_only=True)
        # 这里用了一个笨方法，重新打开文件进行匹配
        KO_tmp_1.to_csv("KO.txt", header=True, index=True, sep="\t")
        with open("KO.txt", 'r') as f4, open(self.kegg_db, 'r') as f5, open(args.outKO,'w') as f6:
            sample_line = f4.readline()
            sample_line = sample_line.strip().split("\t")
            db_line = f5.readline()
            a_dict = {}
            for line in f4:
                line = line.strip().split('\t')
                a_dict[line[0]] = line[1:]
            b_dict = {}
            for line in f5:
                line = line.strip().split('\t')
                b_dict[line[0]] = line[1:3]
            sam_list = sample_line[1:]
            f6.write("KO_name\tKO_des\tKO"+"\t"+"\t".join(sam_list) + "\n")
            for k,v in a_dict.items():
                KO_ids = k
                if KO_ids in b_dict:
                    res_line = b_dict[KO_ids][0]+"\t"+b_dict[KO_ids][1]+"\t"+KO_ids+"\t"+"\t".join(v)
                    f6.write(res_line+"\n")
                else:
                    continue
            f6.close()
        os.system("mv df_gene_ko.txt tmp_KO.txt KO.txt {}".format(args.tmp))
    
    def path_conut(self,df_gene,df_pathway):
        """
        @功能：生成代谢通路 pathway 组成表
        @df_gene: 宏基因组分析得到的非冗余基因集在样本中的分布表
        @df_pathway: 宏基因组分析得到的非冗余基因集对应的 ko 的通路信息表
        """
        df_gene_pathway = pd.merge(df_gene, df_pathway, on="gene", how="inner")
        if os.path.exists("df_gene_pathway.txt"):
            os.remove("df_gene_pathway.txt")
            df_gene_pathway.to_csv("df_gene_pathway.txt", header=True, index=False, sep="\t")
        if not os.path.exists("df_gene_pathway.txt"):
            df_gene_pathway.to_csv("df_gene_pathway.txt", header=True, index=False, sep="\t")
        with open("df_gene_pathway.txt", 'r') as f1, open(self.kegg_pathway, 'r') as f2, open('tmp_pathway.txt', 'w') as f3:
            sample_line = f1.readline()
            db_line = f2.readline()
            sample_line = sample_line.strip().split("\t")
            a_dict = {}
            for line in f1:
                line = line.strip().split('\t')
                a_dict[line[0]] = line[1:]
            b_dict = {}
            for line in f2:
                line = line.strip().split('\t')
                b_dict[line[0]] = line[1:4]
            # level1 level2 level3 pathway [sample1,sample2...]
            sam_list = sample_line[1:-1]
            # f3.write("meta_gene\tKO_name\tKO_des\tKO"+"\t"+"\t".join(sam_list) + "\n")
            f3.write("level1\tlevel2\tlevel3\tpathway"+"\t"+"\t".join(sam_list) + "\n")
            for k,v in a_dict.items():
                path_ids = v[-1]
                if "," in path_ids:
                    path_ids = path_ids.split(",")
                    for path_id in path_ids:
                        if path_id.startswith("map"):
                            continue
                        else:
                            if path_id in b_dict:
                                res_line = b_dict[path_id][0]+"\t"+b_dict[path_id][1]+"\t"+b_dict[path_id][2]+"\t"+path_id+"\t"+"\t".join(v[0:-1])
                                f3.write(res_line+"\n")
                            else:
                                continue
                                # print("{} not in kegg_db".format(path_id))
                else:
                    if path_ids.startswith("map"):
                        continue
                    else:
                        if path_ids in b_dict:
                            res_line = b_dict[path_ids][0]+"\t"+b_dict[path_ids][1]+"\t"+b_dict[path_ids][2]+"\t"+path_ids+"\t"+"\t".join(v[0:-1])
                            f3.write(res_line+"\n")
                        else:
                            continue
                            # print("{} not in kegg_db".format(path_ids))
            f3.close()
        # tmp_pathway.txt 包含了所有冗余的 pathway 信息，这一步需要对相同的 pathway 求和合并
        pathway_tmp = pd.read_csv("tmp_pathway.txt", header=0, sep="\t")
        cols = pathway_tmp.columns[3:]
        KO_tmp_1 = pathway_tmp.groupby("pathway")[cols].sum(numeric_only=True)
        # 这里用了一个笨方法，重新打开文件进行匹配
        KO_tmp_1.to_csv("kegg.txt", header=True, index=True, sep="\t")
        with open("kegg.txt", 'r') as f4, open(self.kegg_pathway, 'r') as f5, open(args.outPathway,'w') as f6:
            sample_line = f4.readline()
            sample_line = sample_line.strip().split("\t")
            db_line = f5.readline()
            a_dict = {}
            for line in f4:
                line = line.strip().split('\t')
                a_dict[line[0]] = line[1:]
            b_dict = {}
            for line in f5:
                line = line.strip().split('\t')
                b_dict[line[0]] = line[1:]
            sam_list = sample_line[1:]
            f6.write("level1\tlevel2\tlevel3\tpathway"+"\t"+"\t".join(sam_list) + "\n")
            for k,v in a_dict.items():
                path_ids = k
                if path_ids in b_dict:
                    res_line = b_dict[path_ids][0]+"\t"+b_dict[path_ids][1]+"\t"+b_dict[path_ids][2]+"\t"+path_ids+"\t"+"\t".join(v)
                    f6.write(res_line+"\n")
                else:
                    continue
            f6.close()
        os.system("mv df_gene_pathway.txt tmp_pathway.txt kegg.txt {}".format(args.tmp))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="kegg",
                                     usage="python kegg.py -kk KEGG_KO.txt -kp KEGG_PATHWAY.txt -mt merged_file.txt -ok out_KO.csv -op out_pathway.csv",
                                     description="Merge KO/pathway count table from eggnog result.")
    parser.add_argument('-kk', '--kegg_KO', dest="keggKO",
                        type=str, required=True, help="Sample's kegg KO information, such as KEGG_KO.txt")
    parser.add_argument('-kp', '--kegg_pathway', dest="keggPathway",
                        type=str, required=True, help="Sample's kegg pathway information, such as KEGG_PATHWAY.txt")
    parser.add_argument('-mt', '--merge_table', dest="mergeTable",
                        type=str, required=True, help="Sample's merged gene count table, such as merged_file.txt")
    parser.add_argument('-ok', '--out_KO', dest="outKO",
                        type=str, required=True, help="Output KO result, such as out_KO.xls")
    parser.add_argument('-op', '--out_pathway', dest="outPathway",
                        type=str, required=True, help="Output pathway result, such as out_pathway.xls")
    parser.add_argument('-t', '--tmp', dest="tmp",
                        type=str, required=False, default="tmp", help="Tmp files dir")
    args = parser.parse_args()
    df_gene = pd.read_csv(args.mergeTable, header=0, sep="\t")
    df_ko = pd.read_csv(args.keggKO, header=0, sep="\t")
    df_pathway = pd.read_csv(args.keggPathway, header=0, sep="\t")
    run = KeggCount()
    run.KO_count(df_gene=df_gene,df_ko=df_ko)
    run.path_conut(df_gene=df_gene,df_pathway=df_pathway)

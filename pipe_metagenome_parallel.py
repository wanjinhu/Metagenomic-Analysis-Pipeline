#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
'''
@File    :   pipe_metagenome_parallel.py
@Time    :   2024/06/25 14:15:32
@Author  :   Wanjin.Hu 
@Version :   2.0
@Contact :   wanjin.hu@diprobio.com
@Description : Parallel version of metagenome pipeline
'''

import argparse
import os
import pandas as pd
from multiprocessing import Pool, cpu_count
from pathlib import Path
import shutil
import subprocess

class PipeMetagenome(object):
    """宏基因组分析流程 - 并行版本
    """

    def __init__(self):
        super(PipeMetagenome).__init__()
        self.bash_env = "/usr/bin/bash"
        self.py_env = "/root/miniconda3/bin/python"
        self.fastp_cmd = "/root/pipe_script/metagenome/pipe_script/fastp.sh"
        self.ref_remove_cmd = "/root/pipe_script/metagenome/pipe_script/ref_remove.sh"
        self.ref_remove_d_cmd = "/root/pipe_script/metagenome/pipe_script/ref_remove_double.sh"
        self.metaphlan_cmd = "/root/pipe_script/metagenome/pipe_script/metaphlan.sh"
        self.metaphlan_merge_cmd = "/root/pipe_script/metagenome/pipe_script/metaphlan_merge.sh"
        self.metaphlan_merge = "/root/miniconda3/bin/merge_metaphlan_tables.py"
        self.megahit_cmd = "/root/pipe_script/metagenome/pipe_script/megahit.sh"
        self.prodigal_cmd = "/root/pipe_script/metagenome/pipe_script/prodigal.sh"
        self.cdhit_cmd = "/root/pipe_script/metagenome/pipe_script/cd-hit.sh"
        self.emapper_cmd = "/root/pipe_script/metagenome/pipe_script/emapper.sh"
        self.count_cmd = "/root/pipe_script/metagenome/pipe_script/gene_count.sh"
        self.count_bwa_cmd = "/root/pipe_script/metagenome/pipe_script/gene_count_bwa.sh"
        self.count_bowtie2_cmd = "/root/pipe_script/metagenome/pipe_script/gene_count_bowtie2.sh"
        self.count_salmon_cmd = "/root/pipe_script/metagenome/pipe_script/gene_count_salmon.sh"
        self.index_bwa_cmd = "/root/pipe_script/metagenome/pipe_script/index_bwa.sh"
        self.index_bowtie2_cmd = "/root/pipe_script/metagenome/pipe_script/index_bowtie2.sh"
        self.index_salmon_cmd = "/root/pipe_script/metagenome/pipe_script/index_salmon.sh"
        self.kegg_cmd = "/root/pipe_script/metagenome/pipe_script/kegg/kegg.py"
        self.cazy_cmd = "/root/pipe_script/metagenome/pipe_script/cazy/cazy_anno.py"
        self.diamond_cmd = "/root/pipe_script/general/diamond_blastp.sh"
        self.cazy_diamond_index = "/root/database/CAZy_dbCAN2_V12/CAZyDB.07262023"
        self.rgi_cmd = "/root/pipe_script/metagenome/pipe_script/card/card_rgi.sh"
        self.card_cmd = "/root/pipe_script/metagenome/pipe_script/card/card_trim.py"
        # 设置并行进程数，默认为CPU核心数的90%
        self.num_processes = max(1, int(cpu_count() * 0.90))
        # 设置环境
        self.env = os.environ.copy()  # 复制当前环境变量
        self.map_method = args.mapMethod  

    def run_in_mamba_env(self, cmd, env_name="base"):
        """在指定的 mamba 环境中运行命令"""
        full_cmd = f"""
            eval "$(micromamba shell hook --shell=bash)"
            micromamba activate {env_name}
            {cmd}
        """
        process = subprocess.Popen(
            full_cmd,
            shell=True,
            executable='/bin/bash',
            env=self.env,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )
        stdout, stderr = process.communicate()
        
        if process.returncode != 0:
            print(f"Error executing command: {cmd}")
            print(f"stderr: {stderr.decode()}")
            raise RuntimeError(f"Command failed with return code {process.returncode}")
        
        return stdout.decode()

    def run_in_conda_env(self, cmd, env_name="base"):
        """在指定的 conda 环境中运行命令"""
        full_cmd = f"""
            eval "$(conda shell.bash hook)"
            conda activate {env_name}
            {cmd}
        """
        process = subprocess.Popen(
            full_cmd,
            shell=True,
            executable='/bin/bash',
            env=self.env,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )
        stdout, stderr = process.communicate()
        
        if process.returncode != 0:
            print(f"Error executing command: {cmd}")
            print(f"stderr: {stderr.decode()}")
            raise RuntimeError(f"Command failed with return code {process.returncode}")
        
        return stdout.decode()

    def run_command(self, cmd):
        """在常规环境中运行命令"""
        process = subprocess.Popen(
            cmd,
            shell=True,
            executable='/bin/bash',
            env=self.env,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )
        stdout, stderr = process.communicate()
        
        if process.returncode != 0:
            print(f"Error executing command: {cmd}")
            print(f"stderr: {stderr.decode()}")
            raise RuntimeError(f"Command failed with return code {process.returncode}")
        
        return stdout.decode()
  
    def load_sample_info(self, args):
        """加载样本信息
        """
        self.raw_dict = {}
        with open(args.fqList) as f:
            for line in f:
                if not line.startswith("#"):
                    line = line.strip().split("\t")
                    self.raw_dict[line[0]] = line[1:]
    
    def process_single_fastp(self, args):
        """单个样本的fastp处理
        """
        sam_name, sam_info = args
        sam_r1, sam_r2 = sam_info
        sam_fastp_path = os.path.join(self.fastp_path, sam_name)
        
        if os.path.exists(sam_fastp_path):
            print(f"Fastp already trimming finished: {sam_name}")
            return
            
        os.makedirs(sam_fastp_path, 0o755)
        sam_r1_clean = f"{sam_fastp_path}/{sam_name}_clean.1.fastq.gz"
        sam_r2_clean = f"{sam_fastp_path}/{sam_name}_clean.2.fastq.gz"
        sam_html = f"{sam_fastp_path}/{sam_name}.html"
        sam_json = f"{sam_fastp_path}/{sam_name}.json"
        
        cmd = f"{self.bash_env} {self.fastp_cmd} " \
              f"{sam_r1} {sam_r1_clean} {sam_r2} {sam_r2_clean} " \
              f"{sam_html} {sam_json}"
        try:
            self.run_command(cmd)
            print(f"Fastp finished successfully for sample: {sam_name}")
        except RuntimeError as e:
            print(f"Error fastp processing for sample {sam_name}: {str(e)}")
        
        # delete raw data
        # if os.path.exists(sam_r1_clean) and os.path.exists(sam_r2_clean):
        #     os.system(command=f"rm -rf {sam_r1} {sam_r2}")

    def fastp_trim(self, args):
        """并行运行fastp质控
        """
        self.outDir = args.outDir
        self.fastp_path = os.path.join(args.outDir, "01-fastp_trim")
        if not os.path.exists(self.fastp_path):
            os.makedirs(self.fastp_path, 0o755)
        
        with Pool(self.num_processes) as pool:
            pool.map(self.process_single_fastp, self.raw_dict.items())

    def process_single_ref_remove(self, args):
        """单个样本的去宿主处理
        """
        sam_name, sam_info = args
        sam_ref_remove_path = os.path.join(self.ref_remove_path, sam_name)
        
        if os.path.exists(sam_ref_remove_path):
            print(f"Ref sequence remove already finished: {sam_name}")
            return
            
        os.makedirs(sam_ref_remove_path, 0o755)
        sam_r1_clean = f"{self.fastp_path}/{sam_name}/{sam_name}_clean.1.fastq.gz"
        sam_r2_clean = f"{self.fastp_path}/{sam_name}/{sam_name}_clean.2.fastq.gz"
        sam_ref_sam = f"{sam_ref_remove_path}/{sam_name}.sam"
        sam_ref_sam_log = f"{sam_ref_remove_path}/{sam_name}.mapping.log"
        sam_r1_unmap = f"{sam_ref_remove_path}/{sam_name}.unmap.1.fastq.gz"
        sam_r2_unmap = f"{sam_ref_remove_path}/{sam_name}.unmap.2.fastq.gz"
        sam_unmap_single = f"{sam_ref_remove_path}/{sam_name}.unmap.single.fastq.gz"

        if self.ref1:
            cmd = f"{self.bash_env} {self.ref_remove_d_cmd} " \
                  f"{self.ref} {sam_r1_clean} {sam_r2_clean} " \
                  f"{sam_ref_sam} {sam_ref_sam_log} {sam_r1_unmap} {sam_r2_unmap} {sam_unmap_single} {self.ref1}"
        else:
            cmd = f"{self.bash_env} {self.ref_remove_cmd} " \
                  f"{self.ref} {sam_r1_clean} {sam_r2_clean} " \
                  f"{sam_ref_sam} {sam_ref_sam_log} {sam_r1_unmap} {sam_r2_unmap} {sam_unmap_single}"
        cmd1 = f"rm -rf {sam_r1_clean} {sam_r2_clean} {sam_ref_sam} {sam_unmap_single} process*"
        
        try:
            self.run_command(cmd)
            self.run_command(cmd1)
            print(f"Reference sequence removal finished successfully for sample: {sam_name}")
        except RuntimeError as e:
            print(f"Error removing reference for sample {sam_name}: {str(e)}")
        
    def ref_remove(self, args):
        """并行运行去宿主
        """
        self.ref = args.ref
        self.ref1 = args.ref1
        self.ref_remove_path = os.path.join(args.outDir, "02-ref_remove")
        if not os.path.exists(self.ref_remove_path):
            os.makedirs(self.ref_remove_path, 0o755)
            
        with Pool(self.num_processes) as pool:
            pool.map(self.process_single_ref_remove, self.raw_dict.items())

    def process_single_metaphlan(self, args):
        """单个样本的metaphlan处理[暂不使用]
        """
        sam_name, sam_info = args
        sam_metaphlan_path = os.path.join(self.metaphlan_path, sam_name)
        
        if os.path.exists(sam_metaphlan_path):
            print(f"Metaphlan already finished: {sam_name}")
            return
            
        os.makedirs(sam_metaphlan_path, 0o755)
        sam_r1_unmap = f"{self.ref_remove_path}/{sam_name}/{sam_name}.unmap.1.fastq.gz"
        sam_r2_unmap = f"{self.ref_remove_path}/{sam_name}/{sam_name}.unmap.2.fastq.gz"
        sam_metaphlan_bz2 = f"{sam_metaphlan_path}/{sam_name}_bowtie2.bz2"
        sam_metaphlan_out = f"{sam_metaphlan_path}/{sam_name}.tsv"
        
        cmd = f"{self.bash_env} {self.metaphlan_cmd} " \
              f"{sam_r1_unmap} {sam_r2_unmap} {sam_metaphlan_bz2} {sam_metaphlan_out}"
        cmd1 = f"rm -rf {sam_metaphlan_bz2}"
        self.run_in_conda_env(cmd)
        self.run_command(cmd1)
        # os.system(command=cmd)
        # os.system(command=f"rm -rf {sam_metaphlan_bz2}")

    def metaphlan_parallel(self, args):
        """并行运行metaphlan物种注释[并发会出问题，暂不使用]
        """
        self.metaphlan_path = os.path.join(args.outDir, "03-metaphlan")
        if not os.path.exists(self.metaphlan_path):
            os.makedirs(self.metaphlan_path, 0o755)
            
        with Pool(self.num_processes) as pool:
            pool.map(self.process_single_metaphlan, self.raw_dict.items())
            
        # 合并metaphlan结果
        metaphlan_merge = f"{self.metaphlan_path}/00_merged_abundance_table.txt"
        if not os.path.exists(metaphlan_merge):
            cmd1 = f"{self.metaphlan_merge} {self.metaphlan_path}/*/*.tsv > {metaphlan_merge}"
            self.run_in_conda_env(cmd1)
            # os.system(command=cmd1)
            
            metaphlan_files = {
                'phylum': f"{self.metaphlan_path}/01_metaphlan_phylum.txt",
                'class': f"{self.metaphlan_path}/02_metaphlan_class.txt",
                'order': f"{self.metaphlan_path}/03_metaphlan_order.txt",
                'family': f"{self.metaphlan_path}/04_metaphlan_family.txt",
                'genus': f"{self.metaphlan_path}/05_metaphlan_genus.txt",
                'species': f"{self.metaphlan_path}/06_metaphlan_species.txt"
            }
            cmd2 = f"{self.bash_env} {self.metaphlan_merge_cmd} " \
                   f"{metaphlan_merge} {' '.join(metaphlan_files.values())}"
            self.run_command(cmd2)
            # os.system(command=cmd2)

    def metaphlan(self,args):
        """metaphlan 物种注释"""
        self.metaphlan_path = os.path.join(args.outDir,"03-metaphlan")
        if not os.path.exists(self.metaphlan_path):
            os.makedirs(self.metaphlan_path, 0o755)
        for k,v in self.raw_dict.items():
            sam_name = k
            sam_metaphlan_path = os.path.join(self.metaphlan_path,sam_name)
            
            if os.path.exists(sam_metaphlan_path):
                print(f"Metaphlan already finished: {sam_name}")
                return
                
            os.makedirs(sam_metaphlan_path, 0o755)
            sam_r1_unmap = f"{self.ref_remove_path}/{sam_name}/{sam_name}.unmap.1.fastq.gz"
            sam_r2_unmap = f"{self.ref_remove_path}/{sam_name}/{sam_name}.unmap.2.fastq.gz"
            sam_metaphlan_bz2 = f"{sam_metaphlan_path}/{sam_name}_bowtie2.bz2"
            sam_metaphlan_out = f"{sam_metaphlan_path}/{sam_name}.tsv"
            
            cmd = f"{self.bash_env} {self.metaphlan_cmd} " \
                f"{sam_r1_unmap} {sam_r2_unmap} {sam_metaphlan_bz2} {sam_metaphlan_out}"
            cmd1 = f"rm -rf {sam_metaphlan_bz2}"
            try:
                self.run_in_conda_env(cmd)
                self.run_command(cmd1)
                print(f"Metaphlan finished successfully for sample: {sam_name}")
            except RuntimeError as e:
                print(f"Error metaphlan for sample {sam_name}: {str(e)}")

        # 合并metaphlan结果
        metaphlan_merge = f"{self.metaphlan_path}/00_merged_abundance_table.txt"
        if not os.path.exists(metaphlan_merge):
            cmd1 = f"{self.metaphlan_merge} {self.metaphlan_path}/*/*.tsv > {metaphlan_merge}"
            try:
                self.run_in_conda_env(cmd1)
                print(f"Metaphlan merge finished successfully: {metaphlan_merge}")
            except RuntimeError as e:
                print(f"Error metaphlan merge: {str(e)}")  
            metaphlan_files = {
                'phylum': f"{self.metaphlan_path}/01_metaphlan_phylum.txt",
                'class': f"{self.metaphlan_path}/02_metaphlan_class.txt",
                'order': f"{self.metaphlan_path}/03_metaphlan_order.txt",
                'family': f"{self.metaphlan_path}/04_metaphlan_family.txt",
                'genus': f"{self.metaphlan_path}/05_metaphlan_genus.txt",
                'species': f"{self.metaphlan_path}/06_metaphlan_species.txt"
            }
            cmd2 = f"{self.bash_env} {self.metaphlan_merge_cmd} " \
                   f"{metaphlan_merge} {' '.join(metaphlan_files.values())}"
            try:
                self.run_command(cmd2)
                print(f"Metaphlan split finished successfully: {metaphlan_merge}")
            except RuntimeError as e:
                print(f"Error metaphlan split: {str(e)}")
        else:
            print(f"Metaphlan merge finished successfully: {metaphlan_merge}")   

    def process_single_megahit(self, args):
        """单个样本的megahit处理
        """
        sam_name, sam_info = args
        sam_megahit_path = os.path.join(self.megahit_path, sam_name)
        
        if os.path.exists(sam_megahit_path):
            print(f"Megahit already finished: {sam_name}")
            os.system(f"rm -rf {sam_megahit_path}/intermediate_contigs")
            return
            
        sam_r1_unmap = f"{self.ref_remove_path}/{sam_name}/{sam_name}.unmap.1.fastq.gz"
        sam_r2_unmap = f"{self.ref_remove_path}/{sam_name}/{sam_name}.unmap.2.fastq.gz"
        sam_contig = f"{sam_megahit_path}/{sam_name}.contigs.fa"
        sam_contig_500 = f"{sam_megahit_path}/{sam_name}.contigs_500.fa"
        
        cmd = f"{self.bash_env} {self.megahit_cmd} " \
              f"{sam_r1_unmap} {sam_r2_unmap} {sam_megahit_path} " \
              f"{sam_name} {sam_contig} {sam_contig_500}"
        cmd1 = f"rm -rf {sam_megahit_path}/intermediate_contigs"

        try:
            self.run_command(cmd)
            self.run_command(cmd1)
            print(f"Megahit finished successfully for sample: {sam_name}")
        except RuntimeError as e:
            print(f"Error megahit for sample {sam_name}: {str(e)}")

    def megahit(self, args):
        """并行运行megahit序列组装
        """
        self.megahit_path = os.path.join(args.outDir, "04-megahit")
        if not os.path.exists(self.megahit_path):
            os.makedirs(self.megahit_path, 0o755)
            
        with Pool(self.num_processes) as pool:
            pool.map(self.process_single_megahit, self.raw_dict.items())

    def process_single_prodigal(self, args):
        """单个样本的prodigal处理
        """
        sam_name, sam_info = args
        sam_prodigal_path = os.path.join(self.prodigal_path, sam_name)
        
        if os.path.exists(sam_prodigal_path):
            print(f"Prodigal already finished: {sam_name}")
            return
            
        os.makedirs(sam_prodigal_path, 0o755)
        sam_prot = f"{sam_prodigal_path}/{sam_name}_prot.faa"
        sam_nucl = f"{sam_prodigal_path}/{sam_name}_nucl.fna"
        sam_gff = f"{sam_prodigal_path}/{sam_name}_genes.gff"
        sam_stat = f"{sam_prodigal_path}/{sam_name}.stat"
        sam_contig_500 = f"{self.megahit_path}/{sam_name}/{sam_name}.contigs_500.fa"
        
        cmd = f"{self.bash_env} {self.prodigal_cmd} " \
              f"{sam_prot} {sam_nucl} {sam_gff} {sam_stat} {sam_contig_500}"
        try:
            self.run_command(cmd)
            print(f"Prodigal finished successfully for sample: {sam_name}")
        except RuntimeError as e:
            print(f"Error prodigal for sample {sam_name}: {str(e)}")

    def prodigal(self, args):
        """并行运行prodigal基因预测
        """
        self.prodigal_path = os.path.join(args.outDir, "05-prodigal")
        if not os.path.exists(self.prodigal_path):
            os.makedirs(self.prodigal_path, 0o755)
            
        with Pool(self.num_processes) as pool:
            pool.map(self.process_single_prodigal, self.raw_dict.items())

    def cdhit(self, args):
        """cdhit去冗余，并对非冗余基因集构建索引
        """
        self.cdhit_path = os.path.join(args.outDir, "06-cdhit")
        if not os.path.exists(self.cdhit_path):
            os.makedirs(self.cdhit_path, 0o755)
        
        prot_cat = f"{self.cdhit_path}/prot.faa"
        nucl_cat = f"{self.cdhit_path}/nucl.fna"
        prot_nonerude = f"{self.cdhit_path}/prot_nonerude.faa"
        nucl_nonerude = f"{self.cdhit_path}/nucl_nonerude.fna"
        nonerude_list = f"{self.cdhit_path}/prot_nonerude.list"
        geneset_length = f"{self.cdhit_path}/geneset_length.txt"
        geneset = f"{self.cdhit_path}/index_geneset"
        
        if os.path.exists(geneset_length) and os.path.exists(geneset):
            print("cd-hit and geneset index already finished.")
            return
            
        # 合并各个样本的基因集
        os.system(command=f"cat {self.prodigal_path}/*/*_prot.faa > {prot_cat}")
        os.system(command=f"cat {self.prodigal_path}/*/*_nucl.fna > {nucl_cat}")
        
        # 运行cd-hit
        cmd = f"{self.bash_env} {self.cdhit_cmd} " \
              f"{prot_cat} {prot_nonerude} {nonerude_list} " \
              f"{nucl_cat} {nucl_nonerude} {geneset_length}"
        try:
            self.run_command(cmd)
            print(f"cd-hit finished successfully.")
        except RuntimeError as e:
            print(f"Error cd-hit: {str(e)}")
        
        # 构建非冗余基因集的索引
        if self.map_method == "bwa":
            cmd_index = f"{self.bash_env} {self.index_bwa_cmd} {nucl_nonerude} {geneset}"
            try:
                self.run_command(cmd_index)
                print(f"Geneset index for bwa is successfully.")
            except RuntimeError as e:
                print(f"Error geneset index for bwa: {str(e)}")
            
        if self.map_method == "bowtie2":
            cmd_index = f"{self.bash_env} {self.index_bowtie2_cmd} {nucl_nonerude} {geneset}"
            try:
                self.run_command(cmd_index)
                print(f"Geneset index for bowtie2 is successfully.")
            except RuntimeError as e:
                print(f"Error geneset index for bowtie2: {str(e)}")            

        if self.map_method == "salmon":
            cmd_index = f"{self.bash_env} {self.index_salmon_cmd} {nucl_nonerude} {geneset}"
            try:
                self.run_in_mamba_env(cmd_index,env_name="salmon")   
                print(f"Geneset index for salmon is successfully.")
            except RuntimeError as e:
                print(f"Error geneset index for salmon: {str(e)}")
        
    def emapper(self, args):
        """emapper功能注释，生成KO/KEGG信息表
        """
        self.emapper_path = os.path.join(args.outDir, "07-emapper")
        if not os.path.exists(self.emapper_path):
            os.makedirs(self.emapper_path, 0o755)
        prot_nonerude = f"{self.cdhit_path}/prot_nonerude.faa"
        eggnog_profile = f"{self.emapper_path}/eggnog"
        eggnog_ann = f"{self.emapper_path}/eggnog.emapper.annotations"
        eggnog_KO = f"{self.emapper_path}/KEGG_KO.txt"
        eggnog_path = f"{self.emapper_path}/KEGG_PATHWAY.txt"
        
        if os.path.exists(eggnog_path):
            print("eggnog-mapper already finished.")
            return
            
        cmd = f"{self.bash_env} {self.emapper_cmd} " \
              f"{prot_nonerude} {eggnog_profile} {eggnog_ann} " \
              f"{eggnog_KO} {eggnog_path}"
        try:
            self.run_in_conda_env(cmd)
            print(f"eggnog-mapper is successfully.")
        except RuntimeError as e:
            print(f"Error eggnog-mapper: {str(e)}")           

    def gene_count_bwa(self,args):
        """bwa
        """
        sam_name, sam_info = args
        sam_count_path = os.path.join(self.samcount_path, sam_name)
        
        if os.path.exists(sam_count_path):
            print(f"Gene count already finished: {sam_name}")
            return
            
        os.makedirs(sam_count_path, 0o755)
        geneset_bwa = f"{self.cdhit_path}/index_geneset"
        sam_r1_unmap = f"{self.ref_remove_path}/{sam_name}/{sam_name}.unmap.1.fastq.gz"
        sam_r2_unmap = f"{self.ref_remove_path}/{sam_name}/{sam_name}.unmap.2.fastq.gz"
        sam_count_bam = f"{sam_count_path}/{sam_name}_mapping_geneset.bam"
        sam_count_txt = f"{sam_count_path}/{sam_name}.count"
        
        cmd = f"{self.bash_env} {self.count_bwa_cmd} " \
              f"{geneset_bwa} {sam_r1_unmap} {sam_r2_unmap} " \
              f"{sam_count_bam} {sam_name} {sam_count_txt}"
        
        if os.path.exists(sam_count_txt):
            print("bwa mapping already finished.")
            return      
        try:
            self.run_command(cmd)
            print(f"bwa mapping is successfully.")
        except RuntimeError as e:
            print(f"Error bwa mapping: {str(e)}")  

    def gene_count_bowtie2(self,args):
        """bowtie2
        """
        sam_name, sam_info = args
        sam_count_path = os.path.join(self.samcount_path, sam_name)
        
        if os.path.exists(sam_count_path):
            print(f"Gene count already finished: {sam_name}")
            return
            
        os.makedirs(sam_count_path, 0o755)
        geneset_bowtie2 = f"{self.cdhit_path}/index_geneset"
        sam_r1_unmap = f"{self.ref_remove_path}/{sam_name}/{sam_name}.unmap.1.fastq.gz"
        sam_r2_unmap = f"{self.ref_remove_path}/{sam_name}/{sam_name}.unmap.2.fastq.gz"
        sam_count_sam = f"{sam_count_path}/{sam_name}_mapping_geneset.sam"
        sam_count_bam = f"{sam_count_path}/{sam_name}_mapping_geneset.bam"
        sam_count_txt = f"{sam_count_path}/{sam_name}.count"
        
        cmd = f"{self.bash_env} {self.count_bowtie2_cmd} " \
              f"{geneset_bowtie2} {sam_r1_unmap} {sam_r2_unmap} {sam_count_sam} " \
              f"{sam_count_bam} {sam_name} {sam_count_txt}"
        if os.path.exists(sam_count_txt):
            print("bowtie2 mapping already finished.")
            return  
        try:
            self.run_command(cmd)
            print(f"bowtie2 mapping is successfully.")
        except RuntimeError as e:
            print(f"Error bowtie2 mapping: {str(e)}") 

    def gene_count_salmon(self,args):
        """salmon
        """
        sam_name, sam_info = args
        sam_count_path = os.path.join(self.samcount_path, sam_name)
        
        if os.path.exists(sam_count_path):
            print(f"Gene count already finished: {sam_name}")
            return
            
        os.makedirs(sam_count_path, 0o755)
        geneset_salmon = f"{self.cdhit_path}/index_geneset"
        sam_r1_unmap = f"{self.ref_remove_path}/{sam_name}/{sam_name}.unmap.1.fastq.gz"
        sam_r2_unmap = f"{self.ref_remove_path}/{sam_name}/{sam_name}.unmap.2.fastq.gz"
        salmon_out = f"{sam_count_path}/salmon_out"
        sam_quant = f"{sam_name}.quant"
        sam_count = f"{sam_count_path}/{sam_name}.count"
        sam_tpm = f"{sam_count_path}/{sam_name}.tpm"

        cmd = f"{self.bash_env} {self.count_salmon_cmd} " \
              f"{geneset_salmon} {sam_r1_unmap} {sam_r2_unmap} {salmon_out} " \
              f"{sam_quant} {sam_name} {sam_tpm} {sam_count}"
        if os.path.exists(sam_tpm):
            print("salmon mapping already finished.")
            return
        try:
            self.run_in_mamba_env(cmd,env_name="salmon")
            print(f"salmon mapping is successfully.")
        except RuntimeError as e:
            print(f"Error salmon mapping: {str(e)}") 

    def process_single_gene_count(self, args):
        """单个样本的基因丰度计算[暂不使用]
        """
        sam_name, sam_info = args
        sam_count_path = os.path.join(self.samcount_path, sam_name)
        
        if os.path.exists(sam_count_path):
            print(f"Gene count already finished: {sam_name}")
            return
            
        os.makedirs(sam_count_path, 0o755)
        geneset_bwa = f"{self.cdhit_path}/geneset_bwa"
        sam_r1_unmap = f"{self.ref_remove_path}/{sam_name}/{sam_name}.unmap.1.fastq.gz"
        sam_r2_unmap = f"{self.ref_remove_path}/{sam_name}/{sam_name}.unmap.2.fastq.gz"
        sam_count_bam = f"{sam_count_path}/{sam_name}_mapping_geneset.bam"
        sam_count_txt = f"{sam_count_path}/{sam_name}.count"
        
        cmd = f"{self.bash_env} {self.count_cmd} " \
              f"{geneset_bwa} {sam_r1_unmap} {sam_r2_unmap} " \
              f"{sam_count_bam} {sam_name} {sam_count_txt}"
        self.run_in_conda_env(cmd)
        # os.system(command=cmd)

    def sam_gene_count(self, args):
        """并行运行基因丰度计算，并合并所有样本结果
        """
        self.samcount_path = os.path.join(args.outDir, "08-sam_count")
        if not os.path.exists(self.samcount_path):
            os.makedirs(self.samcount_path, 0o755)
            
        # 并行处理各个样本的基因计数
        if self.map_method == "bwa":
            with Pool(self.num_processes) as pool:
                pool.map(self.gene_count_bwa, self.raw_dict.items())
        if self.map_method == "bowtie2":
            with Pool(self.num_processes) as pool:
                pool.map(self.gene_count_bowtie2, self.raw_dict.items())
        if self.map_method == "salmon":
            with Pool(self.num_processes) as pool:
                pool.map(self.gene_count_salmon, self.raw_dict.items())
        # with Pool(self.num_processes) as pool:
        #     pool.map(self.process_single_gene_count, self.raw_dict.items())
        
        # 合并所有样本的结果[reads数目]
        merge_file = f"{self.samcount_path}/merged_file.txt"
        if os.path.exists(merge_file):
            print(f"Merge files already finished, check file: {merge_file}")
            return
            
        # 收集所有count文件
        count_files = [f"{self.samcount_path}/{sam_name}/{sam_name}.count" 
                      for sam_name in self.raw_dict.keys()]
        
        # 分批处理大量样本
        chunks = [count_files[i:i + 50] for i in range(0, len(count_files), 50)]
        split_files = []
        
        # 处理每个批次
        for i, chunk in enumerate(chunks):
            df_list = [pd.read_csv(f, sep='\t') for f in chunk]
            df_merged = pd.concat(df_list).groupby('gene', as_index=False).sum()
            df_merged.fillna(0, inplace=True)
            
            split_file = f"{self.samcount_path}/splitFile_{i}.txt"
            df_merged.to_csv(split_file, sep='\t', index=False)
            split_files.append(split_file)
        
        # 最终合并
        if len(split_files) == 1:
            os.rename(split_files[0], merge_file)
        else:
            df_merged = pd.read_csv(split_files[0], sep='\t')
            for f in split_files[1:]:
                df = pd.read_csv(f, sep='\t')
                df_merged = pd.merge(df_merged, df, on='gene', how='outer')
            df_merged.fillna(0, inplace=True)
            df_merged.to_csv(merge_file, sep='\t', index=False)
            
            # 清理临时文件
            for f in split_files:
                os.remove(f)

        # salmon再合并一个tpm的结果
        if self.map_method == "salmon":
            # 合并所有样本的结果[tpm]
            merge_file = f"{self.samcount_path}/merged_file_tpm.txt"
            if os.path.exists(merge_file):
                print(f"Merge files already finished, check file: {merge_file}")
                return
                
            # 收集所有count文件
            count_files = [f"{self.samcount_path}/{sam_name}/{sam_name}.tpm" 
                        for sam_name in self.raw_dict.keys()]
            
            # 分批处理大量样本
            chunks = [count_files[i:i + 50] for i in range(0, len(count_files), 50)]
            split_files = []
            
            # 处理每个批次
            for i, chunk in enumerate(chunks):
                df_list = [pd.read_csv(f, sep='\t') for f in chunk]
                df_merged = pd.concat(df_list).groupby('gene', as_index=False).sum()
                df_merged.fillna(0, inplace=True)
                
                split_file = f"{self.samcount_path}/splitFile_{i}.txt"
                df_merged.to_csv(split_file, sep='\t', index=False)
                split_files.append(split_file)
            
            # 最终合并
            if len(split_files) == 1:
                os.rename(split_files[0], merge_file)
            else:
                df_merged = pd.read_csv(split_files[0], sep='\t')
                for f in split_files[1:]:
                    df = pd.read_csv(f, sep='\t')
                    df_merged = pd.merge(df_merged, df, on='gene', how='outer')
                df_merged.fillna(0, inplace=True)
                df_merged.to_csv(merge_file, sep='\t', index=False)
                
                # 清理临时文件
                for f in split_files:
                    os.remove(f)

    def eggnog_kegg(self, args):
        """eggnog-kegg处理
        """
        self.kegg_path = os.path.join(args.outDir, "09-emapper_kegg")
        if not os.path.exists(self.kegg_path):
            os.makedirs(self.kegg_path, 0o755)
        eggnog_KO = f"{self.emapper_path}/KEGG_KO.txt"
        eggnog_path = f"{self.emapper_path}/KEGG_PATHWAY.txt"
        merge_file = f"{self.samcount_path}/merged_file.txt"
        out_ko = f"{self.kegg_path}/KO_samples.txt"
        out_pathway = f"{self.kegg_path}/pathway_samples.txt"
        tmp_dir = f"{self.kegg_path}/tmp"
        
        if os.path.exists(out_pathway):
            print("eggnog split kegg already finished.")
            return
        
        cmd = f"{self.py_env} {self.kegg_cmd} " \
              f"-kk {eggnog_KO} -kp {eggnog_path} " \
              f"-mt {merge_file} -ok {out_ko} " \
              f"-op {out_pathway} -t {tmp_dir}"
        try:
            self.run_in_conda_env(cmd)
            print(f"eggnog kegg split is successfully.")
        except RuntimeError as e:
            print(f"Error eggnog kegg split: {str(e)}")         
        
        # salmon再合并一个tpm的结果
        if self.map_method == "salmon":
            merge_file = f"{self.samcount_path}/merged_file_tpm.txt"
            out_ko = f"{self.kegg_path}/KO_samples_tpm.txt"
            out_pathway = f"{self.kegg_path}/pathway_samples_tpm.txt"
            tmp_dir = f"{self.kegg_path}/tmp_tpm"
            
            if os.path.exists(out_pathway):
                print("eggnog split kegg already finished.")
                return
            
            cmd = f"{self.py_env} {self.kegg_cmd} " \
                f"-kk {eggnog_KO} -kp {eggnog_path} " \
                f"-mt {merge_file} -ok {out_ko} " \
                f"-op {out_pathway} -t {tmp_dir}"
            try:
                self.run_in_conda_env(cmd)
                print(f"eggnog kegg split is successfully.")
            except RuntimeError as e:
                print(f"Error eggnog kegg split: {str(e)}") 

    def cazy_anno(self, args, evalue=0.00001, max_target=1):
        """CAZY功能注释[暂不运行]
        """
        self.cazy_path = os.path.join(args.outDir, "10-cazy")
        if not os.path.exists(self.cazy_path):
            os.makedirs(self.cazy_path, 0o755)

        prot_nonerude = f"{self.cdhit_path}/prot_nonerude.faa"
        cazy_m6 = f"{self.cazy_path}/cazy_m6.txt"
        merge_file = f"{self.samcount_path}/merged_file.txt"
        cazy_out = f"{self.cazy_path}/cazy_out.txt"
        
        if os.path.exists(cazy_out):
            print("CAZY annotation already finished.")
            return
        
        # Diamond blast against CAZY database
        cmd_1 = f"{self.diamond_cmd} {self.cazy_diamond_index} " \
                f"{prot_nonerude} {cazy_m6} {evalue} {max_target}"
        self.run_in_conda_env(cmd_1)
        # os.system(command=cmd_1)
        
        # CAZY annotation
        cmd_2 = f"{self.py_env} {self.cazy_cmd} " \
                f"-d {cazy_m6} -m {merge_file} -o {cazy_out}"
        self.run_in_conda_env(cmd_2)
        # os.system(command=cmd_2)

    def card_anno(self, args):
        """CARD功能注释
        """
        self.card_path = os.path.join(args.outDir, "10-card")
        if not os.path.exists(self.card_path):
            os.makedirs(self.card_path, 0o755)
        
        # 设置输入输出文件路径
        prot_nonerude = f"{self.cdhit_path}/prot_nonerude.faa"
        card_raw = f"{self.card_path}/card"
        card_info = f"{self.card_path}/card.txt"
        card_merge = f"{self.card_path}/card_merge.txt"
        merge_file = f"{self.samcount_path}/merged_file.txt"
                
        if os.path.exists(card_merge):
            print("CARD annotation already finished.")
            return
        
        # RGI分析
        cmd1 = f"{self.bash_env} {self.rgi_cmd} {prot_nonerude} {card_raw}"
        try:
            self.run_in_mamba_env(cmd1)
            print(f"RGI card is successfully.")
        except RuntimeError as e:
            print(f"Error rgi card: {str(e)}")    
        
        # CARD结果整理
        cmd2 = f"{self.py_env} {self.card_cmd} " \
               f"--input {card_info} " \
               f"--merge_table {merge_file} " \
               f"--output {self.card_path}"
        try:
            self.run_command(cmd2)
            print(f"card split is successfully.")
        except RuntimeError as e:
            print(f"Error card split: {str(e)}")
        
        # salmon再输出一个整理tpm的结果
        if self.map_method == "salmon":
            card_tpm = os.path.join(self.card_path,"card_tpm")
            if not os.path.exists(card_tpm):
                os.makedirs(card_tpm, 0o755)
            merge_file = f"{self.samcount_path}/merged_file_tpm.txt"
            
            cmd2 = f"{self.py_env} {self.card_cmd} " \
                f"--input {card_info} " \
                f"--merge_table {merge_file} " \
                f"--output {card_tpm}"
            try:
                self.run_command(cmd2)
                print(f"card split is successfully.")
            except RuntimeError as e:
                print(f"Error card split: {str(e)}")                       

    def result_merge(self, args):
        """整理并复制分析结果到结果目录
        """
        # 创建结果目录
        self.result_path = Path(args.outDir) / "00-result"
        self.result_path.mkdir(parents=True, exist_ok=True, mode=0o755)
        
        # 定义需要复制的文件映射
        files_to_copy = {
            # 源文件路径: (目标文件名, 描述)
            Path(self.metaphlan_path).glob("*.txt"): 
                ("", "物种组成结果"),
            Path(self.cdhit_path) / "prot_nonerude.faa": 
                ("prot_nonerude.faa", "非冗余蛋白序列"),
            Path(self.cdhit_path) / "nucl_nonerude.fna": 
                ("nucl_nonerude.fna", "非冗余核酸序列"), 
            Path(self.cdhit_path) / "geneset_length.txt": 
                ("geneset_length.txt", "基因长度信息"), 
            Path(self.emapper_path) / "eggnog.emapper.annotations":
                ("eggnog.emapper.annotations", "eggNOG注释结果"),
            Path(self.kegg_path).glob("*.txt"):
                ("", "KEGG通路注释结果"),
            Path(self.samcount_path) / "merged_file.txt":
                ("merged_file.txt", "基因丰度矩阵"),
            Path(self.samcount_path) / "merged_file_tpm.txt":
                ("merged_file_tpm.txt", "基因丰度-TPM矩阵"),  
            Path(self.card_path) / "card_merge.txt":
                ("card_merge.txt", "抗生素抗性基因注释"),
            Path(self.card_path) / "card_tpm":
                ("card_tpm", "抗生素抗性基因注释TPM结果")
        }

        # 复制文件
        for src, (dst_name, desc) in files_to_copy.items():
            try:
                if isinstance(src, Path):
                    # 单个文件
                    if not src.exists():
                        print(f"跳过不存在的文件: {src}")
                        continue
                        
                    dst = self.result_path / (dst_name or src.name)
                    if dst.exists():
                        print(f"目标文件已存在，跳过: {dst}")
                        continue
                        
                    if src.is_dir():  # 检查是否为文件夹
                        shutil.copytree(src, dst)  # 复制文件夹
                        print(f"已复制 {desc}: {src.name}")
                    else:
                        shutil.copy2(src, dst)  # 复制文件
                        print(f"已复制 {desc}: {src.name}")
                else:
                    # 通配符匹配的多个文件
                    files = list(src)  # 转换为列表以检查是否为空
                    if not files:
                        print(f"未找到匹配的文件: {src}")
                        continue
                        
                    for f in files:
                        dst = self.result_path / f.name
                        if dst.exists():
                            print(f"目标文件已存在，跳过: {dst}")
                            continue
                            
                        shutil.copy2(f, dst)
                        print(f"已复制 {desc}: {f.name}")
                        
            except Exception as e:
                print(f"复制文件失败 {src}: {str(e)}")

    def main(self):
        """主流程"""
        self.load_sample_info(args)
        self.fastp_trim(args)
        self.ref_remove(args)
        self.metaphlan(args)
        self.megahit(args)
        self.prodigal(args)
        self.cdhit(args)
        self.emapper(args)
        self.sam_gene_count(args)
        self.eggnog_kegg(args)
        # self.cazy_anno(args)
        self.card_anno(args)
        self.result_merge(args)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="PipeMetagenome",
                                     usage="\n"
                                     "python pipe_metagenome_parallel.py --fastq_list fq.list --output_dir result --ref ref_bowtie2_index --map_method salmon --processes 32\n\n"
                                     "ref_bowtie2_index:\n"
                                     "canis: /root/database/Canis_GCF_000002285.5/Canis_GCF_000002285_5\n"
                                     "canis_lupus_dingo: /root/database/Canis_lupus_dingo_GCF_003254725.2/Canis_lupus_dingo\n"
                                     "human: /root/database/hg38_GCF_000001405.40/GCF_000001405.40/hg38\n",
                                     description="Parallel pipeline of metagenome")
    parser.add_argument("-l", "--fastq_list", dest="fqList", required=True, type=str, 
                        help="[required]Input samples raw fq list")
    parser.add_argument("-o", "--output_dir", dest="outDir", required=False, type=str, default="MetagenomeResults",
                        help="[optional]Output results dir [default: MetagenomeResults]")
    parser.add_argument("-r", "--ref", dest="ref", required=False, type=str, default="/root/database/hg38_GCF_000001405.40/GCF_000001405.40/hg38",
                        help="[optional]Required reference genome bowtie2 index [default: hg38]")
    parser.add_argument("-r1", "--ref1", dest="ref1", required=False, type=str, 
                        help="[optional]Secondary reference genome bowtie2 index")
    parser.add_argument("-m", "--map_method", dest="mapMethod", required=True, type=str, 
                        help="[required]Gene conut mapping method [one of bwa, bowtie2, salmon]")   
    parser.add_argument("-p", "--processes", dest="processes", required=False, type=int, default=32, 
                        help="[optional]Number of processes to use [default: 32]")
    
    args = parser.parse_args()

    if not os.path.exists(args.outDir):
        os.makedirs(args.outDir)
    run = PipeMetagenome()
    if args.processes:
        run.num_processes = args.processes
    if args.mapMethod not in ['bwa', 'bowtie2', 'salmon']:
        exit("Please input right mapping method.")
    run.main()
    
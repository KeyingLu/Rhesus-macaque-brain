# -*- coding:utf-8 -*-

import os
from multiprocessing import Pool


blast = "/mnt/xdlab1/tools/ncbi-blast-2.8.0+/bin"
gffread = "/mnt/xdlab1/tools/gffread-0.9.12.Linux_x86_64/gffread"
faSomeRecords = "/mnt/xdlab1/tools/faSomeRecords"

human_dir = "/mnt/xdlab/ref/GRCh38_ensembl"
human_fasta = human_dir + "/Homo_sapiens.GRCh38.dna_sm.fa"
human_gtf = human_dir + "/Homo_sapiens.GRCh38.91.chr.gtf"


work_dir = "/mnt/xdlab1/home/looking/brain_project"
LNCipedia_fa = work_dir + "/ANNOTATION/lncipedia_5_2_hc.fasta"
GENCODE = work_dir + "/ANNOTATION/gencode.v28.long_noncoding_RNAs_ex_chr.gtf"
NONCODE = work_dir + "/ANNOTATION/NONCODEv5_human_hg38_lncRNA_ex_chr.gtf"
lncRNA_trans = work_dir + "/matrix0420/de_novo_TACO_minExpr5.5_res/lncRNA_transcripts_id_2845.txt"
lncRNA_trans_fa = lncRNA_trans + ".fa"

#sed 's/^chr//g' NONCODEv5_human_hg38_lncRNA.gtf > NONCODEv5_human_hg38_lncRNA_ex_chr.gtf


############################################################
## GENCODE
GENCODE_fa = GENCODE + ".fa"
os.system("%s -w %s -g %s %s" %
          (gffread, GENCODE_fa, human_fasta, GENCODE))

GENCODE_fa_origin = work_dir + "/ANNOTATION/gencode.v28.lncRNA_transcripts.fa"
GENCODE_fa = work_dir + "/ANNOTATION/gencode.v28.lncRNA_transcripts_rename.fa"
transcript_id = ""
with open(GENCODE_fa, "w") as f_out:
    with open(GENCODE_fa_origin, "r") as f_in:
        for line in f_in:
            info = line
            if line[0] == ">":
                transcript_id = line.split("|")[0]
                info = transcript_id + "\n"
            f_out.write(info)

GENCODE_fa_blastdb = GENCODE_fa + ".blastdb"
os.system('%s/makeblastdb -in %s -dbtype nucl -parse_seqids '
          '-out %s' %
          (blast, GENCODE_fa, GENCODE_fa_blastdb))

lncRNA_trans_GENCODE = work_dir + "/matrix0420/de_novo_TACO_minExpr5.5_res/lncRNA_GENCODE_blastn_res.txt"
os.system('%s/blastn -num_threads 10 -query %s -db %s '
          '-out %s -outfmt 6 ' %
          (blast, lncRNA_trans_fa, GENCODE_fa_blastdb, lncRNA_trans_GENCODE))

## NONCODE

# NONCODE_fa = NONCODE + ".fa"
# os.system("%s -w %s -g %s %s" %
#           (gffread, NONCODE_fa, human_fasta, NONCODE))


NONCODE_fa = work_dir + "/ANNOTATION/NONCODEv5_human.fa"
NONCODE_fa_blastdb = NONCODE_fa + ".blastdb"
os.system('%s/makeblastdb -in %s -dbtype nucl -parse_seqids '
          '-out %s' %
          (blast, NONCODE_fa, NONCODE_fa_blastdb))

lncRNA_trans_NONCODE = work_dir + "/matrix0420/de_novo_TACO_minExpr5.5_res/lncRNA_NONCODE_blastn_res.txt"
os.system('%s/blastn -num_threads 10 -query %s -db %s '
          '-out %s -outfmt 6 ' %
          (blast, lncRNA_trans_fa, NONCODE_fa_blastdb, lncRNA_trans_NONCODE))


## LNCipedia
LNCipedia_fa_blastdb = LNCipedia_fa + ".blastdb"
os.system('%s/makeblastdb -in %s -dbtype nucl -parse_seqids '
          '-out %s' %
          (blast, LNCipedia_fa, LNCipedia_fa_blastdb))

lncRNA_trans_LNCipedia = work_dir + "/matrix0420/de_novo_TACO_minExpr5.5_res/lncRNA_LNCipedia_blastn_res.txt"
os.system('%s/blastn -num_threads 10 -query %s -db %s '
          '-out %s -outfmt 6 ' %
          (blast, lncRNA_trans_fa, LNCipedia_fa_blastdb, lncRNA_trans_LNCipedia))



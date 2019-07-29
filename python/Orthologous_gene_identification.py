# -*- coding:utf-8 -*-

import os
from multiprocessing import Pool


TransDecoder = "/mnt/data1/Tools/TransDecoder-TransDecoder-v5.2.0"
samtools = "/mnt/data1/Tools/samtools1.35/samtools"
blast = "/mnt/data1/Tools/ncbi-blast-2.8.0+/bin"
hmmer = "/mnt/data1/Tools/hmmer-3.1b2-linux-intel-x86_64/binaries"
Inparanoid = "/mnt/data1/Tools/InParanoid/tmp/tmpuBl5dw/inparanoid_4.1"
faSomeRecords = "/mnt/data1/Tools/faSomeRecords"

Rhesus_dir =  "/mnt/data1/Ref/Rhesus_macaque_8.0.1"
Rhesus_fasta = Rhesus_dir + "/Macaca_mulatta.Mmul_8.0.1.dna_sm.fa"
gtf = Rhesus_dir + "/Macaca_mulatta.Mmul_8.0.1.91.chr.gtf"
ref = Rhesus_dir + "/Macaca_mulatta.Mmul_8.0.1.dna_sm_chr_tran"


work_dir = "/mnt/data2/Rhesus_brain"



###################################################################
#               extract the long open reading frames              #
###################################################################
TACO_minExpr_gtf = work_dir + "/TACO_minExpr_5.5/assembly.refcomp.gtf"
coding_transcripts_fa = TACO_minExpr_gtf + "_CPC2_CPAT_coding.fa"
os.chdir(work_dir + "/TACO_minExpr_5.5/")
os.system('%s/TransDecoder.LongOrfs -t %s'%
          (TransDecoder, coding_transcripts_fa))

# Use file: /mnt/data2/Rhesus_brain/TACO_minExpr_5.5
# /assembly.refcomp.gtf_CPC2_CPAT_coding.fa.transdecoder_dir
# /longest_orfs.pep
# for Pfam and/or BlastP searches to enable homology-based coding region identification.


############################################
#          BlastP Search with human        #
############################################
swissprot_fa = work_dir + "/ANNOTATION/Uniprot_Reviewed_Human.fasta"
swissprot = work_dir + "/ANNOTATION/Uniprot_Reviewed_Human.fasta.blastdb"
# os.system('%s/makeblastdb -in %s -dbtype prot -parse_seqids '
#           '-out %s' %
#           (blast, swissprot_fa, swissprot))
transdecoder_dir = coding_transcripts_fa + ".transdecoder_dir"
blastp_out = transdecoder_dir + "/blastp.outfmt6"
os.system('%s/blastp -query %s/longest_orfs.pep '
          '-db %s  -max_target_seqs 1 '
          '-outfmt 6 -evalue 1e-5 -num_threads 10 > %s' %
          (blast, transdecoder_dir, swissprot, blastp_out))
# finish


############################################
#          BlastP Search with Mouse        #
############################################
# swissprot_fa = work_dir + "/ANNOTATION/Uniprot_Reviewed_Mouse.fasta"
# swissprot = work_dir + "/ANNOTATION/Uniprot_Reviewed_Mouse.fasta.blastdb"
# os.system('%s/makeblastdb -in %s -dbtype prot -parse_seqids '
#           '-out %s' %
#           (blast, swissprot_fa, swissprot))
# transdecoder_dir = coding_transcripts_fa + ".transdecoder_dir"
# blastp_out = transdecoder_dir + "/blastp_Mouse.outfmt6"
# os.system('%s/blastp -query %s/longest_orfs.pep '
#           '-db %s  -max_target_seqs 1 '
#           '-outfmt 6 -evalue 1e-5 -num_threads 10 > %s' %
#           (blast, transdecoder_dir, swissprot, blastp_out))
# # finish


############################################
#         Pfam Search with human           #
############################################
transdecoder_dir = coding_transcripts_fa + ".transdecoder_dir"
Pfam_db = work_dir + "/ANNOTATION/Pfam-A.hmm"
# os.system('%s/hmmpress %s' % (hmmer, Pfam_db))
os.system('%s/hmmscan --cpu 15 '
          '--domtblout %s/pfam.domtblout '
          '%s '
          '%s/longest_orfs.pep' %
          (hmmer, transdecoder_dir, Pfam_db, transdecoder_dir))


#####################################################
#   Integrating the Blast and Pfam search results   #
#          into coding region selection             #
#####################################################
os.chdir(work_dir + "/gtf_for_assembled/TACO_minExpr_5.5/")
os.system('%s/TransDecoder.Predict '
          '-t %s '
          '--retain_pfam_hits %s/pfam.domtblout '
          '--retain_blastp_hits %s/blastp.outfmt6' %
          (TransDecoder, coding_transcripts_fa, transdecoder_dir, transdecoder_dir))



############################################
#                Inparanoid                #
############################################
# 文件都在同一目录下
Rhesus_transdecoder_pep = coding_transcripts_fa + ".transdecoder.pep"
swissprot_fa = work_dir + "/ANNOTATION/Uniprot_Reviewed_Human.fasta"
os.system('cp %s %s' % (swissprot_fa, os.path.dirname(Rhesus_transdecoder_pep)))
os.system('%s/inparanoid.pl %s/Uniprot_Reviewed_Human.fasta %s ' %
          (Inparanoid, os.path.dirname(Rhesus_transdecoder_pep), Rhesus_transdecoder_pep))



############################################
#            remove description            #
############################################
pep = "/home/looking/test/assembly.refcomp.gtf_CPC2_CPAT_coding.fa.transdecoder.pep"
pep2 = "/home/looking/test/assembly.refcomp.gtf_CPC2_CPAT_coding.fa.transdecoder_without_description.pep"
f1 = open(pep, "r")
line = f1.readline()
f2 = open(pep2, "w")
while line:
    if line[0] == ">":
        line = line.split()[0]
    f2.write(line)
    f2.flush()
    line = f1.readline()


f1.close()
f2.close()


############################################
#    human protein to human transcripts    #
############################################

##  nucleotide
human_transcripts_fa = human_gtf + ".fa"
human_protein_transcript_id = work_dir + "/gtf_for_assembled/TACO_minExpr_5.5/" \
                                         "CPC2_CPAT_inparanoid/human_protein_transcript_id.txt"
human_protein_transcript_fa = work_dir + "/gtf_for_assembled/TACO_minExpr_5.5/" \
                                         "CPC2_CPAT_inparanoid/human_protein_transcript_id.fa"
os.system('%s %s %s %s ' %
          (faSomeRecords, human_transcripts_fa,
           human_protein_transcript_id, human_protein_transcript_fa))

##### protein to nucleotide
swissprot_fa = work_dir + "/ANNOTATION/Uniprot_Reviewed_Human.fasta"

##### tblastn socre
blast = "/mnt/xdlab1/tools/ncbi-blast-2.8.0+/bin"
human_dir = "/mnt/xdlab/ref/GRCh38_ensembl"
human_gtf = human_dir + "/Homo_sapiens.GRCh38.91.chr.gtf"
human_transcripts_blastdb = human_gtf + ".fa.blastdb"
out = work_dir + "/gtf_for_assembled/TACO_minExpr_5.5/CPC2_CPAT_inparanoid/" \
                 "human_protein_transcript_fa_tblastn_res.txt"
os.system('%s/tblastn -num_threads 10 -query %s -db %s '
          '-out %s -outfmt 6' %
          (blast, swissprot_fa,
           human_transcripts_blastdb, out))



























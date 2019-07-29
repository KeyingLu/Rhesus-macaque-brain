# -*- coding:utf-8 -*-


import os
from multiprocessing import Pool

# read_distribution_tool = "/mnt/data1/Tools//RSeQC-2.6.4/scripts/read_distribution.py"
# CollectRnaSeqMetrics = "/mnt/data1/Tools//picard-tools2-1.119/CollectRnaSeqMetrics.jar"
# gtfToGenePred = "/mnt/data1/Tools//gtfToGenePred"
# gtf2bed = "/home/looking/software/gtf2bed"
qualimap = "/mnt/data1/Tools//qualimap_v2.2.1/qualimap"
samtools = "/mnt/data1/Tools//samtools1.35/samtools"
Rhesus_dir = "/mnt/data1/Ref/Rhesus_macaque_8.0.1"
gtf = Rhesus_dir + "/Macaca_mulatta.Mmul_8.0.1.91.chr.gtf"



work_dir = "/mnt/data2/Rhesus_brain"
bam_dir = work_dir + "/RNA_bam"
anno_dir = work_dir + "/ANNOTATION"


def load_single_line(filenames):
    f = open(filenames, 'r')
    sample_list = []
    line = f.readline()
    while line:
        sample_list.append(line.strip())
        line = f.readline()
    f.close()
    return sample_list


########################################
def reads_distribution(sample_id, out_dir, region_file):
    bam_sort = bam_dir + "/%s_sort.bam" % sample_id
    os.system('%s rnaseq --java-mem-size=4G -bam %s -gtf %s -outdir %s -p strand-specific-forward' %
          (qualimap, bam_sort, region_file, out_dir))


sample_list = load_single_line(work_dir + "/sample_list/Rhesus_sample_list.txt")
pool = Pool(20)
for sample_id in sample_list[156:]:
    out_dir = work_dir + "/read_distribution_by_qualimap/%s" % sample_id
    if not os.path.exists(out_dir):
        os.system('mkdir -p %s' % out_dir)
    pool.apply_async(reads_distribution, args=(sample_id,out_dir, gtf))

pool.close()
pool.join()
del pool



##############################
def feature_info(feature, res_file):
    info = os.popen('grep "%s" %s' % (feature, res_file)).read().split()
    Reads = info[2].replace(',', '')
    Percent = info[3][1:-2]
    return [feature, Reads, Percent]

read_distribution_dir = work_dir + "/read_distribution_by_qualimap"
genomic_origin_file = read_distribution_dir  + "/genomic_origin.txt"
genomic_origin = open(genomic_origin_file, "w")
genomic_origin.write("sample_id\tFeature\tReads\tPercent\n")

sample_list = load_single_line(work_dir + "/sample_list/Rhesus_sample_list.txt")
for sample_id in sample_list:
    res_file = read_distribution_dir + "/%s/rnaseq_qc_results.txt" % sample_id
    genomic_origin.write(sample_id + "\t" + "\t".join(feature_info("exonic", res_file)) + "\n")
    genomic_origin.write(sample_id + "\t" + "\t".join(feature_info("intronic", res_file)) + "\n")
    genomic_origin.write(sample_id + "\t" + "\t".join(feature_info("intergenic", res_file)) + "\n")
    genomic_origin.flush()

genomic_origin.close()





#################################################################
#              the Fraction of reads aligned to rRNA            #
#################################################################
import os
Rhesus_dir = "/mnt/data1/Ref/Rhesus_macaque_8.0.1"
gtf = Rhesus_dir + "/Macaca_mulatta.Mmul_8.0.1.91.chr.gtf"
rRNA_gtf = Rhesus_dir + "/Macaca_mulatta.Mmul_8.0.1.91.rRNA.gtf"
work_dir = "/mnt/data2/Rhesus_brain"


cmd = 'grep "rRNA" %s > %s' % (gtf, rRNA_gtf)
os.system(cmd)


sample_list = load_single_line(work_dir + "/sample_list/Rhesus_sample_list.txt")
pool = Pool(20)
for sample_id in sample_list:
    out_dir = work_dir + "/read_rRNA_by_qualimap/%s" % sample_id
    if not os.path.exists(out_dir):
        os.system('mkdir -p %s' % out_dir)
    pool.apply_async(reads_distribution, args=(sample_id,out_dir, rRNA_gtf))

pool.close()
pool.join()
del pool


##############################
def feature_info(feature, res_file):
    info = os.popen('grep "%s" %s' % (feature, res_file)).read().split()
    Reads = info[2].replace(',', '')
    Percent = info[3][1:-2]
    return [feature, Reads, Percent]

read_distribution_dir = work_dir + "/read_rRNA_by_qualimap"
genomic_origin_file = read_distribution_dir  + "/reads_rRNA_fraction.txt"
genomic_origin = open(genomic_origin_file, "w")
genomic_origin.write("sample_id\tFeature\tReads\tPercent\n")

sample_list = load_single_line(work_dir + "/sample_list/Rhesus_sample_list.txt")
for sample_id in sample_list:
    res_file = read_distribution_dir + "/%s/rnaseq_qc_results.txt" % sample_id
    genomic_origin.write(sample_id + "\t" + "\t".join(feature_info("exonic", res_file)) + "\n")
    genomic_origin.write(sample_id + "\t" + "\t".join(feature_info("intronic", res_file)) + "\n")
    genomic_origin.write(sample_id + "\t" + "\t".join(feature_info("intergenic", res_file)) + "\n")
    genomic_origin.flush()

genomic_origin.close()






#################################################################
#            the Fraction of reads aligned to novel gtf         #
#################################################################
import os
Rhesus_dir = "/mnt/data1/Ref/Rhesus_macaque_8.0.1"
gtf = Rhesus_dir + "/Macaca_mulatta.Mmul_8.0.1.91.chr.gtf"
novel_gtf = "/mnt/data2/Rhesus_brain/TACO_minExpr_5.5/assembly.gtf"
work_dir = "/mnt/data2/Rhesus_brain"


sample_list = load_single_line(work_dir + "/sample_list/Rhesus_sample_list.txt")
pool = Pool(20)
for sample_id in sample_list:
    out_dir = work_dir + "/read_novel_gtf_by_qualimap/%s" % sample_id
    if not os.path.exists(out_dir):
        os.system('mkdir -p %s' % out_dir)
    pool.apply_async(reads_distribution, args=(sample_id,out_dir, novel_gtf))

pool.close()
pool.join()
del pool



read_distribution_dir = work_dir + "/read_novel_gtf_by_qualimap"
genomic_origin_file = read_distribution_dir  + "/reads_novel_gtf_fraction.txt"
genomic_origin = open(genomic_origin_file, "w")
genomic_origin.write("sample_id\tFeature\tReads\tPercent\n")

sample_list = load_single_line(work_dir + "/sample_list/Rhesus_sample_list.txt")
for sample_id in sample_list:
    res_file = read_distribution_dir + "/%s/rnaseq_qc_results.txt" % sample_id
    genomic_origin.write(sample_id + "\t" + "\t".join(feature_info("exonic", res_file)) + "\n")
    genomic_origin.write(sample_id + "\t" + "\t".join(feature_info("intronic", res_file)) + "\n")
    genomic_origin.write(sample_id + "\t" + "\t".join(feature_info("intergenic", res_file)) + "\n")
    genomic_origin.flush()

genomic_origin.close()






















##########################
# gtf_bed = Rhesus_dir + "/Macaca_mulatta.Mmul_8.0.1.91.chr.gtf.bed"
# feature_stat = work_dir + "/4_6_feature_stat.txt"
# os.system('%s -i %s -r %s > %s' % (read_distribution_tool, bam, gtf_bed, feature_stat))

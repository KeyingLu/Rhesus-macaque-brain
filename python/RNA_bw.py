# -*- coding:utf-8 -*-

import os
import subprocess
from multiprocessing import Pool


bedtools = "/mnt/xdlab1/tools/bedtools2/bedtools2/bin/bedtools"
samtools = "/mnt/xdlab1/tools/samtools1.35/samtools"
bedGraphToBigWig = "/home/looking/software/bedGraphToBigWig"
fai_file = "/mnt/xdlab/ref/M.fascicularis/Mmul_8.0.1.91_ensembl/" \
           "Macaca_mulatta.Mmul_8.0.1.dna_sm.fa.fai"

work_dir = "/mnt/xdlab1/home/looking/brain_project"
bam_dir = os.path.join(work_dir, "bam_new")
fastq_dir = os.path.join(work_dir, "fastq")
bw_dir = os.path.join(work_dir, "bigwig")
bedgraph_dir = os.path.join(bw_dir, "bedgraph")

dir_list = [bedgraph_dir, bw_dir]


for dirs in dir_list:
    if not os.path.exists(dirs):
        os.system('mkdir {dirs}'.format(dirs=dirs))


def load_single_line(filenames):
    f = open(filenames, 'r')
    sample_list = []
    line = f.readline()
    while line:
        sample_list.append(line.strip())
        line = f.readline()
    f.close()
    return sample_list


def get_bigwig_with_scaled(sample_id):
    bam = os.path.join(bam_dir, "%s_sort.bam" % sample_id)
    bedgraph = os.path.join(bedgraph_dir, "%s.bdg" % sample_id)
    reads_count = int(os.popen("%s view -bf 0x2 %s |wc -l" %
                               (samtools, bam)).read().split(" ")[0])
    scale = 1000000.0 * 2 / reads_count
    os.system('{samtools} view -bf 0x2 {bam}|'
              '{bedtools} genomecov -ibam stdin -bg -split -scale {scale}> {bedgraph}'.format(
               samtools=samtools, bam=bam, bedtools=bedtools,
               scale=scale, bedgraph=bedgraph))
    bw = os.path.join(bw_dir, "%s.bw" % sample_id)
    os.system('sort -k1,1 -k2,2n %s -o %s' % (bedgraph, bedgraph))
    os.system('%s %s %s %s' % (bedGraphToBigWig, bedgraph, fai_file, bw))


# sample_list = load_single_line(os.path.join(fastq_dir, "sample_list"))
sample_list = ["11_24"]
pools = Pool(4)
for sample_id in sample_list:
    pools.apply_async(get_bigwig_with_scaled, args=(sample_id,))

pools.close()
pools.join()
del pools















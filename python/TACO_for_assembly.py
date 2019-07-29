# -*- coding:utf-8 -*-


import os
from multiprocessing import Pool



TACO = "/mnt/xdlab1/tools/taco-v0.7.3.Linux_x86_64/taco_run"
Refcomp = "/mnt/xdlab1/tools/taco-v0.7.3.Linux_x86_64/taco_refcomp"
gffcompare = "/mnt/xdlab1/tools/gffcompare-0.10.4.Linux_x86_64/gffcompare"
Rhesus_dir = "/mnt/xdlab/ref/M.fascicularis/Mmul_8.0.1.91_ensembl"
gtf = Rhesus_dir + "/Macaca_mulatta.Mmul_8.0.1.91.chr.gtf"


work_dir = "/mnt/xdlab1/home/looking/brain_project"
# gtf_dir = work_dir + "/training_stringTie_cov"

def load_single_line(filenames):
    f = open(filenames, 'r')
    sample_list = []
    line = f.readline()
    while line:
        sample_list.append(line.strip())
        line = f.readline()
    f.close()
    return sample_list

sample_list = load_single_line(work_dir + "/sample_list/Rhesus_sample_list.txt")


# for cov in xrange(5, 18):
#     cov = float(cov) / 2
#     cov_dir = gtf_dir + "/cov_%s" % cov
#     mergelist_file = cov_dir + "/mergelist.txt"
#     TACO_dir = cov_dir + "/TACO_res"
#     os.system('%s -p 8 -o %s %s' % (TACO, TACO_dir, mergelist_file))
#     TACO_gtf = TACO_dir  + "/assembly.gtf"
#     os.system("%s -p 8 -o %s -r %s -t %s" % (Refcomp, TACO_dir, gtf, TACO_gtf))
#     os.system('%s -o %s/gffcmp -r %s  %s' % (gffcompare, TACO_dir, gtf, TACO_gtf))
#

########################  minExpr
gtf_dir = work_dir + "/gtf_for_assembled"
for minExpr in [21,22,23,24,26,27,28,29]:
    minExpr = float(minExpr) / 10
    TACO_dir = gtf_dir + "/TACO_minExpr_%s" % minExpr
    mergelist_file = gtf_dir + "/mergelist.txt"
    os.system('%s -p 8 --filter-min-expr %s -o %s %s' % (TACO, minExpr, TACO_dir, mergelist_file))
    TACO_gtf = TACO_dir  + "/assembly.gtf"
    os.system("%s -p 8 -o %s -r %s -t %s" % (Refcomp, TACO_dir, gtf, TACO_gtf))
    os.system('%s -o %s/gffcmp -r %s  %s' % (gffcompare, TACO_dir, gtf, TACO_gtf))



# ########################  minExpr and gffcmp_QMN
# gtf_dir = work_dir + "/gtf_for_assembled"
# for minExpr in xrange(1, 17):
#     minExpr = float(minExpr) / 2
#     TACO_dir = gtf_dir + "/TACO_minExpr_%s" % minExpr
#     TACO_gtf = TACO_dir  + "/assembly.gtf"
#     os.system('%s -o %s/gffcmp_QMN -r %s -Q -M -N %s' % (gffcompare, TACO_dir, gtf, TACO_gtf))


#########################
def FRAC_test(gtf_dir, FRAC):
    TACO_dir = gtf_dir + "/TACO_FRAC_%s" % FRAC
    mergelist_file = gtf_dir + "/mergelist.txt"
    os.system('%s -p 8 --isoform-frac %s -o %s %s' % (TACO, FRAC, TACO_dir, mergelist_file))
    TACO_gtf = TACO_dir + "/assembly.gtf"
    os.system('%s -o %s/gffcmp -r %s  %s' % (gffcompare, TACO_dir, gtf, TACO_gtf))


gtf_dir = work_dir + "/gtf_for_assembled"
pool = Pool(5)
for FRAC in xrange(56, 100):
    FRAC = float(FRAC) / 100
    pool.apply_async(FRAC_test, args=(gtf_dir, FRAC, ))

pool.close()
pool.join()
del pool

##############################################################################
#                                 TACO final                                 #
##############################################################################
FRAC = 0.96
minExpr = 5.5
gtf_dir = work_dir + "/gtf_for_assembled"
TACO_dir = gtf_dir + "/TACO_FRAC%s_minExpr%s" % (FRAC, minExpr)
mergelist_file = gtf_dir + "/mergelist.txt"
os.system('%s -p 8 --isoform-frac %s --filter-min-expr %s -o %s %s' %
          (TACO, FRAC, minExpr, TACO_dir, mergelist_file))
TACO_gtf = TACO_dir + "/assembly.gtf"
os.system('%s -o %s/gffcmp -r %s  %s' % (gffcompare, TACO_dir, gtf, TACO_gtf))



# #########################
# gtf_dir = work_dir + "/gtf_for_assembled"
# for FRAC in xrange(1, 100):
#     FRAC = float(FRAC) / 100
#     TACO_dir = gtf_dir + "/TACO_FRAC_%s" % FRAC
#     TACO_gtf = TACO_dir  + "/assembly.gtf"
#     os.system('%s -o %s/gffcmp_QMN -r %s -Q -M -N %s' % (gffcompare, TACO_dir, gtf, TACO_gtf))







##############################################################################
#              statistic for TACO minExpr and gffcmp_QMN                     #
##############################################################################
work_dir = "/mnt/xdlab1/home/looking/brain_project"
gtf_dir = work_dir + "/gtf_for_assembled"
stat_for_plot_1 = gtf_dir + "/TACO_minExpr_gffcmp_QMN.stats_sensitivity_precision.txt"
stat_for_plot_2 = gtf_dir + "/TACO_minExpr_gffcmp_QMN.stats_Novel_Missed.txt"
stat1 = open(stat_for_plot_1, "w")
stat2 = open(stat_for_plot_2, "w")
stat1.write("minExpr\tType\tSensitivity_Precision\tPercent\n")
stat2.write("minExpr\tType\tNovel_or_Missed\tTotal\tPercent\n")

for minExpr in xrange(1, 17):
    minExpr = float(minExpr) / 2
    TACO_dir = gtf_dir + "/TACO_minExpr_%s" % minExpr
    gffcmp_stat_file = TACO_dir + "/gffcmp_QMN.stats"
    gffcmp_stat = open(gffcmp_stat_file, "r")
    line = gffcmp_stat.readline()
    while line:
        if line[0] == "#" or line == "\n":
            line = gffcmp_stat.readline()
            continue
        info = line.split()
        if info[1][:-1] == "level":
            Type = info[0] + "_" + info[1][:-1]
            Sensitivity = info[2]
            Precision = info[4]
            stat1.write("\t".join([str(minExpr), Type, "Sensitivity", Sensitivity]) + "\n")
            stat1.write("\t".join([str(minExpr), Type, "Precision", Precision]) + "\n")
            stat1.flush()
        if len(info[2].split("/")) == 2:
            Type = info[0] + "_" + info[1][:-1]
            Novel_or_Missed = info[2].split("/")[0]
            Total = info[2].split("/")[1]
            Percent = info[4][:-2]
            stat2.write("\t".join([str(minExpr), Type, Novel_or_Missed, Total, Percent]) + "\n")
            stat2.flush()
        line = gffcmp_stat.readline()
    gffcmp_stat.close()

stat1.close()
stat2.close()

##############################################################################
#                      statistic for TACO minExpr                            #
##############################################################################
work_dir = "/mnt/xdlab1/home/looking/brain_project"
gtf_dir = work_dir + "/gtf_for_assembled"
stat_for_plot_1 = gtf_dir + "/TACO_minExpr_gffcmp.stats_sensitivity_precision.txt"
stat_for_plot_2 = gtf_dir + "/TACO_minExpr_gffcmp.stats_Novel_Missed.txt"
stat_for_plot_3 = gtf_dir + "/TACO_minExpr_gffcmp.stats_transcripts_count.txt"
stat1 = open(stat_for_plot_1, "w")
stat2 = open(stat_for_plot_2, "w")
stat3 = open(stat_for_plot_3, "w")
stat1.write("minExpr\tType\tSensitivity_Precision\tPercent\n")
stat2.write("minExpr\tType\tNovel_or_Missed\tTotal\tPercent\n")
stat3.write("minExpr\ttranscripts_count\n")

# for minExpr in [0.5, 1.0, 1.5, 2, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6,
#                 2.7, 2.8, 2.9, 3.0,3.5, 4.0, 4.5, 5.0, 5.5, 6.0,
#                 6.5, 7.0, 7.5, 8, 8.5, 9, 9.5]:
for minExpr in [0.5, 1.0, 1.5, 2, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0,
                6.5, 7.0, 7.5, 8, 8.5, 9, 9.5]:
    minExpr = float(minExpr)
    TACO_dir = gtf_dir + "/TACO_minExpr_%s" % minExpr
    gffcmp_stat_file = TACO_dir + "/gffcmp.stats"
    gffcmp_stat = open(gffcmp_stat_file, "r")
    line = gffcmp_stat.readline()
    while line:
        if line[0] == "#" or line == "\n":
            line = gffcmp_stat.readline()
            continue
        info = line.split()
        if info[1][:-1] == "level":
            Type = info[0] + "_" + info[1][:-1]
            Sensitivity = info[2]
            Precision = info[4]
            F1 = 2*float(Sensitivity)*float(Precision)/(float(Sensitivity) + float(Precision))
            stat1.write("\t".join([str(minExpr), Type, "Sensitivity", Sensitivity]) + "\n")
            stat1.write("\t".join([str(minExpr), Type, "Precision", Precision]) + "\n")
            stat1.write("\t".join([str(minExpr), Type, "F1_score", str(F1)]) + "\n")
            stat1.flush()
        if len(info[2].split("/")) == 2:
            Type = info[0] + "_" + info[1][:-1]
            Novel_or_Missed = info[2].split("/")[0]
            Total = info[2].split("/")[1]
            Percent = info[4][:-2]
            stat2.write("\t".join([str(minExpr), Type, Novel_or_Missed, Total, Percent]) + "\n")
            stat2.flush()
        line = gffcmp_stat.readline()
    gffcmp_stat.close()
    gffcmp_stat = open(gffcmp_stat_file, "r")
    lines = gffcmp_stat.readlines()
    final_line = lines[len(lines)-1]
    trans_count = final_line.split()[0]
    stat3.write("\t".join([str(minExpr),trans_count]) + "\n")


stat1.close()
stat2.close()
stat3.close()


##############################################################################
#                        statistic for TACO FRAC                             #
##############################################################################
work_dir = "/mnt/xdlab1/home/looking/brain_project"
gtf_dir = work_dir + "/gtf_for_assembled"
stat_for_plot_1 = gtf_dir + "/TACO_FRAC_gffcmp.stats_sensitivity_precision.txt"
stat_for_plot_2 = gtf_dir + "/TACO_FRAC_gffcmp.stats_Novel_Missed.txt"
stat1 = open(stat_for_plot_1, "w")
stat2 = open(stat_for_plot_2, "w")
stat1.write("FRAC \tType\tSensitivity_Precision\tPercent\n")
stat2.write("FRAC \tType\tNovel_or_Missed\tTotal\tPercent\n")


for FRAC in xrange(1, 100):
    FRAC = float(FRAC) / 100
    TACO_dir = gtf_dir + "/TACO_FRAC_%s" % FRAC
    gffcmp_stat_file = TACO_dir + "/gffcmp.stats"
    gffcmp_stat = open(gffcmp_stat_file, "r")
    line = gffcmp_stat.readline()
    while line:
        if line[0] == "#" or line == "\n":
            line = gffcmp_stat.readline()
            continue
        info = line.split()
        if info[1][:-1] == "level":
            Type = info[0] + "_" + info[1][:-1]
            Sensitivity = info[2]
            Precision = info[4]
            F1 = 2*float(Sensitivity)*float(Precision)/(float(Sensitivity) + float(Precision))
            stat1.write("\t".join([str(FRAC), Type, "Sensitivity", Sensitivity]) + "\n")
            stat1.write("\t".join([str(FRAC), Type, "Precision", Precision]) + "\n")
            stat1.write("\t".join([str(FRAC), Type, "F1_score", str(F1)]) + "\n")
            stat1.flush()
        if len(info[2].split("/")) == 2:
            Type = info[0] + "_" + info[1][:-1]
            Novel_or_Missed = info[2].split("/")[0]
            Total = info[2].split("/")[1]
            Percent = info[4][:-2]
            stat2.write("\t".join([str(FRAC), Type, Novel_or_Missed, Total, Percent]) + "\n")
            stat2.flush()
        line = gffcmp_stat.readline()
    gffcmp_stat.close()

stat1.close()
stat2.close()



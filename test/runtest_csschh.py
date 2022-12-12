#!/usr/bin/env python

import subprocess as sbc
import sys
import argparse

############################################################

ni_values = ["5000000", "50000000", "500000000", "1000000000"]
ni_default = "50000000"

sk_values = ["1.0", "1.4", "1.8", "2.2"]
sk_default = "1.1"

phi1_values = ["0.001", "0.01", "0.1"]
phi1_default = "0.01"
phi2_values = ["0.01", "0.02", "0.04", "0.08", "0.1"]
phi2_default = "0.05"

k1_values = ["4200", "8400", "16800", "33600"]
k1_default = "8400"
k2_values = ["63000", "252000", "1008000", "4032000"]
k2_default = "756000"

pdi = "32768"
sdi = "32768"

seed_base = 16033099
seed_step = 50
step_num = 10

############################################################

header = "ni,Distinct_Prim,Distinct_Corr,seed,hseed,zipfpar,hwzpar,phi1,phi2," \
		 "eps1,eps2,TrueHH_Prim,TrueHH_Corr,Cand_Prim,Cand_Corr,MatchHH_Prim," \
		 "MatchHH_Corr,Recall_Prim,Precision_Prim,Recall_Corr,Precision_Corr," \
		 "meanAbsErrCan_Corr,maxAbsErrCan_Corr,meanRelErrCan_Corr,maxRelErrCan_Corr," \
		 "Measured_eps1,Measured_eps2,updatePerMs,queryTime(ms),space(Kb),k1,k2\n"

parser = argparse.ArgumentParser()

parser.add_argument("cmd", help="executable name")
parser.add_argument("test_file", help="name of file that stores results")

options = parser.parse_args()

exec_name = "./" + options.cmd


def print_to_stderr(msg):
    sys.stderr.write(msg)
    sys.stderr.flush()
    return


def test_on_size():
    outfile = open(options.cmd + "_" + options.test_file + "_on_size.csv", "w")
    outfile.write(header)
    outfile.flush()
    for i in range(4):

        print_to_stderr('Test k1 = ' + k1_values[i] + '\t k2 = ' + k2_values[i] + '\n')

        for seed in range(seed_base, seed_base + (step_num * seed_step), seed_step):
            sbc.call([exec_name,"-pdi", pdi, "-sdi", sdi, "-k1", k1_values[i], "-k2", k2_values[i], "-seed", str(seed), "-ni", ni_default,
            		 "-phi1", phi1_default, "-phi2", phi2_default, "-sk", sk_default], stderr=outfile, stdout=sbc.DEVNULL)
            print_to_stderr("#")

        print_to_stderr("\n")

        outfile.write("\n")
        outfile.flush()

    outfile.close()
    return


def test_on_thresh(phi1):
    outfile = open(options.cmd + "_" + options.test_file + "_on_thresh_phi" + phi1 + ".csv", "w")
    outfile.write(header)
    outfile.flush()
    for phi2 in phi2_values:

        print_to_stderr('Test phi2 = ' + phi2 + '\t\n')

        for seed in range(seed_base, seed_base + (step_num * seed_step), seed_step):
            sbc.call([exec_name,"-pdi", pdi, "-sdi", sdi,"-k1", k1_default, "-k2", k2_default, "-seed", str(seed), "-ni", ni_default,
            		 "-phi1", phi1, "-phi2", phi2, "-sk", sk_default], stderr=outfile, stdout=sbc.DEVNULL)
            print_to_stderr("#")

        print_to_stderr("\n")

        outfile.write("\n")
        outfile.flush()

    outfile.close()
    return


def test_on_skew():
    outfile = open(options.cmd + "_" + options.test_file + "_on_skew.csv", "w")
    outfile.write(header)
    outfile.flush()
    for sk in sk_values:

        print_to_stderr('Test sk = ' + sk + '\t\n')

        for seed in range(seed_base, seed_base + (step_num * seed_step), seed_step):
            sbc.call([exec_name,"-pdi", pdi, "-sdi", sdi,"-k1", k1_default, "-k2", k2_default, "-seed", str(seed), "-ni", ni_default,
            		 "-phi1", phi1_default, "-phi2", phi2_default, "-sk", sk], stderr=outfile, stdout=sbc.DEVNULL)
            print_to_stderr("#")

        print_to_stderr("\n")

        outfile.write("\n")
        outfile.flush()

    outfile.close()
    return


def test_on_ni():
    outfile = open(options.cmd + "_" + options.test_file + "_on_ni.csv", "w")
    outfile.write(header)
    outfile.flush()
    for ni in ni_values:

        print_to_stderr('Test ni = ' + ni + '\t\n')

        for seed in range(seed_base, seed_base + (step_num * seed_step), seed_step):
            sbc.call([exec_name,"-pdi", pdi, "-sdi", sdi,"-k1", k1_default, "-k2", k2_default, "-seed", str(seed), "-ni", ni,
            		 "-phi1", phi1_default, "-phi2", phi2_default, "-sk", sk_default], stderr=outfile, stdout=sbc.DEVNULL)
            print_to_stderr("#")

        print_to_stderr("\n")

        outfile.write("\n")
        outfile.flush()

    outfile.close()
    return


test_on_size()

test_on_ni()

for phi1 in phi1_values:
	test_on_thresh(phi1)

test_on_skew()

print_to_stderr("\nTest completed!!\n")

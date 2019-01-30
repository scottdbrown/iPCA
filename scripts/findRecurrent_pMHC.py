'''
Find Recurrent pMHC
Within TCGA data, look for recurrent pMHC that meet 500 nM threshold

Date: October 13, 2016
@author: sbrown

Edited March 16, 2017:
    - Take into account the expression data.
    - Also output a matrix of recurrent x subjects

Edited May 30, 2017:
    - Count ENSP-mut pairs for each pMHC.
'''

## Import Libraries
import sys
import argparse
import os

DEBUG = False
VERB = False


if __name__ == "__main__":

    ## Deal with command line arguments
    parser = argparse.ArgumentParser(description = "Find Recurrent pMHC")
    ## add_argument("name", "(names)", metavar="exampleOfValue - best for optional", type=int, nargs="+", choices=[allowed,values], dest="nameOfVariableInArgsToSaveAs")
    parser.add_argument("pMHC_file", help = "File with pMHC", type = str)
    parser.add_argument("pMHC_sample_file", help = "File with summary stats for pMHC to see set of samples analyzed", type = str)
    parser.add_argument("peptideOutputFile", help = "output file for binding peptides (ignore MHC)", type = str)
    parser.add_argument("peptideOutputMatrix", help = "output matrix for binding peptides (ignore MHC)", type = str)
    parser.add_argument("pmhcOutputFile", help = "output file for pMHCs", type = str)
    parser.add_argument("pmhcOutputMatrix", help = "output matrix for pMHCs", type = str)

    parser.add_argument("-d", "--debug", action = "store_true", dest = "DEBUG", help = "Flag for setting debug/test state.")
    parser.add_argument("-v", "--verbose", action = "store_true", dest = "VERB", help = "Flag for setting verbose output.")
    args = parser.parse_args()

    ## Set Global Vars
    DEBUG = args.DEBUG
    VERB = args.VERB


    epis = {}
    ## key will be epitope, with three sets: HLAs, GENEs, Subjects
    ## This will allow counts of why this epitope is frequent.

    pmhcs = {}
    ## key will be epitope-MHC pair, with two sets: GENEs, Subjects

    ## Get subjects analyzed.
    print("Reading in subjects analyzed...")
    subs = []
    HEADER = True
    for line in open(args.pMHC_sample_file, "r"):
        if HEADER:
            HEADER = False
        else:
            line = line.split("\t")
            if line[0][:12] not in subs:
                subs.append(line[0][:12])


    print("Reading in peptides...")
    HEADER = True
    for line in open(args.pMHC_file, "r"):
        if HEADER:
            HEADER = False
        else:
            barcode, ensp, mut, hla, pep, varpos, mutic, wtic, exp = line.rstrip().split("\t")
            mutic = float(mutic)
            if exp != "None":
                exp = float(exp)
            else:
                exp = None

            if mutic <= 500 and (exp is None or exp > 0.7):
                if pep not in epis:
                    epis[pep] = [set(),set(),set()]
                if (pep,hla) not in pmhcs:
                    pmhcs[(pep,hla)] = [set(),set()]

                epis[pep][0].add(hla)
                epis[pep][1].add("{}-{}".format(ensp,mut))
                epis[pep][2].add(barcode[:12])

                pmhcs[(pep,hla)][0].add("{}-{}".format(ensp,mut))
                pmhcs[(pep,hla)][1].add(barcode[:12])


    print("Writing and tallying recurrent binding peptides...")
    out = open(args.peptideOutputFile, "w")
    mat = open(args.peptideOutputMatrix, "w")
    out.write("peptide\tnumHLA\tnumMut\tnumSub\n")
    mat.write("peptide\t{}\n".format("\t".join(subs)))
    for pep in epis:
        out.write("{}\t{}\t{}\t{}\n".format(pep, len(epis[pep][0]), len(epis[pep][1]), len(epis[pep][2])))
        if len(epis[pep][2]) > 1:
            ## is recurrent
            epiline = [0 for i in range(len(subs))]
            for i in range(len(subs)):
                if subs[i] in epis[pep][2]:
                    ## subject has this pep
                    epiline[i] = 1
            mat.write("{}\t{}\n".format(pep, "\t".join(str(x) for x in epiline)))
    out.close()
    mat.close()

    print("Writing and tallying recurrent pMHCs...")
    out = open(args.pmhcOutputFile, "w")
    mat = open(args.pmhcOutputMatrix, "w")
    out.write("pMHC\tnumMut\tnumSub\n")
    mat.write("pMHC\t{}\n".format("\t".join(subs)))
    for pMHC in pmhcs:
        out.write("{}\t{}\t{}\n".format("_".join(pMHC), len(pmhcs[pMHC][0]), len(pmhcs[pMHC][1])))
        if len(pmhcs[pMHC][1]) > 1:
            ## is recurrent
            pmhcline = [0 for i in range(len(subs))]
            for i in range(len(subs)):
                if subs[i] in pmhcs[pMHC][1]:
                    ## subject has this pMHC
                    pmhcline[i] = 1
            mat.write("{}\t{}\n".format("_".join(pMHC), "\t".join(str(x) for x in pmhcline)))
    out.close()
    mat.close()

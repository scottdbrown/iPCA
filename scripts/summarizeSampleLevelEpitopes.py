'''
Summarize Sample Level Epitopes
Get number of nonsyn muts, binders, mut prots with at least one binder.

Date: September 1, 2016
@author: sbrown

Edited February 28, 2017:
    - Better sample-level stats queried.

Edited March 16, 2017:
    - Use expression threshold.
    - Keep samples that had no mutatons meeting filters (use HLA file as reference)
'''

## Import Libraries
import sys
import argparse

DEBUG = False
VERB = False

IC50_THRESH = 500
EXPRESSION_THRESH = 0.7

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

if __name__ == "__main__":

    ## Deal with command line arguments
    parser = argparse.ArgumentParser(description = "Summarize Sample Level Epitopes")
    ## add_argument("name", "(names)", metavar="exampleOfValue - best for optional", type=int, nargs="+", choices=[allowed,values], dest="nameOfVariableInArgsToSaveAs")
    parser.add_argument("epitopeListFile", help = "File of all epitopes", type = str)
    parser.add_argument("whitelist_file", help = "File with whitelist subjects", type = str)
    parser.add_argument("blacklist_file", help = "File with blacklist aliquots", type = str)
    parser.add_argument("outputFile", help = "File to write output", type = str)
    parser.add_argument("-d", "--debug", action = "store_true", dest = "DEBUG", help = "Flag for setting debug/test state.")
    parser.add_argument("-v", "--verbose", action = "store_true", dest = "VERB", help = "Flag for setting verbose output.")
    args = parser.parse_args()

    ## Set Global Vars
    DEBUG = args.DEBUG
    VERB = args.VERB

    whitelist_sub = set()   ## hold subjects that are okay to use.
    blacklist_ali = set()   ## hold aliquots that are not okay to use.

    HEADER = True
    for line in open(args.whitelist_file, "r"):
        if HEADER:
            HEADER = False
        else:
            whitelist_sub.add(line.rstrip().replace('"',''))

    for line in open(args.blacklist_file, "r"):
        blacklist_ali.add(line.rstrip())

    ## want the sample ID, number of mutations, number of binding epitopes, and number of genes with at least one binding epitopes.

    numMutations = {}
    numImmMutation = {}
    numPeptides = {}
    numBindingPMHC = {}
    numBindingExpressedPMHC = {}


    print("Processing epitope file...")
    HEADER = True
    for line in open(args.epitopeListFile, "r"):
        if HEADER:
            HEADER = False
        else:
            samp, ensp, mut, hla, mutPep, varPos, mutIC50, refIC50, exp = line.rstrip().split("\t")
            mutIC50 = float(mutIC50)
            if is_number(exp):
                exp = float(exp)
            else:
                exp = None



            if samp not in blacklist_ali and samp[:12] in whitelist_sub:
                ## if not, blacklisted, or not whitlisted
                if samp not in numMutations:
                    numMutations[samp] = set()
                    numImmMutation[samp] = set()
                    numPeptides[samp] = 0
                    numBindingPMHC[samp] = 0
                    numBindingExpressedPMHC[samp] = 0


                numMutations[samp].add("{}_{}".format(ensp, mut))
                numPeptides[samp] += 1
                if mutIC50 < IC50_THRESH:
                    numImmMutation[samp].add("{}_{}".format(ensp, mut))
                    numBindingPMHC[samp] += 1
                    if exp is None or exp > EXPRESSION_THRESH:
                        numBindingExpressedPMHC[samp] += 1
            else:
                print("{} is not being used.".format(samp))

    print("Writing...")
    out = open(args.outputFile, "w")
    out.write("barcode\tnumberOfNonSynonymousSNP\tnumberOfImmunogenicMutation\tnumberOfPeptideTested\tnumberOfBindingPMHC\tnumberOfBindingExpressedPMHC\n")
    for samp in numMutations:
        out.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(samp, len(numMutations[samp]), len(numImmMutation[samp]), numPeptides[samp], numBindingPMHC[samp], numBindingExpressedPMHC[samp]))
    out.close()

    print("done.")

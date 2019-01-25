'''
Parse NetMHCpan Output
Go through all directories and get output

Date: September 1, 2016
@author: sbrown

Edited February 28, 2017:
	- Make a bit more RAM efficient.
Edited March 16, 2017:
    - Allow multiple directories.
'''

## Import Libraries
import sys
import argparse
import os

DEBUG = False
VERB = False

resI = {}       ## hold results
pepInfo = {}    ## hold submitted 21(ish)mer info.

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

def parsePeptideInfo(filename):
    HEADER = True
    for line in open(filename, "r"):
        if HEADER:
            HEADER = False
        else:
            tumor, peptideID, ENSP, mutation, varPos, ref_mut, seq = line.rstrip().split("\t")
            pepInfo[int(peptideID)] = [ref_mut, ENSP, mutation, int(varPos)]

def parseClassI(filename):
    sample = filename.split("/")[-2]    ## second to last is the sample name folder

    if sample not in resI:
        resI[sample] = {}

    INPREDICTIONS = False
    for line in open(filename, "r"):
        if INPREDICTIONS:
            if line.strip() == "":
                ## finished this chunk
                INPREDICTIONS = False
            elif not line.startswith("---"):
                try:
                    line = line.strip().rstrip().split()
                    pos = int(line[0])
                    hla = line[1]
                    peptide = line[2]
                    core = line[3]
                    Of = line[4]
                    Gp = line[5]
                    Gl = line[6]
                    Ip = line[7]
                    Il = line[8]
                    Icore = line[9]
                    identity = int(line[10])
                    score = line[11]
                    ic50 = float(line[12])
                    rank = line[13]
                except:
                    print("Error reading line")
                    print("File: {}".format(filename))
                    print("Line: {}".format(line))
                    sys.exit()


                ## Look up identity.
                ref_mut, ENSP, mut, varPos = pepInfo[identity]

                try:
                    protPos = int(mut[1:-1])
                except ValueError:
                    ## Likely selenocysteine (Sec) as the amino acid.
                    mutList = list(mut)
                    if is_number(mutList[1]):
                        ## reference is single amino acid
                        if "".join(mutList[len(mutList)-3:]) == "Sec":
                            mut = "".join(mutList[:len(mutList)-3]) + "U"
                        else:
                            sys.exit("Unknown amino acid: {}".format(pepInfo[identity]))
                    else:
                        ## reference is 3 letter amino acid name
                        if "".join(mutList[:3]) == "Sec":
                            mut = "U" + "".join(mutList[3:])
                        else:
                            sys.exit("Unknown amino acid: {}".format(pepInfo[identity]))


                ## make spot for data
                if ENSP not in resI[sample]:
                    resI[sample][ENSP] = {}
                if mut not in resI[sample][ENSP]:
                    resI[sample][ENSP][mut] = {}
                if hla not in resI[sample][ENSP][mut]:
                    resI[sample][ENSP][mut][hla] = {}
                #if ref_mut not in resI[sample][ENSP][mut][hla]:
                #    resI[sample][ENSP][mut][hla][ref_mut] = []

                ## have mutation protein position, variant position in submitted peptide, and short peptide position. Need to calculate if mutant position is in each short peptide.

                pepVarPos = varPos - pos + 1
                pepLen = len(peptide)

                if pepVarPos > 0 and pepVarPos <= pepLen:
                    ## variant is within peptide.

                    ## check that amino acid is correct (control for Selenocysteine U/Sec converting to X)
                    if (ref_mut == "ref" and peptide[pepVarPos-1] != mut[0] and mut[0] != "U") or (ref_mut == "mut" and peptide[pepVarPos-1] != mut[-1] and mut[-1] != "U"):
                        print("Mutation does not match peptide!")
                        print("NetMHC output: {}".format(line))
                        print("PepInfo: {}".format(pepInfo[identity]))
                        sys.exit()

                    if pos not in resI[sample][ENSP][mut][hla]:
                        resI[sample][ENSP][mut][hla][pos] = {}
                    if pepLen not in resI[sample][ENSP][mut][hla][pos]:
                        resI[sample][ENSP][mut][hla][pos][pepLen] = {}
                    resI[sample][ENSP][mut][hla][pos][pepLen][ref_mut] = [peptide, pepVarPos, ic50]


        elif line.strip().startswith("Pos"):
            INPREDICTIONS = True

def writeOutput(out):
    for samp in resI:
        for ensp in resI[samp]:
            for mut in resI[samp][ensp]:
                for hla in resI[samp][ensp][mut]:
                    for pos in resI[samp][ensp][mut][hla]:
                        for pepLen in resI[samp][ensp][mut][hla][pos]:
                            mutPepSeq = resI[samp][ensp][mut][hla][pos][pepLen]["mut"][0]
                            pepVarPos = resI[samp][ensp][mut][hla][pos][pepLen]["mut"][1]
                            mutIC50 = resI[samp][ensp][mut][hla][pos][pepLen]["mut"][2]
                            refIC50 = resI[samp][ensp][mut][hla][pos][pepLen]["ref"][2]

                            out.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(samp, ensp, mut, hla, mutPepSeq, pepVarPos, mutIC50, refIC50))

if __name__ == "__main__":

    ## Deal with command line arguments
    parser = argparse.ArgumentParser(description = "Parse NetMHCpan Output")
    ## add_argument("name", "(names)", metavar="exampleOfValue - best for optional", type=int, nargs="+", choices=[allowed,values], dest="nameOfVariableInArgsToSaveAs")
    parser.add_argument("--peptideInfoFiles", help = "File(s) containing details for submitted peptides", type = str, nargs = "+")
    parser.add_argument("--result_dirs", help = "Path to results", type = str, nargs = "+")
    parser.add_argument("outputFile", help = "File to write output to", type = str)
    parser.add_argument("-d", "--debug", action = "store_true", dest = "DEBUG", help = "Flag for setting debug/test state.")
    parser.add_argument("-v", "--verbose", action = "store_true", dest = "VERB", help = "Flag for setting verbose output.")
    args = parser.parse_args()

    ## Set Global Vars
    DEBUG = args.DEBUG
    VERB = args.VERB

    out = open(args.outputFile, "w")
    out.write("barcode\tensp\tmutation\thla\tmut_peptide\tpeptide_variant_position\tmut_ic50\tref_ic50\n")

    for peptideInfoFile, result_dir in zip(args.peptideInfoFiles, args.result_dirs):
        ## reset these for each variable
        resI = {}
        pepInfo = {}

        ## Parse peptide info file
        if VERB: print("Parsing peptide info file...")
        parsePeptideInfo(peptideInfoFile)

        ## Parse NetMHCpan output and write
        #parseClassI()
        out = open(args.outputFile, "a")


        if VERB: print("Parsing NetMHCpan results...")
        for root, dirs, filenames in os.walk(result_dir):
            for f in filenames:
                if f.endswith(".pMHC"):
                    if VERB: print("I: Parsing {}".format(os.path.join(root,f)))
                    parseClassI(os.path.join(root,f))

                    ## Write results
                    if VERB: print("Writing output...")
                    writeOutput(out)

                    ## clear resholder
                    resI = {}

        out.close()
    
    print("done.")
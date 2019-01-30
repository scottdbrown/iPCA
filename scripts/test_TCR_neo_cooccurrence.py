'''
Test TCR neoantigen co-occurrence. 
Look for neoantigens and TCR characteristics that appear together

Date: Marhc 20, 2017
@author: sbrown
'''

## Import Libraries
import sys
import argparse
import os

DEBUG = False
VERB = False

EXPRESSION_THRESH = 0.7

if __name__ == "__main__":

    ## Deal with command line arguments
    parser = argparse.ArgumentParser(description = "Test TCR neoantigen co-occurrence ")
    ## add_argument("name", "(names)", metavar="exampleOfValue - best for optional", type=int, nargs="+", choices=[allowed,values], dest="nameOfVariableInArgsToSaveAs")
    parser.add_argument("neoantigen_file", help = "File with neoantigens", type = str)
    parser.add_argument("tcr_file", help = "File with TCR info", type = str)
    parser.add_argument("tcr_mapping_file", help = "File with TCR sample mapping", type = str)
    parser.add_argument("output_file", help = "File to write output", type = str)
    parser.add_argument("-d", "--debug", action = "store_true", dest = "DEBUG", help = "Flag for setting debug/test state.")
    parser.add_argument("-v", "--verbose", action = "store_true", dest = "VERB", help = "Flag for setting verbose output.")
    args = parser.parse_args()

    ## Set Global Vars
    DEBUG = args.DEBUG
    VERB = args.VERB

    #print(args)
    ## arguments accessible as args.cmdline_arg or args.cmdflg or args.destName
    ## can test using parser.parse_args("-cmdflg value other_value".split())

    tcr_samps = set()

    print("Reading in TCR sample mapping...")
    cghub2bar = {}
    HEADER = True
    for line in open(args.tcr_mapping_file, "r"):
        if HEADER:
            HEADER = False
        else:
            line = line.split("\t")
            cghub2bar[line[5]] = line[4]
            tcr_samps.add(line[4])


    pep = {}
    pmhc = {}

    pmhcsamp = {}

    cdr3 = {}
    cdr3a = {}
    cdr3b = {}
    vgene = {}
    jgene = {}

    tcrsamp = {}
    cdr3asamp = {}
    cdr3bsamp = {}

    pep_cdr3 = {}
    pep_vgene = {}
    pep_jgene = {}

    pmhc_cdr3 = {}
    pmhc_vgene = {}
    pmhc_jgene = {}

    cdr3a_cdr3b = {}

    neo_samps = set()
    sampCount = set()
    print("Reading in neoantigens...")
    HEADER = True
    for line in open(args.neoantigen_file, "r"):
        if HEADER:
            HEADER = False
        else:
            line = line.rstrip().split("\t")
            if line[0][:16] in tcr_samps:
                sampCount.add(line[0][:16])
                if line[8] != "None":
                    exp = float(line[8])
                else:
                    exp = None

                if exp is None or exp > EXPRESSION_THRESH:
                    samp = line[0][:16]
                    neo_samps.add(samp)
                    hla = line[3]
                    peptide = line[4]

                    ## track peptide
                    if peptide not in pep:
                        pep[peptide] = set()
                    pep[peptide].add(samp)

                    ## track pmhc
                    pepmhc = "{}_{}".format(peptide, hla)
                    if pepmhc not in pmhc:
                        pmhc[pepmhc] = set()
                    pmhc[pepmhc].add(samp)

                    ## track this sample
                    if samp not in pmhcsamp:
                        pmhcsamp[samp] = []
                    pmhcsamp[samp].append(pepmhc)

    print("Reading in TCRs and measuring co-occurrence...")
    HEADER = True
    for line in open(args.tcr_file, "r"):
        if HEADER:
            HEADER = False
        else:
            ## get info from TCR file.
            line = line.split("\t")
            if line[0] in cghub2bar:
                samp = cghub2bar[line[0]]
                chain = line[2]
                vseg = line[4]
                aaSeq = line[5]
                jseg = line[10]

                if samp in neo_samps:

                    if samp not in cdr3asamp:
                        cdr3asamp[samp] = set()
                    if chain=="alpha": cdr3asamp[samp].add(aaSeq)

                    if samp not in cdr3bsamp:
                        cdr3bsamp[samp] = set()
                    if chain=="beta": cdr3bsamp[samp].add(aaSeq)

                    ## save information on recurrent cdr3/vseg/jseg
                    cdr3_chain = "{}_{}".format(aaSeq, chain)
                    if cdr3_chain not in cdr3:
                        cdr3[cdr3_chain] = set()
                    cdr3[cdr3_chain].add(samp)

                    if chain == "alpha":
                        if aaSeq not in cdr3a:
                            cdr3a[aaSeq] = set()
                        cdr3a[aaSeq].add(samp)
                    elif chain == "beta":
                        if aaSeq not in cdr3b:
                            cdr3b[aaSeq] = set()
                        cdr3b[aaSeq].add(samp)


                    if len(vseg.split(",")) == 1:
                        ## only one v gene reported
                        if vseg not in vgene:
                            vgene[vseg] = set()
                        vgene[vseg].add(samp)
                    else:
                        vseg = None

                    if len(jseg.split(",")) == 1:
                        ## only one j gene reported
                        if jseg not in jgene:
                            jgene[jseg] = set()
                        jgene[jseg].add(samp)
                    else:
                        jseg = None



                    ## look at co-occurrence. ##
                    for pepmhc in pmhcsamp[samp]:
                        peptide = pepmhc.split("_")[0]

                        # pep-cdr3
                        if (peptide, cdr3_chain) not in pep_cdr3:
                            pep_cdr3[(peptide, cdr3_chain)] = set()
                        pep_cdr3[(peptide, cdr3_chain)].add(samp)
                        # pep-v
                        if vseg != None:
                            if (peptide, vseg) not in pep_vgene:
                                pep_vgene[(peptide, vseg)] = set()
                            pep_vgene[(peptide, vseg)].add(samp)
                        # pep-j
                        if jseg != None:
                            if (peptide, jseg) not in pep_jgene:
                                pep_jgene[(peptide, jseg)] = set()
                            pep_jgene[(peptide, jseg)].add(samp)


                        # pmhc-cdr3
                        if (pepmhc, cdr3_chain) not in pmhc_cdr3:
                            pmhc_cdr3[(pepmhc, cdr3_chain)] = set()
                        pmhc_cdr3[(pepmhc, cdr3_chain)].add(samp)
                        # pmhc-v
                        if vseg != None:
                            if (pepmhc, vseg) not in pmhc_vgene:
                                pmhc_vgene[(pepmhc, vseg)] = set()
                            pmhc_vgene[(pepmhc, vseg)].add(samp)
                        # pmhc-j
                        if jseg != None:
                            if (pepmhc, jseg) not in pmhc_jgene:
                                pmhc_jgene[(pepmhc, jseg)] = set()
                            pmhc_jgene[(pepmhc, jseg)].add(samp)
                
    ## alpha-beta co-occurence:
    for samp in cdr3asamp:
        for a in cdr3asamp[samp]:
            for b in cdr3bsamp[samp]:
                if (a,b) not in cdr3a_cdr3b:
                    cdr3a_cdr3b[(a,b)] = set()
                cdr3a_cdr3b[(a,b)].add(samp)


    ## now print if there is overlap, and how significant it is.

    print("{} total samples.".format(len(sampCount)))

    ## each time more than one subject shares two things, print it, and also print how many subjects those things are found in individually.
    out = open(args.output_file, "w")
    out.write("type\tmutation\tTCR\tfound_together\tmutation_count\tTCR_count\n")
    if VERB: print("Co-occurrences of peptide - CDR3:")
    for (p,c) in pep_cdr3:
        if len(pep_cdr3[(p,c)]) > 1:
            if VERB: print("Peptide {} and CDR3 {} is shared by {} samples.".format(p, c, len(pep_cdr3[(p,c)])))
            if VERB: print("Peptide {} is found in {} samples.".format(p, len(pep[p])))
            if VERB: print("CDR3 {} is found in {} samples.".format(c, len(cdr3[c])))
            out.write("peptide-cdr3\t{}\t{}\t{}\t{}\t{}\n".format(p, c, len(pep_cdr3[(p,c)]), len(pep[p]), len(cdr3[c])))

    if VERB: print("Co-occurrences of peptide - Vgene:")
    for (p,v) in pep_vgene:
        if len(pep_vgene[(p,v)]) > 1:
            if VERB: print("Peptide {} and Vgene {} is shared by {} samples.".format(p, v, len(pep_vgene[(p,v)])))
            if VERB: print("Peptide {} is found in {} samples.".format(p, len(pep[p])))
            if VERB: print("Vgene {} is found in {} samples.".format(v, len(vgene[v])))
            out.write("peptide-vgene\t{}\t{}\t{}\t{}\t{}\n".format(p, v, len(pep_vgene[(p,v)]), len(pep[p]), len(vgene[v])))

    if VERB: print("Co-occurrences of peptide - Jgene:")
    for (p,j) in pep_jgene:
        if len(pep_jgene[(p,j)]) > 1:
            if VERB: print("Peptide {} and Jgene {} is shared by {} samples.".format(p, j, len(pep_jgene[(p,j)])))
            if VERB: print("Peptide {} is found in {} samples.".format(p, len(pep[p])))
            if VERB: print("Jgene {} is found in {} samples.".format(j, len(jgene[j])))
            out.write("peptide-jgene\t{}\t{}\t{}\t{}\t{}\n".format(p, j, len(pep_jgene[(p,j)]), len(pep[p]), len(jgene[j])))


    if VERB: print("Co-occurrences of pMHC - CDR3:")
    for (p,c) in pmhc_cdr3:
        if len(pmhc_cdr3[(p,c)]) > 1:
            if VERB: print("pMHC {} and CDR3 {} is shared by {} samples.".format(p, c, len(pmhc_cdr3[(p,c)])))
            if VERB: print("pMHC {} is found in {} samples.".format(p, len(pmhc[p])))
            if VERB: print("CDR3 {} is found in {} samples.".format(c, len(cdr3[c])))
            out.write("pmhc-cdr3\t{}\t{}\t{}\t{}\t{}\n".format(p, c, len(pmhc_cdr3[(p,c)]), len(pmhc[p]), len(cdr3[c])))

    if VERB: print("Co-occurrences of pMHC - Vgene:")
    for (p,v) in pmhc_vgene:
        if len(pmhc_vgene[(p,v)]) > 1:
            if VERB: print("pMHC {} and Vgene {} is shared by {} samples.".format(p, v, len(pmhc_vgene[(p,v)])))
            if VERB: print("pMHC {} is found in {} samples.".format(p, len(pmhc[p])))
            if VERB: print("Vgene {} is found in {} samples.".format(v, len(vgene[v])))
            out.write("pmhc-vgene\t{}\t{}\t{}\t{}\t{}\n".format(p, v, len(pmhc_vgene[(p,v)]), len(pmhc[p]), len(vgene[v])))

    if VERB: print("Co-occurrences of pMHC - Jgene:")
    for (p,j) in pmhc_jgene:
        if len(pmhc_jgene[(p,j)]) > 1:
            if VERB: print("pMHC {} and Jgene {} is shared by {} samples.".format(p, j, len(pmhc_jgene[(p,j)])))
            if VERB: print("pMHC {} is found in {} samples.".format(p, len(pmhc[p])))
            if VERB: print("Jgene {} is found in {} samples.".format(j, len(jgene[j])))
            out.write("pmhc-jgene\t{}\t{}\t{}\t{}\t{}\n".format(p, j, len(pmhc_jgene[(p,j)]), len(pmhc[p]), len(jgene[j])))

    if VERB: print("Co-occurrences of CDR3a - CDR3b:")
    for (a,b) in cdr3a_cdr3b:
        if len(cdr3a_cdr3b[(a,b)]) > 1:
            if VERB: print("CDR3a {} and CDR3b {} is shared by {} samples.".format(a, b, len(cdr3a_cdr3b[(a,b)])))
            if VERB: print("CDR3a {} is found in {} samples.".format(a, len(cdr3a[a])))
            if VERB: print("CDR3b {} is found in {} samples.".format(b, len(cdr3b[b])))
            out.write("cdr3a-cdr3b\t{}\t{}\t{}\t{}\t{}\n".format(a, b, len(cdr3a_cdr3b[(a,b)]), len(cdr3a[a]), len(cdr3b[b])))

    out.close()
    print("done.")
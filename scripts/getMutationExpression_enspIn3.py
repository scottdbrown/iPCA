'''
Get Mutation Expression
Get expression of gene-containing mutations.

Date:
@author: sbrown
'''

## Import Libraries
import sys
import argparse
import os
from biomart import BiomartServer
import numpy
import time

biomartURL = "www.ensembl.org/biomart"
biomartDataset = "hsapiens_gene_ensembl"

DEBUG = False
VERB = False


if __name__ == "__main__":

    ## Deal with command line arguments
    parser = argparse.ArgumentParser(description = "Get Mutation Expression")
    ## add_argument("name", "(names)", metavar="exampleOfValue - best for optional", type=int, nargs="+", choices=[allowed,values], dest="nameOfVariableInArgsToSaveAs")
    parser.add_argument("pMHC_file", help = "File with predicted pMHC", type = str)
    parser.add_argument("expression_file", help = "File with expression (EB++)", type = str)
    parser.add_argument("out_file", help = "File to write updated pMHC to", type = str)
    parser.add_argument("id_conv_file_new", help = "File to write ID conversions in the case of script failure", type = str)
    parser.add_argument("--id_conv_file_existing", dest = "ids", help = "File containing ensp-entrez conversions")
    parser.add_argument("-d", "--debug", action = "store_true", dest = "DEBUG", help = "Flag for setting debug/test state.")
    parser.add_argument("-v", "--verbose", action = "store_true", dest = "VERB", help = "Flag for setting verbose output.")
    args = parser.parse_args()

    ## Set Global Vars
    DEBUG = args.DEBUG
    VERB = args.VERB

    ensp2entrez = {}
    ## if have partial conversions already, read in.
    if args.ids:
        for line in open(args.ids, "r"):
            line = line.rstrip().split("\t")
            if line[1] != "None":
                ensp2entrez[line[0]] = int(line[1])
            else:
                ensp2entrez[line[0]] = None

        ## now write to new file (might be the same file)
        e2e = open(args.id_conv_file_new, "w")
        for ensp in ensp2entrez:
            ## write existing mapping to new file.
            e2e.write("{}\t{}\n".format(ensp, ensp2entrez[ensp]))
        e2e.close()

    ## connect to biomart
    server = BiomartServer(biomartURL)
    hge = server.datasets[biomartDataset]


    ## read in expression data
    print("Reading in expression data...")
    exp = {}
    subColMap = []
    genes = set()
    HEADER = True
    for line in open(args.expression_file, "r"):
        line = line.strip().split("\t")
        if HEADER:
            for item in line:
                item = item.replace('"','')
                exp[item[:16]] = {}
                exp[item[:16]][None] = [None]    ## set value if no Entrez can be found.
                subColMap.append(item[:16])
            HEADER = False
        else:
            gene = int(line[0].replace('"','').split("|")[1])
            genes.add(gene)
            for i in range(1,len(line)):
                if gene not in exp[subColMap[i]]:
                    exp[subColMap[i]][gene] = []
                if line[i] == "NA":
                    val = None
                else:
                    val = float(line[i])
                exp[subColMap[i]][gene].append(val)  ## handles if multiple values, can average later.

    print("Validate: Subject TCGA-OR-A5J1-01A Gene 100133144 Value: {} = 3.2661".format(exp["TCGA-OR-A5J1-01A"][100133144][0]))

    ## read in pMHC data
    out = open(args.out_file, "w")
    print("Reading in pMHC data...")
    HEADER = True
    for line in open(args.pMHC_file, "r"):
        if HEADER:
            out.write(line.rstrip() + "\texpressionEB++\n")
            HEADER = False
        else:
            linesplit = line.rstrip().split("\t")
            samp = linesplit[0][:16]
            ensp = linesplit[3]

            if samp in exp:

                if ensp not in ensp2entrez:
                    ## use biomart to get entrez
                    if VERB: print("Biomarting {}...".format(ensp))
                    GOT_ENTREZ = False
                    while not GOT_ENTREZ:
                        try:
                            entrez = int(hge.search({"filters":{"ensembl_peptide_id":ensp}, "attributes":["entrezgene"]}).text)
                            if entrez not in genes: ## not in EB++ file.
                                entrez = None
                            GOT_ENTREZ = True
                        except ValueError:
                            entrez = None
                            GOT_ENTREZ = True
                        except:
                            print("Error connecting to Biomart, trying again...")
                            time.sleep(10)

                    ensp2entrez[ensp] = entrez
                    e2e = open(args.id_conv_file_new, "a")
                    e2e.write("{}\t{}\n".format(ensp, entrez))
                    e2e.close()

                ## check if non-None values present.
                expSet = set(exp[samp][ensp2entrez[ensp]])
                if None in expSet:
                    expSet.remove(None)
                if len(expSet) > 0:
                    expList = [n for n in exp[samp][ensp2entrez[ensp]] if n != None]
                    expression = numpy.mean(expList)
                else:
                    expression = None

            else:
                expression = None

            out.write("{}\t{}\n".format(line.rstrip(), expression))

    out.close()
    print("done.")
'''
Prepare peptide binding predictions
Using provided Polysolver HLA calls
and MC3 validated mutation calls.

Date: August 12, 2016
@author: sbrown

Edited February 24, 2016:
    - Use OptiType calls
    - Use newest MC3 file with appropriate filters
Edited May 22, 2018:
    - Output genomic coords instead of ENSP.
'''

## Import Libraries
import sys
import argparse
import os

NETMHCPAN_BIN = "/home/sbrown/bin/netMHCpan-3.0/netMHCpan"

DEBUG = False
VERB = False

HLA_DICT = {}
PROT_DICT = {}
whitelist = set()
participantsMeetingFilter = set()
samplesMeetingFilter = set()


if __name__ == "__main__":

    ## Deal with command line arguments
    parser = argparse.ArgumentParser(description = "Prepare peptide binding predictions")
    ## add_argument("name", "(names)", metavar="exampleOfValue - best for optional", type=int, nargs="+", choices=[allowed,values], dest="nameOfVariableInArgsToSaveAs")
    parser.add_argument("hla_file", help = "File containing HLA calls", type = str)
    parser.add_argument("whitelist_file", help = "File containing whitelist sample IDs", type = str)
    parser.add_argument("mc3_file", help = "File containing mutation calls", type = str)
    parser.add_argument("ref_prot", help = "Reference proteome (matches genome version for mutation calls)", type = str)
    parser.add_argument("workingDir", help = "Path to directory to write files.", type = str)
    parser.add_argument("-d", "--debug", action = "store_true", dest = "DEBUG", help = "Flag for setting debug/test state.")
    parser.add_argument("-v", "--verbose", action = "store_true", dest = "VERB", help = "Flag for setting verbose output.")
    args = parser.parse_args()

    ## Set Global Vars
    DEBUG = args.DEBUG
    VERB = args.VERB


    ## Read in HLA calls
    ## {TCGA-BAR-CODE: set(HLA alleles)}

    ## format (no header): CANCER-TCGA-AA-AAAA \t A1 A2 B1 B2 C1 C2
    print("Reading in HLA information...")
    for line in open(args.hla_file, "r"):
        line = line.rstrip().split("\t")
        aliquot = line[0]
        sub = aliquot[:12]
        hlas_temp = line[1:7]
        hlas = set()
        ## fix hla nomenclature
        for h in hlas_temp:
            h = h.replace("*","")
            h = "HLA-" + h
            hlas.add(h)


        if sub not in HLA_DICT:
            HLA_DICT[sub] = hlas
        else:
            print("Subject {} already has HLAs entered.\nExisting: {}\nNew:{}".format(sub, HLA_DICT[sub], hlas))


    print("Reading in whitelist...")
    for line in open(args.whitelist_file, "r"):
        whitelist.add(line.rstrip().replace('"',''))


    ## Read in reference proteome
    ## {ENSP: SEQUENCE}
    print("Reading in reference proteome...")
    for line in open(args.ref_prot, "r"):
        if line.startswith(">"):
            header = line.rstrip()
            prot = header.split(" ")[0]
            ## cleave off ">"
            prot = prot[1:]
            ## remove ".#" at end
            prot = prot.split(".")[0]
            if prot not in PROT_DICT:
                PROT_DICT[prot] = ""
            else:
                sys.exit("Protein {} already has sequence:\n\n{}\n\nExiting.".format(prot, PROT_DICT[prot]))
        else:
            PROT_DICT[prot] += line.rstrip()

    #if DEBUG: print(PROT_DICT)

    ## Read through MC3 file, processing if subject has HLA calls
    ## header (of previous file. may not be exact):
    ## [0]Hugo_Symbol [1]Entrez_Gene_Id [2]Center [3]NCBI_Build [4]Chromosome [5]Start_Position [6]End_Position  [7]Strand [8]Variant_Classification [9]Variant_Type [10]Reference_Allele [11]Tumor_Seq_Allele1 [12]Tumor_Seq_Allele2 [13]dbSNP_RS [14]dbSNP_Val_Status [15]Tumor_Sample_Barcode [16]Matched_Norm_Sample_Barcode [17]Match_Norm_Seq_Allele1 [18]Match_Norm_Seq_Allele2 [19]Tumor_Validation_Allele1 [20]Tumor_Validation_Allele2 [21]Match_Norm_Validation_Allele1 [22]Match_Norm_Validation_Allele2 [23]Verification_Status [24]Validation_Status [25]Mutation_Status [26]Sequencing_Phase [27]Sequence_Source [28]Validation_Method [29]Score [30]BAM_File [31]Sequencer [32]Tumor_Sample_UUID [33]Matched_Norm_Sample_UUID [34]HGVSc [35]HGVSp [36]HGVSp_Short [37]Transcript_ID [38]Exon_Number [39]t_depth [40]t_ref_count [41]t_alt_count [42]n_depth [43]n_ref_count [44]n_alt_count [45]all_effects [46]Allele [47]Gene [48]Feature [49]Feature_type [50]Consequence [51]cDNA_position [52]CDS_position [53]Protein_position [54]Amino_acids [55]Codons [56]Existing_variation [57]ALLELE_NUM [58]DISTANCE [59]STRAND [60]SYMBOL [61]SYMBOL_SOURCE [62]HGNC_ID [63]BIOTYPE [64]CANONICAL [65]CCDS [66]ENSP [67]SWISSPROT [68]TREMBL [69]UNIPARC [70]RefSeq [71]SIFT [72]PolyPhen [73]EXON [74]INTRON [75]DOMAINS [76]GMAF [77]AFR_MAF [78]AMR_MAF [79]ASN_MAF [80]EAS_MAF [81]EUR_MAF [82]SAS_MAF [83]AA_MAF [84]EA_MAF [85]CLIN_SIG [86]SOMATIC [87]PUBMED [88]MOTIF_NAME [89]MOTIF_POS [90]HIGH_INF_POS [91]MOTIF_SCORE_CHANGE [92]IMPACT [93]PICK [94]VARIANT_CLASS [95]TSL [96]HGVS_OFFSET [97]PHENO [98]MINIMISED [99]ExAC_AF [100]ExAC_AF_AFR [101]ExAC_AF_AMR [102]ExAC_AF_EAS [103]ExAC_AF_FIN [104]ExAC_AF_NFE [105]ExAC_AF_OTH [106]ExAC_AF_SAS [107]GENE_PHENO [108]FILTER [109]COSMIC [110]CENTERS [111]CONTEXT [112]DBVS [113]NCALLERS

    print("Reading in mutation information...")

    peps = {}
    pepsInfo = {}
    PEP_NUM = 0


    ## Filters:
    ## Exclude LAML
    ## Retain FILTER %in% PASS, wga, native_wga_mix
    ## NCALLERS > 1
    ## Retain coding mutations, defined as  Variant_Classification IN ("Frame_Shift_Del","Frame_Shift_Ins","In_Frame_Del","In_Frame_Ins","Missense_Mutation","Nonsense_Mutation","Nonstop_Mutation","Splice_Site","Translation_Start_Site")

    HEADER = True
    for line in open(args.mc3_file, "r"):
        line = line.rstrip().split("\t")
        if HEADER:
            i = 0
            for key in line:
                if key == "Chromosome":
                    chromi = i
                elif key == "Start_Position":
                    posi = i
                elif key == "Variant_Classification":
                    varclassi = i
                elif key == "Variant_Type":
                    vartypei = i
                elif key == "Tumor_Sample_Barcode":
                    tumouri = i
                elif key == "Matched_Norm_Sample_Barcode":
                    normali = i
                elif key == "HGVSp_Short":
                    pMuti = i
                elif key == "Protein_position":
                    protPosi = i
                elif key == "Amino_acids":
                    refmuti = i
                elif key == "ENSP":
                    enspi = i
                elif key == "FILTER":
                    filteri = i
                elif key == "NCALLERS":
                    ncallersi = i

                i += 1
            HEADER = False

        else:
            ## If missense variant
            if line[varclassi] == "Missense_Mutation" and line[vartypei] == "SNP" and len(line[refmuti].split("/")) == 2:
                if DEBUG: print(line)
                tumour = line[tumouri]
                normal = line[normali]
                sub = tumour[:12]
                pMut = line[pMuti][2:] ## cleaves off leading "p."
                protPos = int(line[protPosi])
                ref,mut = line[refmuti].split("/")
                ensp = line[enspi]

                ## meets filtering criteria
                if line[filteri] in ["PASS","wga","native_wga_mix"] and int(line[ncallersi]) > 1 and sub in whitelist:

                    participantsMeetingFilter.add(sub)
                    samplesMeetingFilter.add(tumour)

                    ## if we have HLA information
                    if sub in HLA_DICT:

                        ## get refernce protein sequence
                        if ensp in PROT_DICT:
                            ensp_seq = PROT_DICT[ensp]

                            ## check that reference amino acid sequence matches
                            #if DEBUG: print("protPos = {}".format(protPos))
                            protRefAA = ensp_seq[protPos-1:protPos]
                            if ref != protRefAA:
                                print("Reference allele does not match reference sequence: {} - {}: {}/{}".format(tumour, ensp, ref, protRefAA))
                            else:
                                ## get flanking peptide sequence
                                lowerPepBound = protPos - 10 - 1
                                upperPepBound = protPos + 10

                                ## check edge cases
                                if lowerPepBound < 0:
                                    lowerPepBound = 0
                                if upperPepBound > len(ensp_seq):
                                    upperPepBound = len(ensp_seq)

                                pepPos = lowerPepBound + 1

                                refPeptide = ensp_seq[lowerPepBound:upperPepBound]
                                varPos = protPos - lowerPepBound

                                mutPeptide = refPeptide[0:varPos-1] + mut + refPeptide[varPos:]

                                ## add peptide to list.
                                if tumour not in peps:
                                    peps[tumour] = {}
                                    pepsInfo[tumour] = {}
                                ## Reference
                                PEP_NUM += 1
                                peps[tumour][PEP_NUM] = refPeptide
                                pepsInfo[tumour][PEP_NUM] = "{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(line[chromi], line[posi], ensp, pMut, varPos, "ref", refPeptide)
                                ## Mutant
                                PEP_NUM += 1
                                peps[tumour][PEP_NUM] = mutPeptide
                                pepsInfo[tumour][PEP_NUM] = "{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(line[chromi], line[posi], ensp, pMut, varPos, "mut", mutPeptide)

                        else:
                            print("No protein sequence available for {}:\n{}".format(ensp, line))
                    else:
                        if VERB: print("No HLA info available for {}.".format(sub))

                #else:
                    #print("Bypassing tumour {}.\nFILTER: {}\tNCALLERS: {}\tSUB: {}".format(tumour, line[filteri], line[ncallersi], sub))

    print("Number of participant IDs meeting maf filters: {}".format(len(participantsMeetingFilter)))
    print("Number of tumour IDs meeting maf filters: {}".format(len(samplesMeetingFilter)))

    ## make directories and write files
    print("Writing directories and files...")
    transferfof = open(os.path.join(args.workingDir, "filesToTranfer.fof"), "w")
    analysisfile = open(os.path.join(args.workingDir, "analysisToRun.sh"), "w")

    for tum in peps:
        os.makedirs(os.path.join(args.workingDir, tum))
        sub = tum[:12]

        ## get HLA
        hlas = ",".join(HLA_DICT[sub])

        ## write peptides
        out = open(os.path.join(args.workingDir, tum, "peptides.fa"), "w")
        toWrite = ""
        for k in peps[tum]:
            toWrite += ">{}\n{}\n".format(k, peps[tum][k])
        out.write(toWrite)
        out.close()

        transferfof.write("{}\n".format(os.path.join(args.workingDir, tum, "peptides.fa")))
        analysisfile.write("{}\tsource /home/sbrown/bin/pythonvenv/python3/bin/activate;{} -tdir tmpdirXXXXXX -a {} -f peptides.fa > bindingRes.pMHC;\n".format(tum, NETMHCPAN_BIN, hlas))

    transferfof.close()
    analysisfile.close()


    idsfile = open(os.path.join(args.workingDir, "peptideIDs.tsv"), "w")
    idsfile.write("tumour\tminPeptideID\tchromosome\tposition\tENSP\tmutation\tvarPos\tref_mut\tsequence\n")
    for tum in pepsInfo:
        for k in pepsInfo[tum]:
            idsfile.write("{}\t{}\t{}".format(tum,k,pepsInfo[tum][k]))
    idsfile.close()

    print("done.")



# TCGA immunePanCanAtlas

Scott Brown; sbrown@bcgsc.ca

### Contents:

[Neoantigens](#neoantigen-predictions-for-snvs)

[Figures](#figure-stuff)

## Neoantigen Predictions for SNVs

1. **Prepare peptide-MHC binding predictions.**

   **Script:** [prepPeptideBindingPredictions_final_genomicCoord.py](scripts/prepPeptideBindingPredictions_final_genomicCoord.py)

   **Input:** OptiType calls, list of whitelisted samples, maf file, proteome reference.

   **Output:** Directories, files, and shell scripts to perform NetMHCpan 3.0 binding predictions.

   **Brief Description:** For each SNV which meets filtering criteria, the relevant reference protein sequence is determined, and a 21mer peptide sequence with the mutant amino acid at the center is created. 

   **Environment:** Python 3.

   **Example Command:** `python scripts/prepPeptideBindingPredictions_final_genomicCoord.py ../../data/optitype_hla/OptiTypeCallsHLA021917.tsv ../../data/Whitelisted_doNotUseIsFalse_samples_syn4906913.tsv ../../data/mc3.v0.2.8.CONTROLLED.maf /projects/sbrown_prj/selfpeptidome/data/Homo_sapiens.GRCh37.pep.all.fa /current/directory/`

2. **Run peptide-MHC binding predictions.**

   **Script:** [clusterTAS](https://github.com/scottdbrown/bcgsc-scripts/blob/master/clusterTAS)

   **Input:** filesToTransfer.fof and analysisToRun.sh file produced during Step 1.

   **Output:** analysis_results directory with a directory for each sample and NetMHCpan 3.0 output for each sample.

   **Brief Description:** clusterTAS is a custom management script/tool to manage and run submissions of jobs on the Genome Science Centre Genesis cluster. The analysisToRun.sh file contains the specifics for the NetMHCpan analysis, which can be run on your system (perhaps with some modification).

3. **Parse the NetMHCpan output.**

   **Script:** [parseNetMHCpanOutput_withCoord.py](scripts/parseNetMHCpanOutput_withCoord.py)

   **Input:** peptideIDs.tsv from Step 1, analysis_results/ dir from Step 2.

   **Output:** List of peptide-MHC predictions (including non-binders) for all samples.

   **Brief Description:** This reads through all the results files, parses the NetMHCpan output and extract relevant information, and connects this to mutation and peptide information. As 21mers are used as input for NetMHCpan binding predictions, not all resulting peptide-MHC combinations will contain the mutant amino acid (for example, an 8mer starting at the first position of the 21mer will not reach the variant position). This script will check the position of the variant amino acid, and only output peptide-MHC combinations which contain the variant.

   **Environment:** Python 3

   **Example Command:** `python scripts/parseNetMHCpanOutput_withCoord.py TCGA_peptideMHC_binding_MC3_v0.2.8.CONTROLLED_genomicCoord_170403.tsv --peptideInfoFiles peptideBinding_feb2017_gc/peptideIDs.tsv peptideBinding_mar2017_gc/peptideIDs.tsv peptideBinding_mar2017_2_gc/peptideIDs.tsv --result_dirs peptideBinding_feb2017/results/analysis_results/ peptideBinding_mar2017/results/analysis_results/ peptideBinding_mar2017_2/results/analysis_results/`

4. **Connect the above data to expression data.**

   **Script:** [getMutationExpression_enspIn3.py](scripts/getMutationExpression_enspIn3.py)

   **Input:** List of all pMHCs produced in Step 3, EB++ PanCan RNASeqV2 Expression Matrix, ensp2entrez gene name conversion file

   **Output:** List of all pMHCs with expression value for source gene added as a column.

   **Brief Description:** This looks up the expression of the source gene for each pMHC. Biomart is used to convert gene names, but some cases occur where it fails, so a separate list of ensp2entrez identifiers can be provided.

   **Environment:** Python 3

   **Example Command:** `python scripts/getMutationExpression_enspIn3.py TCGA_peptideMHC_binding_MC3_v0.2.8.CONTROLLED_genomicCoord_170403.tsv ../data/expression/EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv TCGA_peptideMHC_binding_MC3_v0.2.8.CONTROLLED_genomicCoord_withExpression_170403.tsv ensp2entrez.tsv --id_conv_file_existing ensp2entrez.tsv`

5. **Parse out the binders.**

   **Example Command:** 

   `$ head -1 TCGA_peptideMHC_binding_MC3_v0.2.8.CONTROLLED_genomicCoord_withExpression_170403.tsv > TCGA_peptideMHC_binding_MC3_v0.2.8.CONTROLLED_genomicCoord_predictedBindersOnly_withExpression_170404.tsv`

   `$ awk '{if($8 < 500){print}}'
   TCGA_peptideMHC_binding_MC3_v0.2.8.CONTROLLED_genomicCoord_withExpression_170403.tsv >> TCGA_peptideMHC_binding_MC3_v0.2.8.CONTROLLED_genomicCoord_predictedBindersOnly_withExpression_170404.tsv`

6. **Create sample-level summary.**

   **Script:** [summarizeSampleLevelEpitopes.py](scripts/summarizeSampleLevelEpitopes.py)

   **Input:** List of pMHCs with expression data from Step 4.

   **Output:** Sample-level statistics and counts for pMHC data.

   **Brief Description:** Provides counts of pMHC per subject.

   **Environment:** Python 3

   **Example Command:** `python scripts/summarizeSampleLevelEpitopes.py TCGA_peptideMHC_binding_MC3_v0.2.8.CONTROLLED_withExpression_170403.tsv ../data/Whitelisted_doNotUseIsFalse_samples_syn4906913.tsv ../data/aliquot_blacklist.txt TCGA_pMHC_SNV_sampleSummary_MC3_v0.2.8.CONTROLLED_170404.tsv`



## Figure Stuff

1. **Calculate TCR/pMHC co-occurrence.**

   **Script:** [test_TCR_neo_cooccurrence.py](scripts/test_TCR_neo_cooccurrence.py)

   **Input:** Neoantigen file from above, TCR file ("CGhub ID", "TCR chain", "V gene", "CDR3 amino acid", "J gene" are the required fields, expected in columns 1, 3, 5, 6, 11 (1-based).), TCR sample file (used to map CGhub ID (column 6) to TCGA barcode (column 5)).

   **Output:** File with TCR-pMHC co-occurrence

   **Brief Description:** For each SNV which meets filtering criteria, the relevant reference protein sequence is determined, and a 21mer peptide sequence with the mutant amino acid at the center is created. 

   **Environment:** Python 3.

   **Example Command:** `python scripts/test_TCR_neo_cooccurrence.py TCGA_peptideMHC_binding_MC3_v0.2.8.CONTROLLED_genomicCoord_predictedBindersOnly_withExpression_170404.tsv ../data/mitcr_complete/TCGA_mitcr_cdr3_result_161008.tsv ../data/mitcr_complete/mitcr_sampleStatistics_20160714.tsv mut_tcr_cooccurrence_170322.tsv`

2. **Plot Figures 5D and 5E (TCR-pMHC co-occurrence).**

   **Script:** [mutation_tcr_cooccurrence.R](scripts/mutation_tcr_cooccurrence.R)

   **Input:** `mut_tcr_cooccurrence_170322.tsv` from Step 1.

3. **Calculate pMHC recurrence.**

   **Script:** [findRecurrent_pMHC.py](scripts/findRecurrent_pMHC.py)

   **Input:** pMHC list and sample summary list from Neoantigen Prediction section above.

   **Output:** Lists and matrices of recurrent binding peptides and pMHCs.

   **Brief Description:**  Determines the number of samples that each presented peptide or pMHC combination are found in.

   **Environment:** Python 3.

   **Example Command:** `python scripts/findRecurrent_pMHC.py ../TCGA_peptideMHC_binding_MC3_v0.2.8.CONTROLLED_genomicCoord_predictedBindersOnly_withExpression_170404.tsv ../TCGA_pMHC_SNV_sampleSummary_MC3_v0.2.8.CONTROLLED_170404.tsv recurrentBindingPeptides_170530.tsv recurrentBindingPeptidesMatrix_170530.tsv recurrentPMHC_170530.tsv recurrentPMHCmatrix_170530.tsv`

4. **Plot Figures 5A and 5B (recurrent pMHCs).**

   **Script:** [neoantigen_recurrence.R](scripts/neoantigen_recurrence.R)

   **Input:** `recurrentPMHC_170530.tsv` from Step 3.
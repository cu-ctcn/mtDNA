# mtDNA

This repository contains R and shell scripts that were used to generate the results, figures, and tables presented in the manuscript [“Characterization of mitochondrial DNA quantity and quality in the human aged and Alzheimer’s disease brain”](https://doi.org/10.1186/s13024-021-00495-8) published in Molecular Neurodegeneration. The data itself is not stored in this repository but can be accessed via the AD Knowledge portal. Please visit the [manuscript’s landing page](https://doi.org/10.7303/syn25618990) for more details about the datasets.

## Analysis, figures, and tables
The R scripts in the main folder were used to generate results, figures and tables shown in the manuscript.
* Figure 1: *The mtDNAcn is reduced in cortical brain regions in AD.* (fig1_mtDNAcnRaw.R)
* Figure 2: *Changes of the mtDNAcn are primarily associated with tau in the DLPFC and TDP-43 in the PCC.* (fig2_mtDNAcnAndPathologies.R)
* Figure 3: *The mtDNAcn is associated with cognitive functioning independent of brain pathologies.* (fig3_cognitionMultivariable.R)
* Figure 4: *Frequency of mtDNA heteroplasmic mutations in cortical regions increases with age.* (fig4_heteroplasmy.R, fig4B_circusPlot.R)
* Figure 5: *Correlates of mitochondrial health demonstrate complex relationship with AD-related traits in the DLPFC.* (fig5_mtContent.R)
* Table 1: *Characteristics of the ROSMAP cohort.* (table1_cohorts.R)
* Table 2: *Effect of ApoE e4 genotype on mtDNAcn.* (table2_apoe4.R)
* Supplementary Figure 1: *Log-transformed mtDNAcn values are approximately normally distributed.* (fig1_mtDNAcnRaw.R)
* Supplementary Figure 2: *Changes of the mtDNAcn are primarily associated with tau in the DLPFC.* (fig2_mtDNAcnAndPathologies.R)
* Supplementary Figure 3: *Significant fraction of the ApoE ε4 effect on mtDNAcn is mediated via AD pathologies.* (table2_apoe4.R)
* Supplementary Figure 4: *Number of mtDNA heteroplasmic mutations in cortical regions is associated with age adjusted for sex, pathologic diagnosis, and mtDNAcn.* (fig4_heteroplasmy.R)
* Supplementary Figure 5: *Mitochondrial content score and its relation to mtDNAcn, mtDNA heteroplasmy burden and AD-related phenotypes.* (fig5_mtContent.R)
* Supplementary Table 1: *Characteristics of the Mayo cohort.* (table1_cohorts.R)
* Supplementary Table 2: *Characteristics of the MSBB cohort.* (table1_cohorts.R)
* Supplementary Table 3: *Detailed results from the univariate regression analyses shown in Fig. 2A.* (fig2_mtDNAcnAndPathologies.R) 
* Supplementary Table 4: *Association between cell type proportions and mtDNAcn in the DLPFC.* (fig2_mtDNAcnAndPathologies.R)
* Supplementary Table 5: *Association between number of mtDNA heteroplasmies and age adjusted for pathologic diagnosis.* ( fig4_heteroplasmy.R)
* Supplementary Excel File 2: *Genetic associations with mtDNAcn* (table2_apoe4.R)

## mtDNAcn estimation
The mtDNAcn was calculated using the R script located in the mtDNAcn folder.

## mtDNA variant calling
mtDNA homoplasmic and heteroplasmic variants were detected using the scripts in the subfolder mtVariantsPipeline. The pipeline replicates the [GATK mitochondria pipeline](https://github.com/gatk-workflows/gatk4-mitochondria-pipeline). See the info file in the folder for more information how to invoke the scripts.

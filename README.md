Script for Long et al. (2025) Evolutionary genomics unravels the responses and adaptation to climate change in a key alpine forest tree species. Molecular Biology and Evolution, accepted.


## **1.Population structure** 
This section includes population structure analyses (e.g., Admixture and PCA), gene flow analyses (e.g., Identity-by-Descent and Treemix), and genetic divergence and diversity assessments.

- ***Recent Demographic history***

    1.Population_structure/Demographic_history_PSMC/smc_run_split.sh was used to simulate the population demographic history and the results were visualized by smc_split_plot.sh script.
- ***Geome-wide genetic parameters calculation***

    The scripts under 1.Population_structure/Genome_wide_dxy_fst_pi_ld/ directory were used for estimating the LD, FST, DXY and genetic diversity.
- ***Identity by descent***

    The shared segments that have been inherited from a common ancestor among individuals were measured by 1.Population_structure/Identity_by_descent/1_IBD.sh and 2_IBD.R was used to visualize the results.

- ***Population structure***

    Admixture analysis was performed by 1.Population_structure/Population_structure/Admixture.sh and PCA analyis was run by PCA.sh.

- ***Gene Flow***

    VCF file was transfered to input format for Treemix using 1_plink2treemix.py and Treemix was for estimating the geneflow among lineages using 1.Population_structure/Treemix/2_treemix.sh

## **2.Identifying selection regions**

Local PCA is used to identify genomic regions with abnormal relatedness patterns. This is complemented with analyses such as XP-CLR, BetaScan, dXY, FST, and measures of genetic diversity to provide evidence for both balancing and divergent selection.

## **3.GEA and haploblocks**

GEA (Genotype-Environment Association) analyses are conducted using LFMM and GEMMA. Local PCA is also employed to identify haploblocks.

- ***GEA Analysis***

    3.Haploblocks_GEA/GEA_LFMM_GEMMA/ENV_values_extract.R was used to extract the enviromental values for the populations

    LFMM.sh and GEMMA.sh were used to perform GEA analysis and the results were visualized by GEA_mahhaton_plot.R
- ***Gradient-Forest Analysis***

    3.Haploblocks_GEA/Gradient_forest/1.get_variaton_frq.py was used to calculate the allele frequency for each population and culmu_importance.R was used to estimate the response of genetic variation to the bioclimatic gradient.
- ***Haploblock_detection***
    
    3.Haploblocks_GEA/Haploblocks/Haploblocks.sh was used to detect Haploblocks and pca_mds.R was used to identify high-quality Haploblocks by using function in pca_adaptive.R and vcf2het.pl scripts


## **4.Differential Expression Gene Analyses**

    4.Plasticity_divergence/DESeq.R was used for the identification and analysis of differentially expressed genes (DEGs).

## **5.Local+Forward_Reverse Genetic offest**
- ***Mutation load***

    5.Gnomic_prediction/Mutation_load/SIFT_annotation/1_get_config.sh was used to create the config file and then 2_sift.mkdb.sh was for making the database. 3_annotation.sh was used tpo predict whether an amino acid substitution affects protein function based on sequence homology and the physical properties of amino acids.
- ***Genetic offset***

    5.Gnomic_prediction/Genetic_offset/Forwrad_genetic_offset.R was used for the forwad genetic offset by defining different geographic distances.

    Local_Reverse_genetic_offset.R was used to calculate local genetic offset, which assumes in situ tolerance and reverse genetic offset to assess the maladaptation of populations when simultaneously considering the contributions of in situ adaptation and migration to future shifting climates.

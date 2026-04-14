# Medeiros_Almeida_Biogeography_Disjunct_Polypodiaceae

GENERAL INFORMATION

Title of Dataset: Supporting files for the article 'Medeiros MB, Almeida TE. 2026. Evolutionary and Biogeographic History of Disjunct Species of Polypodiaceae Between Tropical America and the Afrotropics. Acta Botanica Brasilica 40: XXXXX. doi: xxxxx'.

Authors: Matheus Bento Medeiros, M.Sc & Thaís Elias Almeida, PhD

Contact Email: matheusbento@live.com and thais.elias@ufpe.br

Data Collection Period: 2020-2024

Institutions: Universidade Federal do Oeste do Pará (UFOPA) and Universidade Federal de Pernambuco (UFPE)

Keywords: Afrotropic, Biogeography, Neotropic, Phylogenetic Analysis, Polypodiaceae

Funding Source: CAPES

Overview of Data and Files:
This dataset accompanies the article 'Evolutionary and Biogeographic History of Disjunct Species of Polypodiaceae Between Tropical America and the Afrotropics'. The provided files contain essential information for the biogeographic and phylogenetic analysis of Polypodiaceae species in the studied bioregions.

Included Files:

Concatenated_alignment_rbcl+rps4+trnlf+atpb.nex
Description: Alignment of DNA sequences from five plastid regions: 1) the rbcL gene, 2) the rps4 gene and the rps4-trnS intergenic spacer, 3) the trnL intron with short portions of the trnL and trnF genes, 4) the trnL-trnF intergenic spacer, and 4) the atpB gene.
Format: Text file (.nex)
Function: Used as an input for Bayesian inference analysis conducted with BEAST software.
Creation Date: 2025-10-12

Distribution map of the study species_ Input data for the Infomap bioregion.tiff
Description: Image showing the global distribution of the target species.
Format: Image file (.tiff)
Function: Data used to delimit bioregions with the Infomap Bioregions web application.
Creation Date: 2023-06-18

Polypodiaceae_estimate_time_divergence.trees
Description: Raw output file from the Bayesian inference analysis using BEAST; all 10,000 trees before burnin.
Format: Phylogenetic tree file (.trees)
Function: Used to generate a Maximum Credibility Tree used an input for biogeographic analysis in BioGeoBEARS
Creation Date: 2026-03-19

Polypodiaceae_MaximumCredibilityTree.tre
Description
Format: Phylogenetic tree file (.tre)
Function: Used as an input for biogeographic analyses in BioGeoBEARS.
Creation Date: 2026-03-19

Polypodiaceae_Multispecies_ocurrence_Table.txt
Description: Presence-absence matrix of all target species occurrences in each coded bioregion, downloaded from GBIF. The DOI of each download can be found as supplementary material of the article.
Format: Text file (.txt)
Function: Used as an input for biogeographic analysis in BioGeoBEARS.
Creation Date: 2023-09-11

resBAYAREALIKE+J_50BSMs_v1.pdf
Description: The 50 resulting trees from the Biogeographical Stochastic Mapping (BSM) analysis.
Format: Portable Document Format (.pdf)
Function: Output of the Biogeographical Stochastic Mapping (BSM) analysis, for visualization.
Creation Date: 2026-03-17

resBAYAREALIKE+J_50BSMs_v1.gif
Description: It is the animated visualization of the 50 trees from the file resBAYAREALIKE+J_50BSMs_v1.pdf.
Format: Graphic Interchange Format (.gif)
Function: Illustrative, for visualization.
Creation Date: 2026-03-17

Script_BioGeoBEARs_BSM
Description: R code.
Format: R Script (.R)
Function: R code used for run BioGeoBEARS (BioGeography with Bayesian (and likelihood) Evolutionary Analysis in R Scripts) and Biogeographical Stochastic Mapping (BSM) analyses.
Creation Date: 2026-03-17

Citation:
'edeiros, M.B.& Almeida T.E.(2025). Evolutionary and Biogeographic History of Disjunct Species of Polypodiaceae between the Neotropics and the Afrotropics. SciELO Data. DOI: https://doi.org/10.48331/scielodata.YHVLLZ'

Methodological Information:

Data Collection:
The dataset was delineated based on the floristic hypotheses of disjunct species proposed by Moran and Smith (2001), which describe patterns of morphological similarity among species of the Polypodiaceae family in neotropical and afrotropical regions. We obtained genetic sequences from five plastid regions (rbcL, rps4, rps4-trnS, trnL-trnF, and atpB) from the GenBank database to test these hypotheses. Geographic occurrence data were extracted from the Global Biodiversity Information Facility (GBIF) platform and processed to remove duplicates and records with incomplete or incorrect information.

Data Processing:
DNA sequences were individually aligned using MEGA XI software with the MUSCLE algorithm. Manual adjustments were made to the alignments, and the best substitution models were selected using jModelTest v2.1.10. Subsequently, the sequences were concatenated using SequenceMatrix v1.9. For phylogenetic inferences, we utilized BEAUti to configure the alignments, substitution models, and time calibrations based on two fossils. Inferences were conducted with BEAST, employing the relaxed clock model and Birth-Death tree. After running 10 million generations, the generated trees were evaluated using Tracer v1.7.1, discarding the first 10% as burn-in. The maximum credibility topology was then obtained with TreeAnnotator, and this topology was used in subsequent biogeographic inferences.

For the biogeographic analysis, species occurrence data were refined using the Wallace v2.0.5 package in R and QGIS. After processing, a presence-absence matrix of species by bioregion was generated and used as input for the biogeographic analysis in BioGeoBEARS. We employed the DEC, DIVALIKE, and BAYAREALIKE models, with comparisons based on the corrected Akaike Information Criterion (AICc), and biogeographic events were estimated through Stochastic Biogeographic Mapping (BSM).

Required Software:

MEGA XI (Koichiro et al. 2021) – Used for aligning genetic sequences with the MUSCLE algorithm.
jModelTest v2.1.10 (Darriba et al. 2012) – Used to select the best nucleotide substitution models.
SequenceMatrix v1.9 (Vaidya et al. 2011) – Used to concatenate aligned sequences.
BEAST v2.7.3 (Bouckaert et al. 2019) – Used for phylogenetic inference and divergence time estimation.
Tracer v1.7.1 (Rambaut et al. 2018) – Used to assess MCMC run convergence.
TreeAnnotator v2.7.3 (Bouckaert et al. 2019) – Used to generate the final phylogenetic tree topology.
R v4.1.0 (R Core Team 2020) and RStudio v1.4.1106 (RStudio Team 2022) – Used to execute the biogeographic analysis with the BioGeoBEARS package.
BioGeoBEARS (Matzke 2013) – R package used to conduct biogeographic analyses (DEC, DIVALIKE, BAYAREALIKE).
Wallace v2.0.5 (Kass et al. 2022) – Used to filter and refine species occurrence data.
QGIS v3.28 (QGIS Development Team 2023) – Used to visualize and adjust geographic distribution data.
CIPRES Science Gateway (Miller et al. 2010) – Platform used to run BEAST analyses remotely.
Specific Data Information

Additional supplementary files are available with the published article.


References:

Bouckaert R, Vaughan TG, Barido-Sottani J, et al. 2019. BEAST 2.5: An advanced software platform for Bayesian evolutionary analysis. PLoS Computational Biology 15(4): e1006650. DOI: 10.1371/journal.pcbi.1006650

Benton MJ, Wilf P, Sauquet H. 2022. The Angiosperm Terrestrial Revolution and the origins of modern biodiversity. New Phytologist 233: 2017–2035. DOI: 10.1111/nph.17822

Darriba D, Taboada GL, Doallo R, Posada D. 2012. jModelTest 2: more models, new heuristics and parallel computing. Nature Methods 9(8): 772. DOI: 10.1038/nmeth.2109

Kass JM, Pinilla-Buitrago GE, Paz A, et al. 2022. Wallace 2: a shiny app for modeling species niches and distributions redesigned to facilitate expansion

Matzke NJ. 2016. Stochastic mapping under biogeographical models. Available at: PhyloWiki BioGeoBEARS. http://phylo.wikidot.com/biogeobears#stochastic_mapping. 11 Nov. 2023 (Date of last successful access)

Rambaut A, Drummond AJ, Xie D, Baele G, Suchard MA. 2018. Posterior summarization in Bayesian phylogenetics using Tracer 1.7. Systematic Biology 67(5): 901–904. DOI: 10.1093/sysbio/syy032

R Core Team. 2020. R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. Available at: https://www.R-project.org/

Tamura K, Stecher G, Kumar S. 2021. MEGA11: Molecular Evolutionary Genetics Analysis Version 11. Molecular Biology and Evolution 38: 3022–3027. DOI: 10.1093/molbev/msab120

Vaidya G, Lohman DJ, Meier R. 2011. SequenceMatrix: concatenation software for the fast assembly of multi‐gene datasets with character set and codon information. Cladistics 27: 171–180. DOI: 10.1111/j.1096-0031.2010.00329.x

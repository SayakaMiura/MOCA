MOCA_v0.1.2  
(Copyright 2021, Authors and Temple University; see license below)

Updated January 28, 2022
==================

The MOCA (Multi Omics Concordance Analysis) tool has been developed by Jared Huzar. It is written in R (version 3.6.2) on Windows platform (64-bit). You are free to download, modify, and expand this code under a permissive license similar to the BSD 2-Clause License (see below). 

MOCA tests concordance between genetic evolution and transcriptomic evolution using tumor single cell sequencing data. Users can easily expand the capabilities by providing additional cell annotation information, such as spatial, metastatic, etc. See Huzar et al. (ref. 1) for further details on the power and functionality of MOCA.   

Dependencies
==================
1. R (v3.6.2 was tested)
 R package: 
    monocle
    ape
    phytools
    ggtree
    stringr
    ggplot2
    apTreeshape
    dplyr
    lsr
    stats
 
 Please make sure Rscript command is functional.

How to use 
==================
1. Download MOCA.zip

2. Load the MOCA functions into your R environment

Open the scripts pertaining to each function, or use the source command to load the functions into your R environment. 

3. Run the desired functions

Scripts can be run either interactively in RStudio, or using the RScript command in terminal/command prompt.

Please see the MOCA manual for more details regarding installation and usage of MOCA and its tools.

Output files
==================
MOCA's functions produce several output files and figures. They will appear in the working directory after the functions are run.

Example
==================
An example dataset is provided (MOCA\Example\) to run MOCA tools. Provided is an expression matrix (MOCA\Example\expression_matrix1.txt), an annotation file (MOCA\Example\alternates\annotate_file1.txt), and a phylogenetic tree (MOCA/Example/tree1.txt). The expression matrix, along with the sequence data which was used to generate the phylogenetic tree, are from Patel et al. From the tree users can produce the annotation file provided, however if they are only interested in the concordance functions the annotation file is provided. Please see the MOCA manual (Manual.pdf) for an example workflow.

Data
==================
The provided phylogenetic tree (MOCA\Example\tree1.nwk) is in standard newick format and has been inferred using BEAM. The example data stems from a glioblastoma tumor dataset, see Patel et al (see ref 3) for more information on the dataset. Any tree in a standard newick format can be used in MOCA. The provided annotation file (MOCA\Example\annotate_file1.txt) is the result of running BalancedAnnotation on the phylogenetic tree and removing cells which had no ancestry designation and/or cells which were excluded from the expression matrix (MOCA\Example\expression_matrix1.txt). The annotation file must contain two columns. The first column should contain the cell IDs which must match with the column names of the expression matrix, and the second column must contain the cells’ ancestry annotations. The row names must also be the cell IDs. For MOCA, the expression matrix must be a data frame with each row corresponding to a gene, and each column corresponding to a cell. The gene IDs must be the row names and the cell IDs must be the column names. Please see the example data (MOCA\Example\) and the MOCA manual (MOCA_Manual.pdf) for further reference on how to format input data.

All data pertaining to the results in the manuscript is available (Data\). There is a folder containing all of the phylogenetic trees (Data\Trees\), a folder with all of the expression matrices (Data\Expression\), and a folder with all of the annotation files (Data\Annotations\). The data comes from two different studies. Data marked with "Hou" comes from Hou et al (ref 7), and data marked with "MGH26" or "MGH31" is from Patel et al (ref 3). From the Patel et al. study there are 2 different tumors, MGH26 and MGH31. For each tumor there are 4 different trees, expression matrices, and annotation files. The different data files stem from reconstruction using different methods, BEAM (ref 4) and SCITE (ref 5), and different base assignment cutoffs (bases with >60% and >70% of cells were retained). Expression matrices, trees, and annotation files are labeled based on their reconstruction method and which tumor they are from, e.g., MGH26BEAM7Annotation.txt. In addition, there is one expression matrix and one annotation file based on CNV annotations which were originally reported in (ref 6). For more information regarding the data processing and reconstruction please see Huzar et al (ref. 1).

==================
Reference:
[1] Huzar, J., Kim H., Kumar S., Miura S. MOCA for integrated analysis of gene expression and genetic variation in single cells. XXX (2021)

[2] Trapnell, C., et al., The dynamics and regulators of cell fate decisions are revealed by pseudotemporal ordering of single cells. Nat Biotechnol, 2014. 32(4): p. 381-386.

[3] Patel, A.P., et al., Single-cell RNA-seq highlights intratumoral heterogeneity in primary glioblastoma. Science, 2014. 344(6190): p. 1396-401.

[4] Miura S, Huuki LA, Buturla T, Vu T, Gomez K, Kumar S. Computational enhancement of single-cell sequences for inferring tumor evolution. Bioinformatics. 2018;34(17):i917-i26.

[5] Jahn K, Kuipers J, Beerenwinkel N. Tree inference for single-cell data. Genome Biol. 2016;17:86.

[6] Serin Harmanci A, Harmanci AO, Zhou X. CaSpER identifies and visualizes CNV events by integrative analysis of single-cell or bulk RNA-sequencing data. Nat Commun. 2020;11(1):89.

[7] Hou, Y., Guo, H., Cao, C., Li, X., Hu, B., Zhu, P., . . . Peng, J. (2016). Single-cell triple omics sequencing reveals genetic, epigenetic, and transcriptomic heterogeneity in hepatocellular carcinomas. Cell Res, 26(3), 304-319. https://doi.org/10.1038/cr.2016.23

--------
Copyright 2021, Authors and Temple University
BSD 3-Clause "New" or "Revised" License, which is a permissive license similar to the BSD 2-Clause License except that that it prohibits others from using the name of the project or its contributors to promote derived products without written consent. 
Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

<p align="center">
    <img align="center" width="25%" src="https://raw.githubusercontent.com/DennisSchmitz/Jovian/assets/images/Jovian%20logo%20notext.svg">
    <h1 align="center">
        Jovian
        <br>
        A user-friendly Viromics toolkit
    </h1>
</p>

<p align="center">
    <a href="https://github.com/DennisSchmitz/Jovian/releases">
        <img alt="Github release" src="https://img.shields.io/github/v/release/DennisSchmitz/Jovian?include_prereleases&style=flat-square">
    </a>
    <a href="https://www.gnu.org/licenses/agpl-3.0">
        <img alt="licence" src="https://img.shields.io/badge/License-AGPL%20v3-blue.svg?style=flat-square">
    </a>
    <a href="https://snakemake.readthedocs.io">
        <img alt="Snakemake Version" src="https://img.shields.io/badge/snakemake-≥5.20.1-brightgreen.svg?style=flat-square">
    </a>
</p>
<p align="center">
    <b>For Citations, please use the following DOI:</b><br>
    <a href="https://doi.org/10.5281/zenodo.3666156">
        <img alt="Zenodo DOI" src="https://zenodo.org/badge/DOI/10.5281/zenodo.3666156.svg">
    </a>
</p>

<p align="center">
    <b>See the documentation:</b><br>
    <a href="https://jovian.rivm-bioinformatics.com">
        <img alt="Jovian Docs" src="https://github.com/florianzwagemaker/jovian-docs/workflows/Docs/badge.svg?branch=master&">
    </a>
    <br>
    <b>Or view an example notebook:</b><br>
    <a href="https://mybinder.org/v2/gh/DennisSchmitz/Jovian_binder/master?filepath=Notebook_report.ipynb">
        <img alt="Launch an example notebook" src="https://mybinder.org/badge_logo.svg">
    </a>
</p>


**IMPORTANT: manuscript is in preparation**

___

<h2>Table of contents</h2>

- [About Jovian](#about-jovian)
  - [Key features of Jovian:](#key-features-of-jovian)
- [Commands](#commands)
- [Features](#features)
  - [General features](#general-features)
  - [Metagenomics specific features](#metagenomics-specific-features)
  - [Reference-alignment specific features](#reference-alignment-specific-features)
  - [Visualizations](#visualizations)
  - [Virus typing](#virus-typing)
  - [Audit trail](#audit-trail)
- [Requirements](#requirements)
- [Installation](#installation)
- [Usage instructions](#usage-instructions)
- [FAQ](#faq)
- [Example Jovian report](#example-jovian-report)
- [Acknowledgements](#acknowledgements)
- [Authors](#authors)

___

## About Jovian  

Jovian is a Public Health toolkit to automatically process raw NGS data from human clinical matrices (faeces, serum, etc.) into clinically relevant information. It has three main components:  

- **Illumina based Metagenomics**:  
  Includes (amongst other features) data quality control, assembly, taxonomic classification, viral typing, and minority variant identification (quasispecies).  
  :memo: Please refer to the [documentation page for the Illumina Metagenomics workflow][ILM_meta_wf] for more information.  

- **Illumina based Reference-alignment**:  
  Includes (amongst other features) data quality control, alignment, SNP identification, and consensus-sequence generation.  
  :exclamation: A reference fasta is required.  
  :memo: Please refer to the [documentation page for the Illumina Reference based workflow][ILM_ref_wf] for more information.  

- **Nanopore based Reference-alignment**:  
  Includes (amongst other features) data quality control, alignment, SNP identification, and consensus-sequence generation.  
  :exclamation: A reference fasta is required.  
  :exclamation: A fasta with primer sequences is required.  
  :memo: Please refer to the [documentation page for the Nanopore Reference based workflow][ONT_ref_wf] for more information.  

### Key features of Jovian:  

- **User-friendliness**:  
  Wetlab personnel can start, configure and interpret results via an interactive web-report. [Click here for an example report][example_report].  
  This makes doing Public Health analyses much more accessible and user-friendly since minimal command-line skills are required.


- **Audit trail**:  
  All pipeline parameters, software versions, database information and runtime statistics are logged. See details [below](#audit-trail).  

- **Portable**:  
  Jovian is easily installed on off-site computer systems and at back-up sister institutes. Allowing results to be generated even when the internal grid-computer is down (speaking from experience).  

<br>

---

<br>

## Commands  

:memo: Please see the full [Command Line Reference on the documentation site][CLI_ref] for a more detailed explanation of each command, including [example commands for starting an analysis][CLI_examples_analysis] or [common usage examples][CLI_examples_cmd].

Here, we have a short list of commands and use cases that are used very frequently.

**Use case 1:**   
Metagenomic analylsis based on Illumina data:  
```
bash jovian illumina-metagenomics -i <INPUT DIRECTORY>
```

**Use case 2:**  
Align Illumina data against a user-provided reference to generate a consensus genome:  
```
bash jovian illumina-reference -i <INPUT DIRECTORY> -ref <REFERENCE FASTA>
```

**Use case 3:**  
Align Nanopore (multiplex) PCR data against a user-provided reference, remove overrepresented primer sequences, and generate a consensus genome:
```
bash jovian nanopore-reference -i <INPUT DIRECTORY> -ref <REFERENCE FASTA> -pr <PRIMER FASTA>
```
use `bash jovian -h` to see a full list of commands applicable to the Jovian version that you're using.  
<br>

---

## Features

:memo: Please refer to our [documentation for the full list of features][featurelist]

### General features

- Data quality control and cleaning.  
  - Including library fragment length analysis, useful for sample preparation QC.  
- Removal of human* data (patient privacy). _*<sup><sub>You can use [whichever reference you would like][faq]. However, Jovian is intended for human clinical samples.</sup></sub>_  
- Removal of PCR-duplicates for Illumina data.  

### Metagenomics specific features

- Assembly of short reads into bigger scaffolds (often full viral genomes).  
- Taxonomic classification:  
  - Every nucleic acid containing biological entity (i.e. not only viruses) is determined up to species level.  
  - Lowest Common Ancestor (LCA) analysis is performed to move ambiguous results up to their last common ancestor, which makes results more robust.  
- Viral typing:
  - Several viral families and genera can be taxonomically labelled at the sub-species level as described [here](#virus-typing).  
- Viral scaffolds are cross-referenced against the Virus-Host interaction database and NCBI host database.  
- Scaffolds are annotated in detail:  
  - Depth of coverage.  
  - GC content.  
  - Open reading frames (ORFs) are predicted.  
  - Minority variants (quasispecies) are identified.  
- Importantly, results of all processes listed above are presented via an [interactive web-report][example_report] including an [audit trail](#audit-trail).  

### Reference-alignment specific features

- All cleaned reads are aligned against the user-provided reference fasta.  
- In the case of Nanopore (multiplex) PCR data, the overrepresented primer sequences are removed.  
- SNPs are called and a consensus genome is generated.  
- Consensus genomes are filtered at the following coverage cut-off thresholds: 1, 5, 10, 30 and 100x.  
- A tabular overview of the breadth of coverage (BoC) at the different coverage cut-off thresholds is generated.  
- Alignments and visualized via `IGVjs` and allow manual assessment and validation of consensus genomes.  

### Visualizations

All data are visualized via an interactive web-report, [as shown here][example_report], which includes:  

- A collation of interactive QC graphs via `MultiQC`.  
- Taxonomic results are presented on three levels:  
  - For an entire (multi sample) run, interactive heatmaps are made for non-phage viruses, phages and bacteria. They are stratified to different taxonomic levels.  
  - For a sample level overview, `Krona` interactive taxonomic piecharts are generated.  
  - For more detailed analyses, interactive tables are included. Similar to popular spreadsheet applications (e.g. Microsoft Excel).  
    - Classified scaffolds  
    - Unclassified scaffolds (i.e. "Dark Matter")  
- Virus typing results are presented via interactive spreadsheet-like tables.  
- An interactive scaffold alignment viewer (`IGVjs`) is included, containing:  
  - Detailed alignment information.  
  - Depth of coverage graph.  
  - GC content graph.
  - Predicted open reading frames (ORFs).  
  - Identified minority variants (quasispecies).  
- All SNP metrics are presented via interactive spreadsheet-like tables, allowing detailed analysis.  

### Virus typing

After a Jovian analysis is finished you can perform virus-typing (i.e. sub-species level taxonomic labelling). These analyses can be started by the command `bash jovian -vt [virus keyword]`, where `[virus keyword]` can be:  

Keyword | Taxon used for scaffold selection | Notable virus species
--------|-----------------------------------|----------------------
`NoV` | Caliciviridae   | Norovirus GI and GII, Sapovirus  
`EV`  | Picornaviridae  | Enteroviruses (Coxsackie, Polio, Rhino, etc.), Parecho, Aichi, Hepatitis A 
`RVA` | _Rotavirus A_   | Rotavirus A  
`HAV` | _Hepatovirus A_ | Hepatitis A  
`HEV` | _Orthohepevirus A_ | Hepatitis E  
`PV`  | Papillomaviridae | Human Papillomavirus  
`Flavi` | Flaviviridae    | Dengue (work in progress)
`all` | All of the above | All of the above
  
### Audit trail

An audit trail, used for clinical reproducibility and logging, is generated and contains:  

- A unique methodological fingerprint: allowing to exactly reproduce the analysis, even retrospectively by reverting to old versions of the pipeline code.  
- The following information is also logged:  
  - Database timestamps  
  - (user-specified) Pipeline parameters  

However, it has limitations since several things are out-of-scope for Jovian to control:

- The `virus typing-tools` version  
  - Currently we depend on a public web-tool hosted by the [RIVM][RIVM]. These are developed in close collaboration with - *but independently of* - Jovian. A versioning system for the `virus typing-tools` is being worked on, however, this is not trivial and will take some time.  
- Input files and metadata
  - We only save the names and location of input files at the time the analysis was performed. Long-term storage of the data, and documenting their location over time, is the responsibility of the end-user. Likewise, the end-user is responsible for storing datasets with their correct metadata (e.g. clinical information, database versions, etc.). We recommend using [iRODS][irods] for this as described by [Nieroda et al. 2019][irods_ref]. While we acknowledge that database versions are vital to replicate results, the databases Jovian uses have no official versioning, hence why we include timestamps only.  

___

<details>
  <summary>Jovian Illumina Metagenomics workflow visualization</summary>
  Click the image for a full-sized version  
  <a href="https://raw.githubusercontent.com/DennisSchmitz/Jovian/assets/images/Jovian_Workflow_Illumina_meta.svg">
    <img  alt="Jovian Illumina Metagenomics workflow" src="https://raw.githubusercontent.com/DennisSchmitz/Jovian/assets/images/Jovian_Workflow_Illumina_meta.svg">
  </a>
</details>

<br>

<details>
  <summary>Jovian Illumina Reference alignment workflow visualization</summary>
  Click the image for a full-sized version  
  <a href="https://raw.githubusercontent.com/DennisSchmitz/Jovian/assets/images/Jovian_Workflow_Illumina_ref.svg">
    <img alt="Jovian Illumina Reference workflow" src="https://raw.githubusercontent.com/DennisSchmitz/Jovian/assets/images/Jovian_Workflow_Illumina_ref.svg">
  </a>
</details>

<br>

<details>
  <summary>Jovian Nanopore Reference alignment workflow visualization</summary>
  Click the image for a full-sized version  
  <a href="https://raw.githubusercontent.com/DennisSchmitz/Jovian/assets/images/Jovian_Workflow_Nanopore_ref.svg">
    <img alt="Jovian Nanopore reference workflow" src="https://raw.githubusercontent.com/DennisSchmitz/Jovian/assets/images/Jovian_Workflow_Nanopore_ref.svg">
  </a>
</details>

___

## Requirements

:memo: Please refer to our [documentation][requirements] for a detailed overview of the Jovian requirements [here][requirements]
___

## Installation

:memo: Please refer to our [documentation][Install_instructions] for detailed instructions regarding the installation of Jovian [here][Install_instructions]


## Usage instructions

General usage instructions vary for each workflow that we support.  
Please refer to the link below corresponding to the workflow that you wish to use

- [Illumina Metagenomics workflow - Usage instructions and information][ILM_meta_wf]
- [Illumina Reference alignment workflow - Usage instructions and information][ILM_ref_wf]
- [Nanopore Reference alignment workflow - Usage instructions and information][ONT_ref_wf]

___

## FAQ

Can be found [here][faq].  
___

## Example Jovian report

Can be found [here][example_report].  
___

## Acknowledgements

|Name |Publication|Website|
|:---|:---|:---|
|BBtools|NA|https://jgi.doe.gov/data-and-tools/bbtools/|
|BEDtools|Quinlan, A.R. and I.M.J.B. Hall, BEDTools: a flexible suite of utilities for comparing genomic features. 2010. 26(6): p. 841-842.|https://bedtools.readthedocs.io/en/latest/|
|BLAST|Altschul, S.F., et al., Gapped BLAST and PSI-BLAST: a new generation of protein database search programs. 1997. 25(17): p. 3389-3402.|https://www.ncbi.nlm.nih.gov/books/NBK279690/|
|BWA|Li, H. (2013). Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM. arXiv preprint arXiv:1303.3997.|https://github.com/lh3/bwa|
|BioConda|Grüning, B., et al., Bioconda: sustainable and comprehensive software distribution for the life sciences. 2018. 15(7): p. 475.|https://bioconda.github.io/|
|Biopython|Cock, P. J., Antao, T., Chang, J. T., Chapman, B. A., Cox, C. J., Dalke, A., ... & De Hoon, M. J. (2009). Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics, 25(11), 1422-1423.|https://biopython.org/|
|Bokeh|Bokeh Development Team (2018). Bokeh: Python library for interactive visualization.|https://bokeh.pydata.org/en/latest/|
|Bowtie2|Langmead, B. and S.L.J.N.m. Salzberg, Fast gapped-read alignment with Bowtie 2. 2012. 9(4): p. 357.|http://bowtie-bio.sourceforge.net/bowtie2/index.shtml|
|Conda|NA|https://conda.io/|
|DRMAA|NA|http://drmaa-python.github.io/|
|FastQC|Andrews, S., FastQC: a quality control tool for high throughput sequence data. 2010.|https://www.bioinformatics.babraham.ac.uk/projects/fastqc/|
|gawk|NA|https://www.gnu.org/software/gawk/|
|GNU Parallel|O. Tange (2018): GNU Parallel 2018, March 2018, https://doi.org/10.5281/zenodo.1146014.|https://www.gnu.org/software/parallel/|
|Git|NA|https://git-scm.com/|
|igvtools|NA|https://software.broadinstitute.org/software/igv/igvtools|
|Jupyter Notebook|Kluyver, Thomas, et al. "Jupyter Notebooks-a publishing format for reproducible computational workflows." ELPUB. 2016.|https://jupyter.org/|
|Jupyter_contrib_nbextension|NA|https://github.com/ipython-contrib/jupyter_contrib_nbextensions|
|Jupyterthemes|NA|https://github.com/dunovank/jupyter-themes|
|Krona|Ondov, B.D., N.H. Bergman, and A.M. Phillippy, Interactive metagenomic visualization in a Web browser. BMC Bioinformatics, 2011. 12: p. 385.|https://github.com/marbl/Krona/wiki|
|Lofreq|Wilm, A., et al., LoFreq: a sequence-quality aware, ultra-sensitive variant caller for uncovering cell-population heterogeneity from high-throughput sequencing datasets. 2012. 40(22): p. 11189-11201.|http://csb5.github.io/lofreq/|
|MGkit|Rubino, F. and Creevey, C.J. 2014. MGkit: Metagenomic Framework For The Study Of Microbial Communities. . Available at: figshare [doi:10.6084/m9.figshare.1269288].|https://bitbucket.org/setsuna80/mgkit/src/develop/|
|Minimap2|Li, H., Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics, 2018.|https://github.com/lh3/minimap2|
|MultiQC|Ewels, P., et al., MultiQC: summarize analysis results for multiple tools and samples in a single report. 2016. 32(19): p. 3047-3048.|https://multiqc.info/|
|Nb_conda|NA|https://github.com/Anaconda-Platform/nb_conda|
|Nb_conda_kernels|NA|https://github.com/Anaconda-Platform/nb_conda_kernels|
|Nginx|NA|https://www.nginx.com/|
|Numpy|Walt, S. V. D., Colbert, S. C., & Varoquaux, G. (2011). The NumPy array: a structure for efficient numerical computation. Computing in Science & Engineering, 13(2), 22-30.|http://www.numpy.org/|
|Pandas|McKinney, W. Data structures for statistical computing in python. in Proceedings of the 9th Python in Science Conference. 2010. Austin, TX.|https://pandas.pydata.org/|
|Picard|NA|https://broadinstitute.github.io/picard/|
|Prodigal|Hyatt, D., et al., Prodigal: prokaryotic gene recognition and translation initiation site identification. 2010. 11(1): p. 119.|https://github.com/hyattpd/Prodigal/wiki/Introduction|
|Python|G. van Rossum, Python tutorial, Technical Report CS-R9526, Centrum voor Wiskunde en Informatica (CWI), Amsterdam, May 1995.|https://www.python.org/|
|Qgrid|NA|https://github.com/quantopian/qgrid|
|SAMtools|Li, H., et al., The sequence alignment/map format and SAMtools. 2009. 25(16): p. 2078-2079.|http://www.htslib.org/|
|SPAdes|Nurk, S., et al., metaSPAdes: a new versatile metagenomic assembler. Genome Res, 2017. 27(5): p. 824-834.|http://cab.spbu.ru/software/spades/|
|seqkit|Shen, Wei, et al. "SeqKit: a cross-platform and ultrafast toolkit for FASTA/Q file manipulation." PloS one 11.10 (2016).|https://github.com/shenwei356/seqkit|
|Seqtk|NA|https://github.com/lh3/seqtk|
|Snakemake|Köster, J. and S.J.B. Rahmann, Snakemake—a scalable bioinformatics workflow engine. 2012. 28(19): p. 2520-2522.|https://snakemake.readthedocs.io/en/stable/|
|Tabix|NA|www.htslib.org/doc/tabix.html|
|tree|NA|http://mama.indstate.edu/users/ice/tree/|
|Trimmomatic|Bolger, A.M., M. Lohse, and B. Usadel, Trimmomatic: a flexible trimmer for Illumina sequence data. Bioinformatics, 2014. 30(15): p. 2114-20.|www.usadellab.org/cms/?page=trimmomatic|
|Virus-Host Database|Mihara, T., Nishimura, Y., Shimizu, Y., Nishiyama, H., Yoshikawa, G., Uehara, H., ... & Ogata, H. (2016). Linking virus genomes with host taxonomy. Viruses, 8(3), 66.|http://www.genome.jp/virushostdb/note.html|
|Virus typing tools|Kroneman, A., Vennema, H., Deforche, K., Avoort, H. V. D., Penaranda, S., Oberste, M. S., ... & Koopmans, M. (2011). An automated genotyping tool for enteroviruses and noroviruses. Journal of Clinical Virology, 51(2), 121-125.|https://www.ncbi.nlm.nih.gov/pubmed/21514213|

## Authors

- Dennis Schmitz ([RIVM][RIVM] and [EMC][EMC])  
- Sam Nooij ([RIVM][RIVM] and [EMC][EMC])  
- Robert Verhagen ([RIVM][RIVM])  
- Thierry Janssens ([RIVM][RIVM])  
- Jeroen Cremer ([RIVM][RIVM])  
- Florian Zwagemaker ([RIVM][RIVM])  
- Mark Kroon ([RIVM][RIVM])  
- Erwin van Wieringen ([RIVM][RIVM])  
- Harry Vennema ([RIVM][RIVM])  
- Annelies Kroneman ([RIVM][RIVM])  
- Marion Koopmans ([EMC][EMC])  

____
_This project/research has received funding from the European Union’s Horizon 2020 research and innovation programme under grant agreement No. 643476. and the Dutch working group on molecular diagnostics (WMDI)._
____

[ILM_meta_wf]: https://jovian.rivm-bioinformatics.com/Workflows/Illumina/Metagenomics/
[ILM_ref_wf]: https://jovian.rivm-bioinformatics.com/Workflows/Illumina/Reference/
[ONT_ref_wf]: https://jovian.rivm-bioinformatics.com/Workflows/Nanopore/Reference/

[example_report]: https://mybinder.org/v2/gh/DennisSchmitz/Jovian_binder/master?filepath=Notebook_report.ipynb
[CLI_ref]: https://jovian.rivm-bioinformatics.com/CLI/
[CLI_examples_analysis]: https://jovian.rivm-bioinformatics.com/CLI/#examples
[CLI_examples_cmd]: https://jovian.rivm-bioinformatics.com/CLI/#examples_1


[featurelist]: https://jovian.rivm-bioinformatics.com/#features
[faq]: https://jovian.rivm-bioinformatics.com/FAQ/

[requirements]: https://jovian.rivm-bioinformatics.com/Getting-started/Requirements/
[Install_instructions]: https://jovian.rivm-bioinformatics.com/Getting-started/Installation/

[RIVM]: https://www.rivm.nl/en
[EMC]: https://www6.erasmusmc.nl/viroscience/

[irods]: https://irods.org
[irods_ref]: https://www.ncbi.nlm.nih.gov/pubmed/30646845

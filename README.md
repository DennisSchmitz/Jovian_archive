# Jovian, user-friendly metagenomics     

**IMPORTANT: Do not share the code without my express permission as it is unpublished (manuscript in preparation)**  

___

<img align="right" src="../assets/images/Jovian_logo.png">

## Table of content  
- [Table of content](#table-of-content)
- [Jovian description](#jovian-description)
  - [Features](#features)
  - [Visualizations](#visualizations)
  - [Virus typing](#virus-typing)
  - [Audit trail](#audit-trail)
- [Requirements](#requirements)
  - [System requirements](#system-requirements)
  - [Software](#software)
  - [Databases](#databases)
- [Instructions](#instructions)
  - [Installation](#installation)
  - [How to start/configure a Jovian analysis](#how-to-startconfigure-a-jovian-analysis)
  - [Explanation of output folders](#explanation-of-output-folders)
- [FAQ](#faq)
- [Example Jovian report](#example-jovian-report)
- [Acknowledgements](#acknowledgements)
    - [Authors](#authors)

___

## Jovian description  

The pipeline automatically processes raw Illumina NGS data from human clinical matrices (faeces, serum, etc.) into clinically relevant information such as taxonomic classification, viral typing and minority variant identification (quasispecies).
Wetlab personnel can start, configure and interpret results via an interactive web-report. This makes doing metagenomics analyses much more accessible and user-friendly since minimal command-line skills are required.  

### Features    
- Data quality control (QC) and cleaning.  
  - Including library fragment length analysis, usefull for sample preparation QC.  
- Removal of human* data (patient privacy). _*<sup><sub>You can use [whichever reference you would like](../../wiki/Frequently-Asked-Questions#i-dont-care-about-removing-the-human-data-i-have-samples-that-are-from-other-species-can-i-also-automatically-remove-that). However, Jovian is intended for human clinical samples.</sup></sub>_  
- Assembly of short reads into bigger scaffolds (often full viral genomes).  
- Taxonomic classification:  
  - Every nucleic acid containing biological entity (i.e. not only viruses) is determined up to species level.  
  - Lowest Common Ancestor (LCA) analysis is performed to move ambiguous results up to their last common ancestor, which makes results more robust.  
- Viral typing:
  - Several viral families and genera can be taxonomically labelled at the sub-species level as described [here](#virus-typing).  
- Viral scaffolds are cross-referenced against the Virus-Host interaction database and NCBI host database.  
- Scaffolds are annotated with great detail:  
  - Depth of coverage.  
  - GC content.  
  - Open reading frames (ORFs) are predicted.  
  - Minority variants (quasispecies) are identified.  
- Importantly, results of all processes listed above are presented via an [interactive web-report](#visualizations) including an [audit trail](#audit-trail).  

### Visualizations  
All data are visualized via an interactive web-report, [as shown here](#example-jovian-report), which includes:  
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
`NoV` | Caliciviridae   | #TODO  
`EV`  | Picornaviridae  | #TODO  
`RVA` | _Rotavirus A_   | #TODO  
`HAV` | _Hepatovirus A_ (Hepatitis A) | #TODO  
`HEV` | _Orthohepevirus A_ (Hepatitis E)| #TODO  
`PV`  | Papillomaviridae | #TODO  
`Flavi` | Flaviviridae    | #TODO  
  
### Audit trail  
An audit trail, used for clinical reproducability and logging, is generated and contains:  
- A unique methodological fingerprint of the code is generated and accessible via GitHub: allowing to exactly reproduce the analysis, even retrospectively by reverting to old versions of the pipeline code.  
- The following information is also logged:  
  - Database timestamps  
  - (user-specified) Pipeline parameters  

However, several things are out-of-scope for Jovian logging:
- The `IGVjs` version  
- The `virus typing-tools` version  
  - Currently we depend on a public web-tool hosted by the [RIVM](https://www.rivm.nl/en). These are developed in close collaboration with - *but independently of* - Jovian. A versioning system for the `virus typing-tools` is being worked on, however, this is not trivial and will take some time.  
- The database versions  
  - We only save the timestamps of the database files, this is because the databases used by Jovian have no official versioning. Any versioning scheme is therefore out-of-scope for Jovian and a responsibility of the end-user.  
- Input files
  - We only save the names and location of input files at the time the analysis was performed. Long-term storage of the data, and documenting their location over time, is the responsibility of the end-user.  

___

![Jovian_rulegraph.png](../assets/images/rulegraph_Jovian.png?raw=true)
___

## Requirements
Jovian has minimal [system requirements](#system-requirements), it is intended for powerful servers and/or grid-computers but also works on powerful PCs and laptops. It has two major software dependencies, miniConda and IGVjs. During the installation procedure you will be asked if you want to automatically install these (if Conda is not already available). Additionally, it depends on the [following software](#software), but most Linux systems have these pre-installed already. Any metagenomics analysis, Jovian included, depends on several [public databases that you have to download](#databases).    

### System requirements
We have developed and tested the software for the following Linux distributions ("distro's"): `RHEL`, `CentOS` and `Ubuntu`. We do expect Jovian to work on other Linux distro's, but cannot guarantee stability. We are currently assessing if it feasible to also test our software on Windows computers using [Ubuntu on Windows](https://www.microsoft.com/en-us/p/ubuntu/9nblggh4msv6?activetab=pivot:overviewtab). This is a work-in-progress, but preliminary tests look promising.  

The software does <u><b>not</b></u> require root and/or sudo-rights. It is however necessary to have read and write access to the `/tmp` folder on your system. This won't be a problem most of the time since the `/tmp` folder is usually free to read from and write to. However, it is best to check this with your system administrator(s).  

The databases require up to 400GB of disk-space. The pipeline itself, including all Conda software environments and IGVjs, requires up to 15GB of disk-space.

### Software  
|Software name|Website|  
|:---|:---|  
|`git`| https://git-scm.com/downloads |  
|`curl`| https://curl.haxx.se/ |  
|`which`| http://savannah.gnu.org/projects/which |  
|`bzip2`| http://www.bzip.org/ |  

### Databases  
|Database name|Link|Installation instructions|
|:---|:---|:---|
|`NCBI NT & NR`| ftp://ftp.ncbi.nlm.nih.gov/blast/db/ | [link](../../wiki/Installation-Instructions#database-installation)|
|`NCBI Taxdump`| ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/ | [link](../../wiki/Installation-Instructions#database-installation)|
|`NCBI New_taxdump`| ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/ | [link](../../wiki/Installation-Instructions#database-installation)|
|`Virus-Host interaction database`| http://www.genome.jp/virushostdb/note.html | [link](../../wiki/Installation-Instructions#database-installation)|
|`Latest Human Genome`*| https://support.illumina.com/sequencing/sequencing_software/igenome.html | [link](../../wiki/Installation-Instructions#database-installation)|

_* We suggest the latest human genome because Jovian is intended for clinical samples. You can however use any reference you'd like, as [explained here](../../wiki/Frequently-Asked-Questions#i-dont-care-about-removing-the-human-data-i-have-samples-that-are-from-other-species-can-i-also-automatically-remove-that)._

___

## Instructions 
Below you'll find instructions on how to install and start/configure a Jovian analysis.   

### Installation
Can be found on [this wiki page](../../wiki/Installation-Instructions).

### How to start/configure a Jovian analysis  
Currently, the method to launch analyses via the Jupyter Notebook requires some minor tweaks. So I cannot share it yet, we recommend you to use the command-line method below.  

<b>Jupyter Notebook method:</b>  
- Via your Jupyter Notebook browser connection, go to the `Jovian` folder [created during installation](../../wiki/Installation-Instructions). Then, open `Notebook_portal.ipynb`.  
- Follow the instructions in this notebook to start an analysis.  

<b>Command-line interface method:</b>  
- Make sure that Jovian is completely [installed](../../wiki/Installation-Instructions).  
- Go to the folder where `Jovian` was installed.  
- Configure pipeline parameters by changing the [profile/pipeline_parameters.yaml](profile/pipeline_parameters.yaml) file. Either via Jupyter Notebook or with a commandline text-editor of choice.  
- Optional: We recommended you do a `dry-run` before each analysis to check if there are any typo's, missing files or other errors. This can be done via `bash jovian -i <input_directory> -n`
- If the dry-run has completed without errors, you are ready to start a real analysis with the following command:  
`bash jovian -i <input_directory>` 
- After the pipeline has finished, open `Notebook_report.ipynb` via your browser. Click on `Cell` in the toolbar, then press `Run all` and wait for data to be imported.  
  - N.B. You need to have a Jupyter notebook process running in the background, as described [here](../../wiki/Installation-Instructions#starting-the-jupyter-notebook-server-process).

### Explanation of output folders  
|Folder|Contents|
|:---|:---|
|`bin/` |Contains the scripts required for Jovian to work |
|`data/` | Contains intermediate and detailed data |
|`envs/` | Contains all conda environments for the pipeline |
|`files/` | Contains ancillary files for the pipeline |
|`logs/` | Contains all Jovian log files, use these files to troubleshoot errors |
|`profile/` | Contains the files with Snakemake and pipeline parameters |
|`results/` | This contains all files that are important for end-users and are imported by the Jupyter Report |

Also, a hidden folder named `.snakemake` is generated. Do not remove or edit this folder. Jovian was built via `Snakemake` and this folder contains all the software and file-metadata required for proper pipeline functionality.  
___

## FAQ
Can be found on [this wiki page](../../wiki/Frequently-Asked-Questions).
___

## Example Jovian report  
_Data shown below is based on public data, available on ENA via accession ID `PRJNA491626`. It contains Illumina paired-end data of faeces from people with gastroenteritis._  

**MultiQC is used to summarize many pipeline metrics, including read quality, insert-size distribution, biases, etc.:**  
<br>
![Jovian_QC-report_PRJNA491626-public-dataset.PNG](../assets/images/screenshots/Jovian_QC-report_PRJNA491626-public-dataset.PNG?raw=true)

**A summary barchart overview of the entire dataset is also presented:**
<br>
![Jovian_barcharts_PRJNA491626-public-dataset.PNG](../assets/images/screenshots/Jovian_barcharts_PRJNA491626-public-dataset.PNG?raw=true)

**Metagenomics data is presented through three different visualizations, Krona pie-charts give sample level overview:**  
<br>
![Jovian_Krona-chart_PRJNA491626-public-dataset.PNG](../assets/images/screenshots/Jovian_Krona-chart_PRJNA491626-public-dataset.PNG?raw=true)

**Viral and bacterial heatmaps that are stratified to different taxonomic levels give an overview of the complete dataset:**  
<br>
![Jovian_heatmap_PRJNA491626-public-dataset.PNG](../assets/images/screenshots/Jovian_heatmap_PRJNA491626-public-dataset.png?raw=true)

**All classified scaffolds and their metrics are presented through interactive tables that are functionally similar to popular spreadsheet programs. Allowing filtering for certain metrics, e.g. taxonomic level (species up to superkingdom), length, number of ORFs, percentage GC, depth of coverage, etc. to facilitate in-depth analyses:**
<br>
![Jovian_classified-scaffolds_PRJNA491626-public-dataset.png](../assets/images/screenshots/Jovian_classified-scaffolds_PRJNA491626-public-dataset.png?raw=true)

**Any scaffold that could not be classified ("dark matter") is reported in a similar interactive table for further investigation:**
<br>
![Jovian_dark-matter_PRJNA491626-public-dataset.PNG](../assets/images/screenshots/Jovian_dark-matter_PRJNA491626-public-dataset.PNG?raw=true)

**All classified scaffolds are also cross-referenced against the NCBI host information and the Virus-Host database:**
<br>
![Jovian_host-disease_PRJNA491626-public-dataset.PNG](../assets/images/screenshots/Jovian_host-disease_PRJNA491626-public-dataset.PNG?raw=true)

**The typing-tool output for Caliciviridae, Picornaviridae, Hepatoviruses, Orthohepeviruses and Rotaviruses containing the genotype information are also presented:**
<br>
![Jovian_NoV-TT_PRJNA491626-public-dataset.PNG](../assets/images/screenshots/Jovian_NoV-TT_PRJNA491626-public-dataset.PNG?raw=true)
![Jovian_EV-TT_PRJNA491626-public-dataset.PNG](../assets/images/screenshots/Jovian_EV-TT_PRJNA491626-public-dataset.PNG?raw=true)
![Jovian_HAV-TT_PRJNA491626-public-dataset.PNG](../assets/images/screenshots/Jovian_HAV-TT_PRJNA491626-public-dataset.PNG?raw=true)
![Jovian_RVA-TT_PRJNA491626-public-dataset.PNG](../assets/images/screenshots/Jovian_RVA-TT_PRJNA491626-public-dataset.PNG?raw=true)

**IGVjs is used to visualize genomes, you can zoom in to individual sites to inspect e.g. minority variants in greater detail. It incorporates and shows the depth of coverage, GC contents, predicted ORFs, minority variants (quasispecies) alongside each individual aligning read:**  
<br>
![Jovian_IGVjs_PRJNA491626-public-dataset.PNG](../assets/images/screenshots/Jovian_IGVjs_PRJNA491626-public-dataset.PNG?raw=true)
![Jovian_IGVjs-zoom_PRJNA491626-public-dataset.PNG](../assets/images/screenshots/Jovian_IGVjs-zoom_PRJNA491626-public-dataset.PNG?raw=true)

**The SNP information is also presented through a spreadsheet table for filtering and in-depth analysis:**  
<br>
![Jovian_minority-SNP-table_PRJNA491626-public-dataset.PNG](../assets/images/screenshots/Jovian_minority-SNP-table_PRJNA491626-public-dataset.PNG?raw=true)

**Lastly, the logging of software, databases and pipeline settings are presented to the user. A verbose list containing all software in the current running environment, `Jovian_master`, is reported (not shown). Also, a list containing the timestamps of all used databases are reported (not shown). Via Snakemake a report is created describing exactly what software and which versions were used (shown below), alongside information about how long each step in the pipeline took to complete (not shown). The Git hash is reported, the unique Jovian methodological "fingerprint", which allows exact reproduction of results at a later time (shown below). And pipeline settings for the current analysis are reported (shown below):**  
<br>
![Jovian_Snakemake-report_PRJNA491626-public-dataset.PNG](../assets/images/screenshots/Jovian_Snakemake-report_PRJNA491626-public-dataset.PNG?raw=true)
![Jovian_logging-git-hash-config_PRJNA491626-public-dataset.PNG](../assets/images/screenshots/Jovian_logging-git-hash-config_PRJNA491626-public-dataset.PNG?raw=true)

___

## Acknowledgements

|Name |Publication|Website|
|:---|:---|:---|
|`BBtools`|NA|https://jgi.doe.gov/data-and-tools/bbtools/|
|`BEDtools`|Quinlan, A.R. and I.M.J.B. Hall, BEDTools: a flexible suite of utilities for comparing genomic features. 2010. 26(6): p. 841-842.|https://bedtools.readthedocs.io/en/latest/|
|`BLAST`|Altschul, S.F., et al., Gapped BLAST and PSI-BLAST: a new generation of protein database search programs. 1997. 25(17): p. 3389-3402.|https://www.ncbi.nlm.nih.gov/books/NBK279690/|
|`BWA`|Li, H. (2013). Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM. arXiv preprint arXiv:1303.3997.|https://github.com/lh3/bwa|
|`BioConda`|Grüning, B., et al., Bioconda: sustainable and comprehensive software distribution for the life sciences. 2018. 15(7): p. 475.|https://bioconda.github.io/|
|`Biopython`|Cock, P. J., Antao, T., Chang, J. T., Chapman, B. A., Cox, C. J., Dalke, A., ... & De Hoon, M. J. (2009). Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics, 25(11), 1422-1423.|https://biopython.org/|
|`Bokeh`|Bokeh Development Team (2018). Bokeh: Python library for interactive visualization.|https://bokeh.pydata.org/en/latest/|
|`Bowtie2`|Langmead, B. and S.L.J.N.m. Salzberg, Fast gapped-read alignment with Bowtie 2. 2012. 9(4): p. 357.|http://bowtie-bio.sourceforge.net/bowtie2/index.shtml|
|`Conda`|NA|https://conda.io/|
|`DRMAA`|NA|http://drmaa-python.github.io/|
|`FastQC`|Andrews, S., FastQC: a quality control tool for high throughput sequence data. 2010.|https://www.bioinformatics.babraham.ac.uk/projects/fastqc/|
|`gawk`|NA|https://www.gnu.org/software/gawk/|
|`GNU Parallel`|O. Tange (2018): GNU Parallel 2018, March 2018, https://doi.org/10.5281/zenodo.1146014.|https://www.gnu.org/software/parallel/|
|`Git`|NA|https://git-scm.com/|
|`igvtools`|NA|https://software.broadinstitute.org/software/igv/igvtools|
|`Jupyter Notebook`|Kluyver, Thomas, et al. "Jupyter Notebooks-a publishing format for reproducible computational workflows." ELPUB. 2016.|https://jupyter.org/|
|`Jupyter_contrib_nbextension`|NA|https://github.com/ipython-contrib/jupyter_contrib_nbextensions|
|`Jupyterthemes`|NA|https://github.com/dunovank/jupyter-themes|
|`Krona`|Ondov, B.D., N.H. Bergman, and A.M. Phillippy, Interactive metagenomic visualization in a Web browser. BMC Bioinformatics, 2011. 12: p. 385.|https://github.com/marbl/Krona/wiki|
|`Lofreq`|Wilm, A., et al., LoFreq: a sequence-quality aware, ultra-sensitive variant caller for uncovering cell-population heterogeneity from high-throughput sequencing datasets. 2012. 40(22): p. 11189-11201.|http://csb5.github.io/lofreq/|
|`Minimap2`|Li, H., Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics, 2018.|https://github.com/lh3/minimap2|
|`MultiQC`|Ewels, P., et al., MultiQC: summarize analysis results for multiple tools and samples in a single report. 2016. 32(19): p. 3047-3048.|https://multiqc.info/|
|`Nb_conda`|NA|https://github.com/Anaconda-Platform/nb_conda|
|`Nb_conda_kernels`|NA|https://github.com/Anaconda-Platform/nb_conda_kernels|
|`Nginx`|NA|https://www.nginx.com/|
|`Numpy`|Walt, S. V. D., Colbert, S. C., & Varoquaux, G. (2011). The NumPy array: a structure for efficient numerical computation. Computing in Science & Engineering, 13(2), 22-30.|http://www.numpy.org/|
|`Pandas`|McKinney, W. Data structures for statistical computing in python. in Proceedings of the 9th Python in Science Conference. 2010. Austin, TX.|https://pandas.pydata.org/|
|`Picard`|NA|https://broadinstitute.github.io/picard/|
|`Prodigal`|Hyatt, D., et al., Prodigal: prokaryotic gene recognition and translation initiation site identification. 2010. 11(1): p. 119.|https://github.com/hyattpd/Prodigal/wiki/Introduction|
|`Python`|G. van Rossum, Python tutorial, Technical Report CS-R9526, Centrum voor Wiskunde en Informatica (CWI), Amsterdam, May 1995.|https://www.python.org/|
|`Qgrid`|NA|https://github.com/quantopian/qgrid|
|`SAMtools`|Li, H., et al., The sequence alignment/map format and SAMtools. 2009. 25(16): p. 2078-2079.|http://www.htslib.org/|
|`SPAdes`|Nurk, S., et al., metaSPAdes: a new versatile metagenomic assembler. Genome Res, 2017. 27(5): p. 824-834.|http://cab.spbu.ru/software/spades/|
|`Seqtk`|NA|https://github.com/lh3/seqtk|
|`Snakemake`|Köster, J. and S.J.B. Rahmann, Snakemake—a scalable bioinformatics workflow engine. 2012. 28(19): p. 2520-2522.|https://snakemake.readthedocs.io/en/stable/|
|`Tabix`|NA|www.htslib.org/doc/tabix.html|
|`tree`|NA|http://mama.indstate.edu/users/ice/tree/|
|`Trimmomatic`|Bolger, A.M., M. Lohse, and B. Usadel, Trimmomatic: a flexible trimmer for Illumina sequence data. Bioinformatics, 2014. 30(15): p. 2114-20.|www.usadellab.org/cms/?page=trimmomatic|
|`Virus-Host Database`|Mihara, T., Nishimura, Y., Shimizu, Y., Nishiyama, H., Yoshikawa, G., Uehara, H., ... & Ogata, H. (2016). Linking virus genomes with host taxonomy. Viruses, 8(3), 66.|http://www.genome.jp/virushostdb/note.html|
|`Virus typing-tools`|Kroneman, A., Vennema, H., Deforche, K., Avoort, H. V. D., Penaranda, S., Oberste, M. S., ... & Koopmans, M. (2011). An automated genotyping tool for enteroviruses and noroviruses. Journal of Clinical Virology, 51(2), 121-125.|https://www.ncbi.nlm.nih.gov/pubmed/21514213|

#### Authors
- Dennis Schmitz ([RIVM](https://www.rivm.nl/en) and [EMC](https://www6.erasmusmc.nl/viroscience/))  
- Sam Nooij ([RIVM](https://www.rivm.nl/en) and [EMC](https://www6.erasmusmc.nl/viroscience/))  
- Robert Verhagen ([RIVM](https://www.rivm.nl/en))  
- Thierry Janssens ([RIVM](https://www.rivm.nl/en))  
- Jeroen Cremer ([RIVM](https://www.rivm.nl/en))  
- Florian Zwagemaker ([RIVM](https://www.rivm.nl/en))  
- Mark Kroon ([RIVM](https://www.rivm.nl/en))  
- Erwin van Wieringen ([RIVM](https://www.rivm.nl/en))  
- Harry Vennema ([RIVM](https://www.rivm.nl/en))  
- Annelies Kroneman ([RIVM](https://www.rivm.nl/en))  
- Marion Koopmans ([EMC](https://www6.erasmusmc.nl/viroscience/))  

____
_This project/research has received funding from the European Union’s Horizon 2020 research and innovation programme under grant agreement No. 643476._
____

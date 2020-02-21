.. Jovian documentation master file
=====================================
Jovian, user-friendly metagenomics
=====================================

.. image:: https://github.com/DennisSchmitz/Jovian/raw/assets/images/Jovian_logo.png
   :align: right
   :width: 275px

Jovian is an extensive bioinformatics pipeline which automatically processes raw Illumina NGS data into clinically relevant information.

Some of the included features are taxonomic classification, viral typing and minority variant identification.

Wetlab personnel (as well as experienced bioinformaticians) can start, configure and interpret results via an interactive web-report.
This makes doing metagenomics analyses much more accessible and user-friendly since minimal command-line skills are required.


Features
================

* Data quality control and cleaning
* Removal of background host data
* Assembly of short reads into bigger scaffolds (often whole viral genomes)
* Taxonomic classification
* Viral typing
* Cross referencing against the Virus_Host interaction database and the NCBI host database.
* Detailed annotation of scaffolds
* Full audit trail
* Easy presentation of results through an interactive web report

Visualizations
======================

All data is visualized through an interactive web-report powered by Jupyter-notebooks. An example can be found here_.

* A collection of interactive Quality-control graphs via ``MultiQC``
* Taxonomic results are presented on three levels:

   1. Interactive heatmaps are generated for virusses, phages and bacteria.
   2. Interactive taxonomic pie-charts are generated, displaying detailed taxonomic contents of each sample.
   3. Interactive (and searchable) tables are included containing detailed information regarding classified and unclassified scaffolds.
* Virus-typing results are presented with interactive tables.
* An interactive scaffold alignment viewer (IGVjs) is included.

   1. Detailed alignment information.
   2. Depth of coverage graph per scaffold.
   3. GC content graph.
   4. Predicted open reading frames (ORFs).
   5. Identified minority variants.
* All SNP metrics are presented via interactive tables allowing for easy and detailed analysis. 


Audit trail
==================

An audit trail, used for clinical reproducability and logging, is generated for each run and contains the following information:

* A unique fingerprint representing the Jovian code is generated through Github allowing exact reproductions of the analysis. This also allows for retrospective reversions to older versions of the pipeline.
* A unique identifier is generated for each run that can be used for internal archiving and result tracing.
* Database timestamps are logged
* Pipeline parameters (user specified) are logged.

.. image:: https://raw.githubusercontent.com/DennisSchmitz/Jovian/assets/images/rulegraph_Jovian.png



.. _here: https://mybinder.org/v2/gh/DennisSchmitz/Jovian_binder/master?filepath=Notebook_report.ipynb

.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: Installation

   Requirements/requirements
   Installation/installation
   

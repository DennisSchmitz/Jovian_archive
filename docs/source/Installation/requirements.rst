========================
Requirements
========================

Jovian has fairly minimal system requirements.
Jovian is intended for use on powerful (virtual)servers and/or grid-computing configurations.
However, the pipeline will function on most PC's and laptops.

Make sure there is enough storage on your system available as Jovian depends on several databases that have to be downloaded.


Jovian does **not** require root and/or sudo-rights. It is however necessary to have read and write access to the ``/tmp`` folder on your system. This is however unlikely to cause problems as the ``/tmp`` folder is usually free to read from and write to under any condition. We do however advise to consult with your system administrator(s) beforehand.

Supported operating systems
============================
Jovian works exclusively on Linux. Specifically, ``RHEL 7`` (+ ``CentOS 7``) and ``Ubuntu 18.04 LTS``.

Jovian, and all its functions, also works on the Windows Subsystem for Linux with the ``Ubuntu 18.04`` distro. However, stability and functionality are currently not guaranteed on WSL.

System Requirements
=====================

----------
Required hardware
----------

For hardware requirements, storage is currently the only concern as several databases have to be downloaded.
These databases may require up to **500GB** of disk space.
Installation of Jovian requires a total of **6GB** of free disk-space. This however does not include the additionally required disk space that is necessary for actually running the pipeline. 

We advise to have several Terabytes of free storage space available on your system before installation and running Jovian.

----------
Required software
----------

Jovian has some software requirements. However, most of the software listed below is usually pre-installed on Linux systems.

+---------------+---------+
| Software name | Website |
+===============+=========+
| ``git``       | Git_    |
+---------------+---------+
| ``curl``      | Curl_   |
+---------------+---------+
| ``which``     | Which_  |
+---------------+---------+
| ``bzip2``     | Bzip2_  |
+---------------+---------+

Additionally, Jovian requires miniConda_ in order to run. However, during installation of Jovian you will be asked to install miniconda if it doesn't already exist on your system. 

----------
Required databases
----------

Jovian depends on several databases. Jovian is able to download and install these databases for you if these databases do not exist on your system already. However, this has to be initiated manually.

+-------------------------------------+--------------+
| Database name                       | Links        |
+=====================================+==============+
| ``NCBI BLAST NT``                   | Blast_       |
+-------------------------------------+--------------+
| ``NCBI Taxdump``                    | Taxdump_     |
+-------------------------------------+--------------+
| ``NCBI new_taxdump``                | NewTaxdump_  |
+-------------------------------------+--------------+
| ``Virus-Host interaction database`` | VirushostDB_ |
+-------------------------------------+--------------+
| ``Latest human genome`` *           | HuGo_        |
+-------------------------------------+--------------+

\* *We suggest the latest human genome because Jovian is intended for use with clinical samples. You are however able to use any reference you'd like.*



.. Links and URLS:

.. _Git: https://git-scm.com/downloads
.. _Curl: https://curl.haxx.se/
.. _Which: http://savannah.gnu.org/projects/which
.. _Bzip2: http://www.bzip.org/
.. _miniConda: https://docs.conda.io/en/latest/miniconda.html

.. _Blast: ftp://ftp.ncbi.nlm.nih.gov/blast/db/
.. _Taxdump: ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/
.. _NewTaxdump: ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/
.. _VirushostDB: http://www.genome.jp/virushostdb/note.html
.. _HuGo: https://support.illumina.com/sequencing/sequencing_software/igenome.html
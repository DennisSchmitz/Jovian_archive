==================================
Installation
==================================

.. note:: For compatibility information, and information regarding system/software requirements for Jovian.
          Please see the :doc:`Requirements <requirements>`.

Instructions
============

Installation of Jovian is made as simple as possible for end-users, there are two possibilities for installing Jovian.

* Installing Jovian if the required databases are not yet installed.
* Installing Jovian if the required databases are already installed.

Nearly every step during the installation process that requires manual input will only be required once.
Jovian stores the chosen configurations so you don't have to manually enter this information again when you wish to run another analysis.


------------------------------------------------------------------------
Installing Jovian if the required databases are not yet installed
------------------------------------------------------------------------

Jovian is able to download and install the required databases for you in case these are not yet present on your system.
However, Jovian will **not** download the background reference information for you. This has to be done manually for your background reference of choice.

To install Jovian and download/install the required standard databases, follow the steps underneath:

1. Download Jovian from Github with the following command: ``git clone https://github.com/DennisSchmitz/Jovian.git``
2. Change directory to the folder that was just created with ``cd Jovian``.
3. Type (or copy+paste) the following command to start the installation: ``bash jovian -id``

These three steps will launch the installer process of Jovian.


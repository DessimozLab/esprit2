Esprit 2: Detecting split genes
===============================

Esprit 2 is a software package to detect split genes in a proteome using 
a statistical test based on reconstructed gene trees with related reference
genomes. The approach is described in detail in `http://xxx`_.



Installation:
-------------

Clone the git repository from https://github.com/dessimozlab/esprit2

.. code: sh

    git clone https://github.com/dessimozlab/esprit2

Adjust your file load_env according to your environment, i.e. python environment,
paths to required software,... 


Input Data:
-----------

You need to provide a folder with fasta files, one per gene family and an 
orthoxml file. The simplest way to obtain these files is by running 
`OMA Standalone <http://omabrowser.org/standalone>`_ on your dataset. 
The produced files /HierarchicalGroups.orthoxml/ and the folder /HOGFasta/ 
are the two input files you need for Esprit 2.

How to run Esprit 2:
--------------------

As of now, Esprit 2 needs a SunGridEngine (SGE) scheduler. This will likely
change in the future and will be extended to other HPC schedulers.

Run Esprit by 

.. code: sh

    ./pipeline.sh <ID_prefix> <path_to_family_folder> <path_to_orthoxml>

By calling the pipeline.sh script with -h, it will output additional parameters
you can specify together with their default values.


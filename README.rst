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
are the two input files you need for Esprit 2. Please put them to your esprit2 directory. 

How to run Esprit 2:
--------------------

As of now, Esprit 2 needs a SunGridEngine (SGE) scheduler. This will likely
change in the future and will be extended to other HPC schedulers.

Run Esprit by 

.. code: sh

    ./pipeline.sh <ID_prefix> <path_to_family_folder> <path_to_orthoxml>

``ID_prefix`` is a unique species- or chromosome-specific substring of a gene ID. For example, in case of wheat where all gene IDs follow the format ``Traes_chromosomeArm_string`` (e.g., ``Traes_1BL_6EC2AB17D``) or ``TRAES3Bstring`` (e.g., ``TRAES3BF036000300CFD``) for 3B reference assembly, an ``ID_prefix`` could be:
- ``Traes_chromosomeArm``, e.g., ``Traes_1BL`` - for detecting split genes within a chromosome arm, e.g., long arm of chromosome 1B
- ``Traes_chromosome``, e.g., ``Traes_1B`` - for detecting split genes within a chromosome, e.g., chromosome 1B. This will probably yield candidate pairs where one fragment has been assigned to 1BL (long arm) and the other to 1BS (short arm).
- ``Traes`` - for detecting split genes within the whole wheat genome. Please be aware that the set of candidate pairs will contain fragments coming from different chromosomes.
- ``TRAES3B``- for detecting split genes within 3B reference assembly

By calling the pipeline.sh script with -h, it will output additional parameters
you can specify together with their default values.


Output files:
-------------

.. code: sh

	collapsing_results.txt

columns: gene1, gene2, sister taxa before collapsing (True/False), sister taxa after collapsing (True/False)

.. code: sh

	lrt_summaries.tar.gz, lrt_summary.txt, missing_lk.txt

tar.gz contains a summary per case tested, lrt_summary.txt provides test statistics and p-values for all cases, missing_lk.txt indicates cases where the tree likelihood wasn't computed (please have a look at these computations and investigate what went wrong) 	

.. code: sh
	
	predictions_ambiguous.txt, predictions_unambiguous.txt

contain gene IDs for predictions

.. code: sh
	
	alignment_positions.txt 

columns: HOG ID, gene1, gene2, start position of gene1 in the MSA, end position of gene1 in the MSA, start position of gene2 in the MSA, end position of gene2 in the MSA, overlap start position (or -1 if no overlap), overlap end position (or -1 if no overlap), %overlap of aligned gene1, %overlap of aligned gene2

.. code: sh
	
	cuts.txt

columns: HOG ID, gene1, gene2, their cut/middle position in the alignment

.. code: sh

	mapping.txt

mapping between OMA IDs and IWGSC IDs

.. code: sh

	sequence_lengths.txt

columns: HOG ID, gene1, length of gene1, gene2, length of gene2 
contains also pairs with short sequence(s) which didn't pass the min sequence length criteria

.. code: sh

	aln_c.tar.gz, aln.tar.gz, phy_c.tar.gz, phy.tar.gz

contain aligned families in FASTA format (aln_c, aln) and phylip (phy_c, phy). aln_c and phy_c contain families with n-1 sequences whereas aln and phy contain n sequences

.. code: sh

	hog_aln.tar.gz

alignments of HOGs which contain at least 2 wheat genes from the chromosome of interest

.. code: sh

	bootstrap_aln.tar.gz, bootstrap_s_aln.tar.gz, bootstrap_phy.tar.gz, bootstrap_s_phy.tar.gz

similar as above but for bootstrap samples. bootstrap_aln.tar.gz and bootstrap_phy.tar.gz contain samples with n-1 sequences whereas bootstrap_s_aln.tar.gz and bootstrap_s_phy.tar.gz contain samples with n sequences

.. code: sh
	
	collapsed.tar.gz

contains trees after collapsing

.. code: sh

	n_1_res.tar.gz, n_notop_res.tar.gz, n_top_res.tar.gz, n_1_b_res.tar.gz, n_b_notop_res.tar.gz, n_b_top_res.tar.gz

contain stats output from FastTree

.. code: sh
	
	n_1_trees.tar.gz, n_trees_notop.tar.gz, n_1_b_trees.tar.gz

contain FastTree trees

.. code: sh
	
	n_1_trees_s.tar.gz, n_1_b_trees_s.tar.gz

contain input topologies for tree reconstructions with input topology
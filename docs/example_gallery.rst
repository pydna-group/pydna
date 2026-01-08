Example gallery
===============

.. toctree::
   :maxdepth: 1
   :hidden:

   notebooks/Example_Restriction.ipynb
   notebooks/Example_Gibson.ipynb
   notebooks/Example_CRISPR.ipynb

Below are some examples that show the functionality of pydna in real-world scenarios.

* `Example_Restriction <./notebooks/Example_Restriction.html>`_: PCRing a gene out of the genome, and cloning into a vector using restriction and ligation.
* `Example_Gibson <./notebooks/Example_Gibson.html>`_: Gibson assembly of *R. cellulolyticum* genomic fragments into a plasmid, from the original Gibson assembly paper `doi: 10.1038/nmeth.1318 <https://www.nature.com/articles/nmeth.1318>`_.
* `Example_CRISPR <./notebooks/Example_CRISPR.html>`_: Using CRISPR with homologous recombination to delete genes by making two cuts in the genome, and repair it with an oligo. Used in the industrially relevant *K. phaffi*.

External examples
-----------------

The below examples use the extra dependency `teemi <https://github.com/hiyama341/teemi>`_:

* `Example_HT_cazyme_primer_design <https://github.com/hiyama341/cazyme_primer_design/blob/main/notebooks/00-cazyme_primer_design.ipynb>`_: Design primers for a high-throughput CAZyme library.
* `Example_designing_primers_for_kozak_library <https://github.com/hiyama341/teemi/blob/main/colab_notebooks/09_2_BUILD_Combinatorial_library.ipynb>`_: We explore the combinatorial space of the most abundant kozak sequences and make repair-primers for the experiments.


Bioengineering/SynBio projects where you can apply pydna
--------------------------------------------------------

* `02 Primer directed mutagenesis with pydna <https://github.com/hiyama341/teemi/blob/main/teemi_cad_workflows/notebooks/02_primer_directed_mutagenesis.ipynb>`_
* `03 Investigate plastic-degrading enzymes with pydna <https://github.com/hiyama341/teemi/blob/main/teemi_cad_workflows/notebooks/03_investigate_plastic_degrading_enzymes.ipynb>`_
* `04 Promoter library design in *Saccharomyces cerevisiae* <https://github.com/hiyama341/teemi/blob/main/teemi_cad_workflows/notebooks/04_promoter_library_design_in_s_cerevisiae.ipynb>`_

# Example gallery

Below are some examples that show the functionality of pydna in real-world scenarios.

* [Example_Restriction](./markdown_notebooks/Example_Restriction.md): PCRing a gene out of the genome, and cloning into a vector using restriction and ligation.
* [Example_Gibson](./markdown_notebooks/Example_Gibson.md): Gibson assembly of _R. cellulolyticum_ genomic fragments into a plasmid, from the original Gibson assembly paper [doi: 10.1038/nmeth.1318](https://www.nature.com/articles/nmeth.1318).
* [Example_CRISPR](./markdown_notebooks/Example_CRISPR.md): Using CRISPR with homologous recombination to delete genes by making two cuts in the genome, and repair it with an oligo. Used in the industrially relevant _K. phaffi_.

## External examples

The below examples use the extra dependency [teemi](https://github.com/hiyama341/teemi):

* [Example_HT_cazyme_primer_design](https://github.com/hiyama341/cazyme_primer_design/blob/main/notebooks/00-cazyme_primer_design.ipynb): Design primers for a high-throughput CAZyme library.
* [Example_designing_primers_for_kozak_library](https://github.com/hiyama341/teemi/blob/main/colab_notebooks/09_2_BUILD_Combinatorial_library.ipynb): We explore the combinatorial space of the most abundant kozak sequences and make repair-primers for the experiments.

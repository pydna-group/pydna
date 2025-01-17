## README Example

This notebook contains the example shown in the README file.

<a target="_blank" href="https://colab.research.google.com/github/pydna-group/pydna/blob/master/docs/notebooks/Example_Restriction.ipynb">
  <img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab"/>
</a>


```python
%%capture
# Install pydna (only when running on Colab)
import sys
if 'google.colab' in sys.modules:
    %pip install pydna

```


```python
from pydna.dseqrecord import Dseqrecord
# Let's create a DNA sequence record, and add a feature to it
dsr = Dseqrecord("ATGCAAACAGTAATGATGGATGACATTCAAAGCACTGATTCTATTGCTGAAAAAGATAAT")
dsr.add_feature(x=0, y=60,type="gene", label="my_gene") # We add a feature to highlight the sequence as a gene
dsr.figure()

```




    Dseqrecord(-60)
    [48;5;11mATGCAAACAGTAATGATGGATGACATTCAAAGCACTGATTCTATTGCTGAAAAAGATAAT[0m
    TACGTTTGTCATTACTACCTACTGTAAGTTTCGTGACTAAGATAACGACTTTTTCTATTA




```python
# This is how it would look as a genbank file
print(dsr.format("genbank"))
```

    LOCUS       name                      60 bp    DNA     linear   UNK 01-JAN-1980
    DEFINITION  description.
    ACCESSION   id
    VERSION     id
    KEYWORDS    .
    SOURCE      .
      ORGANISM  .
                .
    FEATURES             Location/Qualifiers
         misc            1..60
                         /type="gene"
                         /label="my_gene"
    ORIGIN
            1 atgcaaacag taatgatgga tgacattcaa agcactgatt ctattgctga aaaagataat
    //



```python
# Now let's design primers to amplify it
from pydna.design import primer_design
# limit is the minimum length of the primer, target_tm is the desired melting temperature of the primer
amplicon = primer_design(dsr, limit=13, target_tm=55)
# Let's print the primers, and a figure that shows where they align with the template sequence
print("forward primer:", amplicon.forward_primer.seq)
print("reverse primer:", amplicon.reverse_primer.seq)
amplicon.figure()
```

    forward primer: ATGCAAACAGTAATGATGGA
    reverse primer: ATTATCTTTTTCAGCAATAGAATCA





    5ATGCAAACAGTAATGATGGA...TGATTCTATTGCTGAAAAAGATAAT3
                            |||||||||||||||||||||||||
                           3ACTAAGATAACGACTTTTTCTATTA5
    5ATGCAAACAGTAATGATGGA3
     ||||||||||||||||||||
    3TACGTTTGTCATTACTACCT...ACTAAGATAACGACTTTTTCTATTA5




```python
# Let's say we don't want to just amplify it, but we want to add restriction sites to it!

from pydna.amplify import pcr
# We add the restriction sites to the primers
forward_primer = "ccccGGATCC" + amplicon.forward_primer
reverse_primer = "ttttGGATCC" + amplicon.reverse_primer

# We do the PCR
pcr_product = pcr(forward_primer, reverse_primer, dsr)
# The PCR product is of class `Amplicon`, a subclass of `Dseqrecord`.
# When doing a figure, it shows where primers anneal.
pcr_product.figure()
```




              5ATGCAAACAGTAATGATGGA...TGATTCTATTGCTGAAAAAGATAAT3
                                      |||||||||||||||||||||||||
                                     3ACTAAGATAACGACTTTTTCTATTACCTAGGtttt5
    5ccccGGATCCATGCAAACAGTAATGATGGA3
               ||||||||||||||||||||
              3TACGTTTGTCATTACTACCT...ACTAAGATAACGACTTTTTCTATTA5




```python
# If we want to see the sequence more clearly, we can turn it into a `Dseqrecord`
pcr_product = Dseqrecord(pcr_product)
pcr_product.figure()
```




    Dseqrecord(-80)
    ccccGGATCC[48;5;11mATGCAAACAGTAATGATGGATGACATTCAAAGCACTGATTCTATTGCTGAAAAAGATAAT[0mGGATCCaaaa
    ggggCCTAGGTACGTTTGTCATTACTACCTACTGTAAGTTTCGTGACTAAGATAACGACTTTTTCTATTACCTAGGtttt




```python
from Bio.Restriction import BamHI # cuts GGATCC
# a, payload, c are the cut fragments
a, payload, c = pcr_product.cut (BamHI)
print(a.figure())
print()
print (payload.figure())
print()
print(c.figure())


```

    Dseqrecord(-9)
    [48;5;11m[0mccccG    
    ggggCCTAG
    
    Dseqrecord(-70)
    GATCC[48;5;11mATGCAAACAGTAATGATGGATGACATTCAAAGCACTGATTCTATTGCTGAAAAAGATAAT[0mG    
        GTACGTTTGTCATTACTACCTACTGTAAGTTTCGTGACTAAGATAACGACTTTTTCTATTACCTAG
    
    Dseqrecord(-9)
    [48;5;11m[0mGATCCaaaa
        Gtttt



```python
# We create a circular vector to insert the amplicon into
vector = Dseqrecord("aatgtttttccctCCCGGGcaaaatAGATCTtgctatgcatcatcgatct", circular=True, name="vect")
vector.figure()
```




    Dseqrecord(o50)
    [48;5;11m[0maatgtttttccctCCCGGGcaaaatAGATCTtgctatgcatcatcgatct
    ttacaaaaagggaGGGCCCgttttaTCTAGAacgatacgtagtagctaga




```python
from Bio.Restriction import BglII # cuts AGATCT
linear_vector_bgl = vector.cut(BglII)[0] # Linearize the vector at BglII (produces only one fragment)

# Ligate the fragment of interest to the vector, and call looped() to circularize it
# synced is used to place the origin coordinate (0) in the same place for rec_vector and vector
rec_vector= (linear_vector_bgl + payload).looped().synced(vector)
rec_vector.figure()

```




    Dseqrecord(o116)
    aatgtttttccctCCCGGGcaaaatAGATCC[48;5;11mATGCAAACAGTAATGATGGATGACATTCAAAGCACTGATTCTATTGCTGAAAAAGATAAT[0mGGATCTtgctatgcatcatcgatct
    ttacaaaaagggaGGGCCCgttttaTCTAGGTACGTTTGTCATTACTACCTACTGTAAGTTTCGTGACTAAGATAACGACTTTTTCTATTACCTAGAacgatacgtagtagctaga




```python
# Let's simulate a Gibson assembly
from pydna.assembly import Assembly

fragments = [
    Dseqrecord('aatgtttttccctCACTACGtgctatgcatcat', name="fragment_A"),
    Dseqrecord('tgctatgcatcatCTATGGAcactctaataatg', name="fragment_B"),
    Dseqrecord('cactctaataatgTTACATAaatgtttttccct', name="fragment_C"),
]

# limit is the min. homology length between fragments in the assembly
asm = Assembly(fragments, limit=10)

# From the assembly object, which can generate all possible products, get a circular
product, *rest = asm.assemble_circular()

# We can print a figure that shows the overlaps between fragments
product.figure()

```




     -|fragment_A|13
    |             \/
    |             /\
    |             13|fragment_B|13
    |                           \/
    |                           /\
    |                           13|fragment_C|13
    |                                         \/
    |                                         /\
    |                                         13-
    |                                            |
     --------------------------------------------




```python
# Or show the final sequence:
Dseqrecord(product).figure()
```




    Dseqrecord(o60)
    [48;5;11m[0maatgtttttccctCACTACGtgctatgcatcatCTATGGAcactctaataatgTTACATA
    ttacaaaaagggaGTGATGCacgatacgtagtaGATACCTgtgagattattacAATGTAT



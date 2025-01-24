genocad.com seems still down

But I managed to download the supplementary files from the [paper](https://www.ncbi.nlm.nih.gov/pubmed/25925571).

I put them [here](https://github.com/BjornFJohansson/pydna/tree/py3dev/scripts/files_for_annotation_function)

A number of sbol files seems to hold most or all of the data:

labhost_All.xml
labhost_Aspergillus_nidulans.xml
labhost_Bacillus_subtilis.xml
labhost_Drosophila_melanogaster.xml
labhost_Escherichia_Coli.xml
labhost_Gram-negative_bacteria.xml
labhost_Insect_Cells.xml
labhost_Kluyveromyces_lactis.xml
labhost_Mammalian_Cells.xml
labhost_Pichia_pastoris.xml
labhost_Plant_Cells.xml
labhost_Saccharomyces_cerevisiae.xml
labhost_Schizosaccharomyces_pombe.xml
labhost_Unspecified.xml

Looking closely, the labhost_All.xml should hold all of the features.

These seem to be the same as the ones I got by mail from Tim.

They contain this kind of element:

<component>
   <DnaComponent ns2:about="http://genocad.org/SBOL/doi:10.6084/m9.figshare.1159038/Part:AmpR_promoter-009">
     <ns2:type ns2:resource="so:0000167"/>
     <displayId>AmpR_promoter-009</displayId>
     <name>AmpR promoter-009</name>
     <description/>
     <dnaSequence>
        <DnaSequence ns2:about="http://genocad.org/SBOL/doi:10.6084/m9.figshare.1159038/Part:AmpR_promoter-009#seq">
           <nucleotides>cgcggaacccctatttgtttatttttctaaatacattcaaatatgtatccgctcatgagacaataaccctgataaatgcttcaataatattgaaaaaggaagagt</nucleotides>
        </DnaSequence>
     </dnaSequence>
    </DnaComponent>
</component>

Which seems to contain all that is necessary for describing a feature.

These files could provide a way to store pydna features, but what is the best way to use them.

A simple way would be to use a [sbol parser](https://pypi.python.org/pypi?%3Aaction=search&term=sbol&submit=search)
and parse into Biopython SeqFeature objects and strings and then use some algorithm to match all features against
a sequence, probably with a new method for the Dseqrecord class.

It might be overkill to pull in pySBOL as a dependecy for this as this is a wrapper around libSBOL which is C++.

There is snekbol which is pure python https://github.com/tjomasc/snekbol

`
(bjorn3) bjorn@bjorn-ThinkPad-T450s:~$ pip install snekbol
Collecting snekbol
  Downloading snekbol-0.1.2.tar.gz
Collecting lxml==3.7.3 (from snekbol)
  Downloading lxml-3.7.3-cp35-cp35m-manylinux1_x86_64.whl (7.1MB)
    100% |████████████████████████████████| 7.1MB 220kB/s
Collecting rdflib==4.2.2 (from snekbol)
  Downloading rdflib-4.2.2-py3-none-any.whl (344kB)
    100% |████████████████████████████████| 348kB 3.4MB/s
Collecting validators==0.11.2 (from snekbol)
  Downloading validators-0.11.2.tar.gz
Collecting isodate (from rdflib==4.2.2->snekbol)
  Downloading isodate-0.5.4.tar.gz
Requirement already satisfied: pyparsing in ./anaconda3/envs/bjorn3/lib/python3.5/site-packages (from rdflib==4.2.2->snekbol)
Requirement already satisfied: six>=1.4.0 in ./anaconda3/envs/bjorn3/lib/python3.5/site-packages (from validators==0.11.2->snekbol)
Requirement already satisfied: decorator>=3.4.0 in ./anaconda3/envs/bjorn3/lib/python3.5/site-packages (from validators==0.11.2->snekbol)
Installing collected packages: lxml, isodate, rdflib, validators, snekbol
  Found existing installation: lxml 3.7.2
    DEPRECATION: Uninstalling a distutils installed project (lxml) has been deprecated and will be removed in a future version. This is due to the fact that uninstalling a distutils project will only partially uninstall the project.
    Uninstalling lxml-3.7.2:
      Successfully uninstalled lxml-3.7.2
  Running setup.py install for isodate ... done
  Running setup.py install for validators ... done
  Running setup.py install for snekbol ... done
Successfully installed isodate-0.5.4 lxml-3.7.3 rdflib-4.2.2 snekbol-0.1.2 validators-0.11.2
(bjorn3) bjorn@bjorn-ThinkPad-T450s:~$
`

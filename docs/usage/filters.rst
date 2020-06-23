Filters panel
=============

Understand the variants annotations
+++++++++++++++++++++++++++++++++++

## Small variants
Annotated for gene consequence using snpEff and then annotated with our collection of regulatory regions

## Structural variants
Annotated for gene consequence as
- intronic_sv: overlap only introns
- exonis_sv: overlap at least 1bp of a gene exon

The final list of candidate variant records is then generated as follows:
- Most severe consequence is selected for gene affecting variants, so each variant is reported with 1 gene consequence
- Each variant can have multiple regulatory annotations if it overlap multiple regulatory regions

So for example the same variant can appear as multiple records in the candidate vars table:
``chr1    10000     G    T    5UTR`` 
``chr1    10000     G    T    enhancer``

The filtering process is then applied on these single variants records to be able to finely tune the desired output

Variants filters
++++++++++++++++

Variants filters are applied separately for specific variant group, so that each group of filters act only on the relevant variants.
Variants group are defined as follows:

1. splicing
     Variants in splicing site or splicing region (10 bp from exon-intron junctions)
2. missense
3. regulatory regions 
     Any variant located in one of the annotated regulatory regions (enhancer, promoter, silencer, insulator)
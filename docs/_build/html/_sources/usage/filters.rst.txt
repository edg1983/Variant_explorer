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




            pass_vars <- as.data.frame(vars_in_regions %>% 
                filter(! (
                    d_score < input$d_score_filter | 
                    consequence %nin% RV$accepted_consequence |
                    max_pop_af > input$MaxPopAF_filter |
                    cohort_af > input$CohortAF_filter |
                    SegDup %nin% var_reg_anno$SegDup |
                    LowComplexity %nin% var_reg_anno$LowComplexity |
                    TopVariableGenes %nin% var_reg_anno$TopVariableGenes |
                    (reg_type == "splicing" & 
                         (SpliceAI_SNP_SpliceAI_max < input$spliceAI_filter & 
                              SpliceAI_INDEL_SpliceAI_max < input$spliceAI_filter) ) |
                    (consequence == "missense_variant" &
                         (CADD_PhredScore < input$CADD_filter & 
                              REVEL_score < input$REVEL_filter &
                              MCAP_score < input$MCAP_filter &
                              DANN_score < input$DANN_filter)) |
                    (consequence %in% reg_vars &
                         (ReMM_score < input$ReMM_filter |
                              LinSight < input$LinSight_filter |
                              PhyloP100 < input$PhyloP100_filter |
                              LoF_tolerance > input$LoFtolerance_filter) &
                         rec_id %nin% accepted_reg_recid &
                         reg_type %nin% RV$accepted_connected_gene &
                         TFBS %nin% NC_reg_anno$TFBS &
                         DNase %nin% NC_reg_anno$DNase &
                         UCNE %nin% NC_reg_anno$UCNE)
                )))
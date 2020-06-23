Configuration
=============

Variant Explorer is fully configurable editing the 2 configuration files present in the app folder.

1. **App_configuration.json**: configure app options for plotting and some expected columns
2. **Filters_settings.json**: configure filters and presets

However, some columns are always expected in the input data:

- rec_id: unique ID of the specific line
- var_id: unique ID of the variant
- consequence: consequence of the variant
- var_type: variant type (like SNV, INDEL, etc)
- chrom,pos,ref,alt: genomic coordinates, ref and alt alleles
- gene: affected gene
- GT_sampleID: genotype for each sample
- GQ_sampleID: GQ for each sample

App_configuration
+++++++++++++++++

The app configuration is a standard json file organized as follows

.. code-block:: JSON

    {
    "var_groups": {
        "consequence_groups": {
        "group1": ["consequence1","consequence2", "..."]
        },
        "var_type_groups": {
        "group1": ["var_type1","var_type2", "..."]
        }
    },
    "segregation_cols": {
        "het_affected": "het_aff", 
        "het_unaffected": "het_unaff", 
        "hom_affected": "hom_aff", 
        "hom_unaffected": "hom_unaff", 
        "comphet_affected": "comphet_aff",
        "dnm_affected": "sup_dnm"
    },
    "plot_axes": {
        "genes_axes_options": {
        "label1": "col_name1", 
        "label2": "col_name2"
        },
        "variants_axes_options": {
        "label1": "col_name1",
        "label2": "col_name2"
        },
        "variants_bar_options": {
        "label1": "col_name1",
        "label2": "col_name2"
        }
    },
    "fill_na": {
        "fill_na_vars" : {
        "col_name1": "replace value1",
        "col_name2": "replace value2"
        },
        "fill_na_genes" : {
        "col_name1": "replace value1",
        "col_name2": "replace value2"
        }
    }
    }

Configuration sections:
-----------------------

  "var_groups":
    "consequence_groups":  
        a series of lists defining variants groups based on variant consequence
        ``"lof_vars": ["frameshift", "exonic_sv","stop lost","stop gain","start loss","splice_acceptor_variant","splice_donor_variant"]``

    "var_type_groups":
        a series of lists defining variants groups based on variant type
        ``"sv_vars": ["DEL","DUP","INV","DEL:ME"]``

  "segregation_cols":
    definition of the columns containing counts of affected/unaffected carriers for each genotype
    
    .. code-block:: JSON

        {
        "het_affected": "het_aff", 
        "het_unaffected": "het_unaff", 
        "hom_affected": "hom_aff", 
        "hom_unaffected": "hom_unaff", 
        "comphet_affected": "comphet_aff",
        "dnm_affected": "sup_dnm"
        }

  "plot_axes":
    This section define the columns used to generate overview and filters plots for genes and variants
    
    "genes_axes_options":
        Columns used in the genes filters scatter-plots, in the formate label: column_name.
        The labels are shown in a drop-down box to configure the corresponding plot
        Example:

        .. code-block:: JSON
            
            {
            "GADO Zscore": "gado_zscore", 
            "Exomiser Pheno score": "exomiser_gene_pheno_score",
            "gnomAD pLI": "pLI_gnomad",
            "GDI phred": "GDI_phred",
            "RVIS intolerance": "RVIS",
            "EDS reg space score": "EDS"
            }

    "variants_axes_options":
        Columns used in the variants filters scatter-plots, in the formate label: column_name.
        The labels are shown in a drop-down box to configure the corresponding plot
        Example:
      
        .. code-block:: JSON

            {
            "Maximum population AF": "max_pop_af",
            "CADD phred": "CADD_PhredScore", 
            "DANN score": "DANN_score",
            "ReMM score": "ReMM_score",
            "SpliceAI score": "SpliceAI_SNP_SpliceAI_max",
            "REVEL score": "REVEL_score",
            "MCAP score": "MCAP_score",
            "LinSight score": "LinSight"
            }

    "variants_bar_options":
        Columns used in the variants filters bar plots, in the formate label: column_name.
        The labels are shown in a drop-down box to configure the corresponding plot
        Example:

        .. code-block:: JSON

            { 
            "Region type": "reg_type",
            "Variant consequence": "consequence",
            "Variant type": "var_type",
            "PanelApp": "PanelApp",
            "Chromosome": "chr"
            }
  "fill_na":
    This section define how to fill NA values for gene and variants scores.
    For each column it is possible to assign a specific value to replace NAs
    
    "fill_na_vars":
        NA replacement values for variants table, in the format col_name: replace_value
        Example

        .. code-block:: JSON

            {
            "max_pop_af": 0,
            "cohort_af": 0,
            "DANN_score": 99,
            "CADD_PhredScore": 99,
            "PhyloP100": -99
            }

    "fill_na_genes":
        NA replacement values for genes table, in the format col_name: replace_value
        Example

        .. code-block:: JSON

            {
            "GDI_phred": -99,
            "RVIS": -99,
            "pLI_gnomad": 99,
            "EDS": 99
            }

Filters configuration
+++++++++++++++++++++

The app configuration is a standard json file organized as follows

.. code-block:: JSON
    
    "GENES": {
        "DEFINITIONS": {
            "numeric_fields": { 
                "col_name1": ["label1", ">=", "default_value"],
                "col_name2": ["label2", "<", "default_value"]
            }
        },
        "GROUPS": {
            "global": {
                "associated_values": ["col_name1", "col_name2"]    
            }
        },
        "PRESETS": {
            "global": {
                "default": "preset1",  
                "configuration": {
                    "preset1": {
                        "logic": "OR",
                        "values" : {
                            "column1": value1,
                            "column2": value1
                        }
                    },
                    "preset2": {
                        "logic": "OR",
                        "values" : {
                            "column1": value2,
                            "column2": value2
                        }
                    }
                }
            }
        }
    },
    "VARIANTS": {
        "DEFINITIONS": {
            "numeric_fields": { 
                "col_name1": ["label1", ">=", "default_value"],
                "col_name2": ["label2", "<", "default_value"]
            },
            "factors_fields": {
                "col_name1": ["label1", "%in%", "default_value"],
                "col_name2": ["label2", "grep", "default_value"]
            },
            "binary_fields": {
                "col_name1": ["label1", "include", "default_true/false"],
                "col_name2": ["label2", "exclude", "default_true/false"]
            }
        },
        "FACTORSMAP": {
            "col_name1": {
                "label1": "factor1", 
                "label2": "factor2", 
            },
            "reg_type": {
                "Closest gene":  "closest_gene", 
                "From database": "reg_db"
            }
        },
        "GROUPS": {
            "global": {
                "associated_values": ["col_name1","col_name2"]
            },
            "group1": {
                "definition": ["consequence", ["value1","value2"]],
                "associated_values": ["col_name1","col_name2"]
            }
        },
        "PRESETS": {
        "global": {
            "default": "strict",
            "configuration": {
                "strict": {
                    "logic": "AND",
                    "values": {
                        "column1": value1,
                        "column2": value1
                    }
                },
                "balanced": {
                    "logic": "AND",
                    "values": {
                        "column1": value2,
                        "column2": value2
                    }
                }
            }
        },
        "group1": {
            "default": "balanced",  
            "configuration": {
                "strict": {
                    "logic": "AND",
                    "values": {
                        "column1": value1,
                        "column2": value1
                    }
                },
                "balanced": {
                    "logic": "AND",
                    "values": {
                        "column1": value2,
                        "column2": value2
                    }
                }
            }
        }
    },
    "TOOLTIPS": {
        "col_name1": "tooltip_text_1",
        "col_name2": "tooltip_text_2"
    }
    
    
# Reference

API and command-line reference for Varlens.

## Unified CLI

::: varlens.commands.main
    options:
      members:
        - make_parser
        - run

::: varlens.commands.help
    options:
      members:
        - HelpArgumentParser
        - HelpFormatter
        - color_enabled
        - style_help
        - style_table_header
        - classify_status
        - style_status
        - enable_argcomplete

## Genomic loci

::: varlens.locus.Locus

::: varlens.loci_util.Loci

## Variant and read utilities

::: varlens.variants_util
    options:
      members:
        - union_with_variant_loci
        - variant_loci

::: varlens.util
    options:
      members:
        - name_sources
        - strip_extensions

::: varlens.errors
    options:
      members:
        - InputError
        - input_errors

## Read support and MHC binding

::: varlens.support
    options:
      members:
        - variant_support

::: varlens.mhc_binding
    options:
      members:
        - binding_affinities
        - fragment_predictions

::: varlens.fusions
    options:
      members:
        - analyze_fusion

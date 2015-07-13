
import argparse

def load_variant_collections_parser():
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--variants", nargs="+", required=True)
    parser.add_argument("--ensembl-version")
    parser.add_argument("--variant-filter")
    return parser

def load_read_sets_parser():
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--reads", nargs="+", default=[])
    return parser


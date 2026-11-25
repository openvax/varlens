import argparse

import pytest
import varcode

from varlens import util, variants_util
from varlens.commands import allele_support, fusions, main, reads, variants
from varlens.commands.arguments import add_verbose_arg, add_version_arg
from varlens.errors import InputError
from varlens.loci_util import Loci
from varlens.locus import Locus
from varlens.version import __version__

from . import data_path


def test_version_and_verbose_args(capsys):
    parser = argparse.ArgumentParser(prog="example")
    add_version_arg(parser)
    add_verbose_arg(parser)
    assert parser.parse_args([]).verbose is False
    assert parser.parse_args(["-v"]).verbose is True
    assert parser.parse_args(["--verbose"]).verbose is True
    with pytest.raises(SystemExit) as exit_info:
        parser.parse_args(["--version"])
    assert exit_info.value.code == 0
    assert capsys.readouterr().out == "example %s\n" % __version__


@pytest.mark.parametrize("command", [main, allele_support, fusions, reads, variants])
def test_every_command_prints_its_version(command, capsys):
    with pytest.raises(SystemExit) as exit_info:
        command.run(["--version"])
    assert exit_info.value.code == 0
    assert capsys.readouterr().out.rstrip("\n").endswith(" " + __version__)


@pytest.mark.parametrize("subcommand", ["variants", "reads", "allele-support", "fusions"])
def test_unified_cli_subcommand_version(subcommand, capsys):
    with pytest.raises(SystemExit) as exit_info:
        main.run([subcommand, "--version"])
    assert exit_info.value.code == 0
    assert capsys.readouterr().out.strip() == f"varlens {subcommand} {__version__}"


def test_unified_cli_requires_subcommand(capsys):
    with pytest.raises(SystemExit) as exit_info:
        main.run([])
    assert exit_info.value.code == 2
    assert "COMMAND" in capsys.readouterr().err


def test_unified_cli_dispatches_variants_reads_and_allele_support(capsys):
    vcf = data_path("CELSR1/vcfs/vcf_1.vcf")
    bam = data_path("CELSR1/bams/bam_1.bam")

    variants.run([vcf])
    legacy_variants = capsys.readouterr().out
    main.run(["variants", vcf])
    assert capsys.readouterr().out == legacy_variants

    reads.run([bam, "--header"])
    legacy_reads = capsys.readouterr().out
    main.run(["reads", bam, "--header"])
    assert capsys.readouterr().out == legacy_reads

    allele_support.run([bam, "--locus", "22:46931061"])
    legacy_support = capsys.readouterr().out
    main.run(["allele-support", bam, "--locus", "22:46931061"])
    assert capsys.readouterr().out == legacy_support


def test_unified_cli_dispatches_fusions(capsys):
    pytest.importorskip("isovar.fusion")
    fusion_input = data_path("fusions/ATP5MG--KMT2A.input.json.gz")
    fusions.run([fusion_input])
    legacy_report = capsys.readouterr().out
    main.run(["fusions", fusion_input])
    assert capsys.readouterr().out == legacy_report


def test_unified_cli_handles_input_error_cleanly(capsys):
    with pytest.raises(SystemExit) as exit_info:
        main.run(["fusions", "not-a-fusion.txt"])
    assert exit_info.value.code == 2
    err = capsys.readouterr().err
    assert "varlens fusions: error:" in err
    assert ".json or .json.gz" in err


@pytest.mark.parametrize("command, arguments", [
    (allele_support, []),
    (reads, [data_path("CELSR1/bams/bam_0.bam")]),
    (variants, []),
])
def test_verbose_flag(command, arguments):
    assert command.parser.parse_args(arguments).verbose is False
    assert command.parser.parse_args(arguments + ["-v"]).verbose is True


def test_union_with_variant_loci(tmp_path):
    snv = varcode.Variant("22", 46931060, "A", "C", genome="GRCh37")
    insertion = varcode.Variant("20", 10008921, "C", "CA", genome="GRCh37")
    assert list(variants_util.union_with_variant_loci(None, [snv], "test")) == [
        Locus("22", 46931059, 46931060)]

    explicit = Loci([Locus("22", 99, 100)])
    combined = variants_util.union_with_variant_loci(
        explicit, iter([snv, insertion]), "test")
    assert list(combined) == [
        Locus("22", 99, 100),
        Locus("22", 46931059, 46931060),
        Locus("20", 10008921, 10008921),
    ]
    assert list(explicit) == [Locus("22", 99, 100)]

    vcf = tmp_path / "deletion.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n##reference=GRCh38\n"
        '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type">\n'
        '##INFO=<ID=END,Number=1,Type=Integer,Description="End">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr7\t200\tdel\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=300\n")
    structural = variants_util.load_as_dataframe(str(vcf)).variant.tolist()
    with pytest.raises(InputError, match="test operation does not support"):
        variants_util.union_with_variant_loci(None, structural, "test operation")


def test_name_sources():
    paths = ["data/tumor.bam", "data/normal.bam"]
    args = ("read_source_name", "read source")
    assert util.name_sources(paths, None, *args) == ["tumor", "normal"]
    assert util.name_sources(paths, [], *args) == ["tumor", "normal"]
    assert util.name_sources(paths, ["t", "n"], *args) == ["t", "n"]
    assert util.name_sources(paths, ["both"], *args) == ["both", "both"]
    assert util.name_sources([], None, *args) == []
    with pytest.raises(ValueError, match="read_source_name"):
        util.name_sources(paths + ["data/other.bam"], ["t", "n"], *args)

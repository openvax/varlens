"""Subprocess argcomplete, lazy-import hygiene, and shared cyan help styling (#54)."""

from __future__ import annotations

import argparse
import io
import os
from pathlib import Path
import re
import subprocess
import sys
import tomllib

import pytest

from varlens.commands import allele_support, fusions, reads, variants
from varlens.commands.arguments import (
    READ_ATTRIBUTE_NAMES,
    iter_read_filter_specs,
    parse_int_list,
    parse_mhc_predictor_names,
)
from varlens.commands.help import (
    HelpArgumentParser,
    color_enabled,
    style_help,
    style_status,
    style_table_header,
)

REPO_ROOT = Path(__file__).resolve().parents[1]

HEAVY_MODULES = (
    "pandas",
    "numpy",
    "pysam",
    "pyensembl",
    "varcode",
    "mhctools",
    "topiary",
    "pyarrow",
    "torch",
    "pyfaidx",
    "isovar",
)


class FakeTTY(io.StringIO):
    """In-memory text stream reporting ``isatty() == True``."""

    def isatty(self) -> bool:
        return True


def build_subprocess_env(**extra: str) -> dict[str, str]:
    """Build an environment with the worktree first on ``PYTHONPATH``."""
    env = os.environ.copy()
    existing = env.get("PYTHONPATH", "")
    env["PYTHONPATH"] = (
        f"{REPO_ROOT}{os.pathsep}{existing}" if existing else str(REPO_ROOT)
    )
    env.update(extra)
    return env


def run_completion_subprocess(
    tmp_path: Path, entry_module: str, comp_line: str
) -> set[str]:
    """Run ``entry_module.run()`` under ``_ARGCOMPLETE=1`` and assert zero heavy imports."""
    out_file = tmp_path / "argcomplete_out.txt"
    env = build_subprocess_env(
        _ARGCOMPLETE="1",
        _ARGCOMPLETE_IFS="\n",
        _ARGCOMPLETE_STDOUT_FILENAME=str(out_file),
        COMP_LINE=comp_line,
        COMP_POINT=str(len(comp_line)),
    )
    script = (
        "import importlib, sys\n"
        "import argcomplete\n"
        f"mod = importlib.import_module({entry_module!r})\n"
        "original = argcomplete.autocomplete\n"
        "def checked_autocomplete(parser, **kwargs):\n"
        "    def exit_method(code):\n"
        f"        loaded = [m for m in {HEAVY_MODULES!r} if m in sys.modules]\n"
        "        assert not loaded, f'Heavy modules loaded during completion: {loaded}'\n"
        "        raise SystemExit(code)\n"
        "    return original(parser, exit_method=exit_method, **kwargs)\n"
        "argcomplete.autocomplete = checked_autocomplete\n"
        "mod.run()\n"
    )
    result = subprocess.run(
        [sys.executable, "-c", script],
        cwd=str(REPO_ROOT),
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert result.stdout == ""
    assert result.stderr == ""
    return {
        line.strip()
        for line in out_file.read_text(encoding="utf-8").splitlines()
        if line.strip()
    }


@pytest.mark.parametrize(
    "entry_module,comp_line,expected_flags",
    [
        (
            "varlens.commands.variants",
            "varlens-variants --include-",
            {
                "--include-effect",
                "--include-gene",
                "--include-context",
                "--include-mhc-binding",
                "--include-read-evidence",
                "--include-metadata",
                "--include-failing-variants",
            },
        ),
        (
            "varlens.commands.variants",
            "varlens-variants --mhc-",
            {
                "--mhc-predictor",
                "--mhc-alleles",
                "--mhc-alleles-file",
                "--mhc-peptide-lengths",
            },
        ),
        (
            "varlens.commands.reads",
            "varlens-reads --",
            {
                "--out",
                "--field",
                "--no-standard-fields",
                "--no-sort",
                "--header",
                "--header-set",
                "--locus",
                "--variants",
                "--min-mapping-quality",
                "--version",
                "--verbose",
            },
        ),
        (
            "varlens.commands.allele_support",
            "varlens-allele-support --",
            {
                "--out",
                "--locus",
                "--variants",
                "--reads",
                "--min-mapping-quality",
                "--version",
                "--verbose",
            },
        ),
        (
            "varlens.commands.fusions",
            "varlens-fusions --",
            {
                "--out",
                "--min-fragments",
                "--mhc-predictor",
                "--mhc-alleles",
                "--mhc-alleles-file",
                "--mhc-peptide-lengths",
                "--peptide-lengths",
                "--hla",
                "--version",
            },
        ),
    ],
)
def test_argcomplete_returns_expected_flags_without_heavy_imports(
    tmp_path: Path,
    entry_module: str,
    comp_line: str,
    expected_flags: set[str],
) -> None:
    completions = run_completion_subprocess(tmp_path, entry_module, comp_line)
    assert expected_flags <= completions, (
        f"Expected {expected_flags} in completions {completions}"
    )


@pytest.mark.parametrize(
    "module_name",
    [
        "varlens",
        "varlens.commands",
        "varlens.commands.help",
        "varlens.commands.arguments",
        "varlens.commands.variants",
        "varlens.commands.reads",
        "varlens.commands.allele_support",
        "varlens.commands.fusions",
    ],
)
def test_package_and_parser_imports_avoid_heavy_dependencies(module_name: str) -> None:
    code = (
        "import importlib, sys\n"
        f"importlib.import_module({module_name!r})\n"
        f"loaded = [m for m in {HEAVY_MODULES!r} if m in sys.modules]\n"
        "assert not loaded, f'Heavy modules loaded on import: {loaded}'\n"
    )
    result = subprocess.run(
        [sys.executable, "-c", code],
        cwd=str(REPO_ROOT),
        env=build_subprocess_env(),
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stdout + result.stderr


@pytest.mark.parametrize(
    "entry_module",
    [
        "varlens.commands.variants",
        "varlens.commands.reads",
        "varlens.commands.allele_support",
        "varlens.commands.fusions",
    ],
)
def test_help_and_version_avoid_heavy_dependencies(entry_module: str) -> None:
    code = (
        "import importlib, sys\n"
        f"mod = importlib.import_module({entry_module!r})\n"
        "for flag in ('--help', '--version'):\n"
        "    try:\n"
        "        mod.run([flag])\n"
        "    except SystemExit as exc:\n"
        "        assert exc.code == 0\n"
        f"loaded = [m for m in {HEAVY_MODULES!r} if m in sys.modules]\n"
        "assert not loaded, f'Heavy modules loaded for --help/--version: {loaded}'\n"
    )
    result = subprocess.run(
        [sys.executable, "-c", code],
        cwd=str(REPO_ROOT),
        env=build_subprocess_env(),
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert "usage:" in result.stdout


def test_varlens_lazy_exports_remain_discoverable_and_resolve() -> None:
    code = (
        "import sys\n"
        "import varlens\n"
        "assert set(varlens.__all__) <= set(dir(varlens))\n"
        "assert not ({'pandas', 'numpy', 'pysam'} & sys.modules.keys())\n"
        "from varlens import util, loci_util, read_evidence\n"
        "assert varlens.util is util\n"
        "assert varlens.loci_util is loci_util\n"
        "assert varlens.read_evidence is read_evidence\n"
        "assert not hasattr(varlens, 'nonexistent_varlens_attr')\n"
    )
    result = subprocess.run(
        [sys.executable, "-c", code],
        cwd=str(REPO_ROOT),
        env=build_subprocess_env(),
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stdout + result.stderr


def test_all_console_script_entry_modules_carry_argcomplete_marker() -> None:
    pyproject = tomllib.loads((REPO_ROOT / "pyproject.toml").read_text(encoding="utf-8"))
    scripts = pyproject["project"]["scripts"]
    assert scripts, "Expected console scripts in pyproject.toml"
    for script_name, target in scripts.items():
        module_path = target.split(":", 1)[0].replace(".", "/") + ".py"
        head = (REPO_ROOT / module_path).read_bytes()[:1024]
        assert b"# PYTHON_ARGCOMPLETE_OK" in head, (
            f"Missing '# PYTHON_ARGCOMPLETE_OK' in first 1024 bytes of "
            f"{module_path} (entry point for {script_name})"
        )


@pytest.mark.parametrize("command", [variants, reads, allele_support, fusions])
def test_every_command_uses_help_argument_parser_with_plain_format_help(
    command, monkeypatch: pytest.MonkeyPatch
) -> None:
    assert isinstance(command.parser, HelpArgumentParser)
    plain = command.parser.format_help()
    assert "\033[" not in plain
    assert "``" not in plain

    monkeypatch.delenv("NO_COLOR", raising=False)
    monkeypatch.setenv("TERM", "xterm-256color")

    non_tty = io.StringIO()
    command.parser.print_help(file=non_tty)
    assert non_tty.getvalue() == plain
    assert "\033[" not in non_tty.getvalue()

    tty = FakeTTY()
    command.parser.print_help(file=tty)
    styled = tty.getvalue()
    assert "\033[1musage:\033[0m" in styled
    assert "\033[36m--version\033[0m" in styled
    assert "\033[1;36m" in styled
    # Only normal cyan (36), bold cyan (1;36), and bold (1) are used — no bright codes (90-97).
    assert not re.search(r"\033\[[0-9;]*9[0-7]m", styled)
    assert re.sub(r"\033\[[0-9;]*m", "", styled) == plain


@pytest.mark.parametrize(
    "env_name,env_value",
    [("NO_COLOR", ""), ("NO_COLOR", "1"), ("TERM", "dumb")],
)
def test_help_and_table_styling_respect_no_color_and_dumb_term(
    monkeypatch: pytest.MonkeyPatch, env_name: str, env_value: str
) -> None:
    monkeypatch.delenv("NO_COLOR", raising=False)
    monkeypatch.setenv("TERM", "xterm-256color")
    monkeypatch.setenv(env_name, env_value)

    tty = FakeTTY()
    assert not color_enabled(tty)
    parser = HelpArgumentParser(prog="varlens-test", description="Plain help.")
    parser.add_argument("--example", help="Example ``flag``.")
    parser.print_help(file=tty)
    assert "\033[" not in tty.getvalue()
    assert "\033[" not in style_help(parser.format_help(), tty)
    assert "\033[" not in style_table_header("COL1  COL2\nval1  val2\n", tty)
    assert "\033[" not in style_status("missing", tty, kind="warning")
    assert "\033[" not in style_status("-", tty, kind="absent")


def test_table_header_and_status_styling(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.delenv("NO_COLOR", raising=False)
    monkeypatch.setenv("TERM", "xterm-256color")

    tty = FakeTTY()
    non_tty = io.StringIO()

    table = "CONTIG  STATUS\n22      OK\n"
    styled_table = style_table_header(table, tty)
    assert styled_table == "\033[1;36mCONTIG  STATUS\033[0m\n22      OK\n"
    assert style_table_header(table, non_tty) == table

    assert style_status("OK", tty, kind="normal") == "OK"
    assert style_status("-", tty, kind="absent") == "\033[2m-\033[0m"
    assert style_status("missing", tty, kind="warning") == "\033[33mmissing\033[0m"
    assert style_status("-", non_tty, kind="absent") == "-"
    assert style_status("missing", non_tty, kind="warning") == "missing"

    with pytest.raises(ValueError, match="Unsupported status kind"):
        style_status("OK", tty, kind="green")


def test_help_argument_parser_disables_stdlib_color_on_python_314() -> None:
    parser = HelpArgumentParser(prog="varlens-test")
    if sys.version_info >= (3, 14):
        assert getattr(parser, "color", None) is False


def test_lightweight_argument_specs_match_runtime_modules() -> None:
    from varlens import reads_util
    from varlens.read_evidence.pileup_collection import PileupCollection

    assert list(READ_ATTRIBUTE_NAMES) == list(PileupCollection._READ_ATTRIBUTE_NAMES)
    assert [name for name, _, _ in iter_read_filter_specs()] == list(
        reads_util.READ_FILTERS.keys()
    )
    assert parse_int_list("8,9,10-12") == [8, 9, 10, 11, 12]
    assert parse_mhc_predictor_names("random") == ["random"]
    with pytest.raises(argparse.ArgumentTypeError, match="Unknown predictor"):
        parse_mhc_predictor_names("definitely_not_a_predictor")

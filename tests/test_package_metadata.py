from importlib.metadata import distributions

from varlens.version import __version__


def test_installed_metadata_matches_source_version():
    # A checkout's egg-info can hide stale dist-info in site-packages. Check
    # every visible distribution, not just the first match returned by version().
    installed = list(distributions(name="varlens"))
    assert installed, "Install the checkout with: python -m pip install -e ."
    assert {dist.version for dist in installed} == {__version__}, (
        "Stale Varlens metadata; refresh with: python -m pip install -e ."
    )

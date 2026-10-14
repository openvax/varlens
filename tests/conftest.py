import pytest
from pyensembl import EnsemblRelease


def pytest_addoption(parser):
    parser.addoption(
        "--run-netmhcpan",
        action="store_true",
        help="Run the optional test using an installed NetMHCpan executable",
    )


def pytest_collection_modifyitems(config, items):
    reference_installed = EnsemblRelease(75).installed()
    for item in items:
        if item.get_closest_marker("ensembl") and not reference_installed:
            item.add_marker(
                pytest.mark.skip(
                    reason="Install reference data with: pyensembl install --release 75 --species human",
                )
            )
        if item.get_closest_marker("netmhcpan") and not config.getoption(
            "--run-netmhcpan"
        ):
            item.add_marker(
                pytest.mark.skip(
                    reason="Use --run-netmhcpan to enable the external predictor test"
                )
            )

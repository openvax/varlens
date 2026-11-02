Development
===========

How to run the tests, check the shared OpenVax dependency stack, and build and
publish these docs. This page is for contributors and maintainers.

Tests
-----

To run the tests from a git checkout:

::

    pip install -e '.[test,fusion]'
    pyensembl install --release 75 --species human
    ./lint.sh
    ./test.sh

External NetMHCpan tests are opt-in: run ``./test.sh --run-netmhcpan`` with
NetMHCpan and its model data installed. The other tests do not require it.

CI tests Python 3.11 through 3.14, the pandas 3.0.0 floor and current 3.x,
and reference-dependent annotation in a separate job. A shared-stack job
also installs Varlens, Vaxrank, MHCTools, and Trufflepig together with normal
dependency resolution and runs ``pip check``. To reproduce that check in a
fresh environment:

::

    python -m pip install -e '.[test,fusion]' -r requirements-openvax.txt
    python -m pip check
    ./test.sh

Refresh editable metadata from the existing checkout with ``pip install -e .``
after changing its requirements. The metadata test checks every visible
Varlens distribution against the source version, including stale editable
metadata in site-packages.

Shared environments
-------------------

In a shared environment, resolve all installed OpenVax consumers together
and preserve existing editable checkouts. Varlens 1.1 requires the coordinated
Varcode 11.2, Isovar 1.78 (for fusions), and Topiary 5.109 stack; installations
with Vaxrank need its compatible 3.49 release. Current Topiary and Trufflepig
1.30 require pandas 3; the older Trufflepig 1.28 pandas cap is incompatible
with this stack. Do not force incompatible versions or silently downgrade
another consumer. Repository-local release environments remain separate from
the shared environment.

Documentation
-------------

These docs are built with Sphinx from ``docs/``. The front page includes the
project README, and the command-line reference embeds each command's ``-h``
output:

::

    pip install -e '.[docs]'
    python -m sphinx -W --keep-going -b html docs docs/_build/html

The docs will be written to ``docs/_build/html``. Maintainers with push access
can publish them using ``./deploy-gh-pages.sh``; it builds in a temporary
directory and updates the existing GitHub Pages branch without rewriting its
history. Set ``PYTHON`` to select a particular environment's Python executable.

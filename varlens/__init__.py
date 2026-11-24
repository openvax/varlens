# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Public Varlens package exports, loaded lazily for fast CLI startup."""

from importlib import import_module
from typing import TYPE_CHECKING

from .version import __version__

if TYPE_CHECKING:
    from . import loci_util as loci_util
    from . import read_evidence as read_evidence
    from . import util as util

LAZY_SUBMODULES = frozenset({"loci_util", "read_evidence", "util"})

__all__ = [
    "__version__",
    "loci_util",
    "read_evidence",
    "util",
]


def __getattr__(name: str):
    """Load public submodules on first access to keep CLI startup lightweight."""
    if name in LAZY_SUBMODULES:
        module = import_module(f".{name}", __name__)
        globals()[name] = module
        return module
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__() -> list[str]:
    """Include lazy public exports in interactive discovery."""
    return sorted(set(globals()) | set(__all__))

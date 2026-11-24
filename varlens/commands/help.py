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

"""Readable argparse help and table styling with optional terminal color."""

from __future__ import annotations

import argparse
import os
import re
import shutil
import sys
from typing import Any

OPTION_PATTERN = re.compile(r"(?<![\w-])(--[a-zA-Z][\w-]*|-[a-zA-Z])\b")


def color_enabled(stream: Any) -> bool:
    """Return whether ``stream`` supports optional terminal styling.

    Styling is enabled only when ``stream`` is a TTY, ``NO_COLOR`` is not set,
    and ``TERM`` is not ``"dumb"``.
    """
    return bool(
        getattr(stream, "isatty", lambda: False)()
        and "NO_COLOR" not in os.environ
        and os.environ.get("TERM") != "dumb"
    )


def style_help(text: str, stream: Any) -> str:
    """Style headings, ``usage:``, and option flags on color-capable terminals.

    Uses normal-intensity ANSI cyan (``\\033[36m``) for option flags, bold cyan
    (``\\033[1;36m``) for section headings, and bold (``\\033[1m``) for
    ``usage:``. Everything else stays plain.
    """
    if not color_enabled(stream):
        return text

    lines: list[str] = []
    for line in text.splitlines(keepends=True):
        if line.startswith("usage:"):
            line = OPTION_PATTERN.sub(
                lambda match: "\033[36m" + match[0] + "\033[0m", line
            )
            line = line.replace("usage:", "\033[1musage:\033[0m", 1)
            lines.append(line)
        elif line.rstrip().endswith(":") and not line[0].isspace():
            ending = "\n" if line.endswith("\n") else ""
            lines.append("\033[1;36m" + line.rstrip("\r\n") + "\033[0m" + ending)
        else:
            line = OPTION_PATTERN.sub(
                lambda match: "\033[36m" + match[0] + "\033[0m", line
            )
            lines.append(line)
    return "".join(lines)


def style_table_header(text: str, stream: Any) -> str:
    """Style the first line of an aligned table in bold cyan on color terminals."""
    if not text or not color_enabled(stream):
        return text
    lines = text.splitlines(keepends=True)
    ending = "\n" if lines[0].endswith("\n") else ""
    lines[0] = "\033[1;36m" + lines[0].rstrip("\r\n") + "\033[0m" + ending
    return "".join(lines)


def style_status(text: str, stream: Any, kind: str = "normal") -> str:
    """Style a table status value on color-capable terminals.

    Parameters
    ----------
    text : str
        Status text to display.
    stream : file-like
        Destination stream checked via ``color_enabled``.
    kind : {"normal", "absent", "warning"}, optional
        Status category. ``"absent"`` renders dim (``\\033[2m``),
        ``"warning"`` renders yellow (``\\033[33m``), and ``"normal"``
        (including ``"OK"``) remains plain.
    """
    if kind not in ("normal", "absent", "warning"):
        raise ValueError(
            f"Unsupported status kind {kind!r}; expected 'normal', 'absent', or 'warning'"
        )
    if not text or not color_enabled(stream):
        return text
    if kind == "absent":
        return "\033[2m" + text + "\033[0m"
    if kind == "warning":
        return "\033[33m" + text + "\033[0m"
    return text


class HelpFormatter(argparse.RawDescriptionHelpFormatter):
    """Limit reading width and separate options without adding table borders."""

    def __init__(self, prog: str) -> None:
        width = max(40, min(100, shutil.get_terminal_size().columns - 2))
        super().__init__(prog, max_help_position=32, width=width)

    def _split_lines(self, text: str, width: int) -> list[str]:
        # RST literals are useful in generated docs but noisy in terminal help.
        return super()._split_lines(text.replace("``", ""), width)

    def _format_action(self, action: argparse.Action) -> str:
        text = super()._format_action(action)
        return text + "\n" if text else text


class HelpArgumentParser(argparse.ArgumentParser):
    """Keep ``format_help()`` plain; add color only when printing to a terminal."""

    def __init__(self, *args: Any, **kwargs: Any) -> None:
        kwargs.setdefault("formatter_class", HelpFormatter)
        if sys.version_info >= (3, 14):
            # Apply the same styling on all supported Python versions.
            kwargs["color"] = False
        super().__init__(*args, **kwargs)

    def print_help(self, file: Any = None) -> None:
        stream = sys.stdout if file is None else file
        self._print_message(style_help(self.format_help(), stream), stream)

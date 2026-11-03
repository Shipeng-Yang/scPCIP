"""Start-up banner shown by the command line and by an interactive Python session.

Set SCPCIP_NO_BANNER=1 to switch it off; colour is used only on a terminal and never when NO_COLOR is set.
"""

from __future__ import annotations

import os
import sys

from . import __version__

ART = r"""
          _____   _____ _____ _____
         |  __ \ / ____|_   _|  __ \
 ___  ___| |__) | |      | | | |__) |
/ __|/ __|  ___/| |      | | |  ___/
\__ \ (__| |    | |____ _| |_| |
|___/\___|_|     \_____|_____|_|
"""

TAGLINE = ("single-cell Plant Cell Identity Programs   v{version}\n"
           "root | leaf | flower | seed  -  cross-species, no reference atlas needed\n")


def text(colour: bool = False) -> str:
    art = ART.strip("\n")
    if colour:
        art = f"\033[1;38;5;107m{art}\033[0m"
    return art + "\n\n" + TAGLINE.format(version=__version__)


def enabled() -> bool:
    return os.environ.get("SCPCIP_NO_BANNER", "").strip().lower() not in ("1", "true", "yes")


def show(stream=None) -> None:
    """Print the banner (to stderr by default, so that command output on stdout is unchanged)."""
    stream = sys.stderr if stream is None else stream
    isatty = getattr(stream, "isatty", lambda: False)
    colour = bool(isatty()) and "NO_COLOR" not in os.environ and os.environ.get("TERM") != "dumb"
    print(text(colour=colour), file=stream)


def interactive_session() -> bool:
    if hasattr(sys, "ps1") or sys.flags.interactive:
        return True
    ipython = sys.modules.get("IPython")
    return bool(ipython and getattr(ipython, "get_ipython", lambda: None)())

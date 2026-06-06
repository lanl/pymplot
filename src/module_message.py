"""
Consistent, colored terminal output for pymplot.

All informational output goes to stdout; warnings and errors go to stderr.
Colors are auto-disabled when the stream is not a TTY or when the
NO_COLOR / PYMPLOT_NO_COLOR environment variable is set.
"""

import sys
import os


# ---------------------------------------------------------------------------
# color helpers
# ---------------------------------------------------------------------------

def _tty(stream):
    return getattr(stream, 'isatty', lambda: False)() and \
           not os.environ.get('NO_COLOR') and \
           not os.environ.get('PYMPLOT_NO_COLOR')

_RST = '\033[0m'
_BLD = '\033[1m'
_DIM = '\033[2m'
_RED = '\033[91m'
_YEL = '\033[93m'
_GRN = '\033[92m'
_CYN = '\033[96m'
_BLU = '\033[94m'
_WHT = '\033[97m'


def _fmt(code, text, stream):
    return f'{code}{text}{_RST}' if _tty(stream) else text


def _out(code, text):
    return _fmt(code, text, sys.stdout)

def _err(code, text):
    return _fmt(code, text, sys.stderr)


# ---------------------------------------------------------------------------
# public API
# ---------------------------------------------------------------------------

def msg_input(path, shape, dmin, dmax):
    """Print a summary of a data file that was read."""
    shape_str = ' × '.join(str(s) for s in shape)
    range_str = f'{dmin:.4e}  →  {dmax:.4e}'
    print(_out(_CYN, f'  ← {path}'))
    print(_out(_DIM, f'     dimensions:  {shape_str}'))
    print(_out(_DIM, f'     value range: {range_str}'))


def msg_clip(cmin, cmax):
    """Print the clip / plot range, indented under the last msg_input call."""
    range_str = f'{cmin:.4e}  →  {cmax:.4e}'
    print(_out(_DIM, f'     clip range:  {range_str}'))


def msg_output(path):
    """Print the output file path."""
    from datetime import datetime
    ts = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    print(_out(_GRN,  f'  → {path}') + '   ' + _out(_DIM, ts))


def msg_done():
    """Print a trailing blank line (replaces the bare print() calls)."""
    print()


def warn(msg):
    """Print a yellow warning to stderr."""
    prefix = _err(_YEL, '  ⚠  warning: ')
    print(f'{prefix}{msg}', file=sys.stderr)


def error(msg):
    """Print a red error to stderr (caller is responsible for exit)."""
    prefix = _err(_RED + _BLD, '  ✗  error: ')
    print(f'{prefix}{msg}', file=sys.stderr)


def fatal(msg):
    """Print a red error to stderr and exit with code 1."""
    error(msg)
    sys.exit(1)

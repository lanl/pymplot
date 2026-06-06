"""
Function API for pymplot plotting functions.

Each function takes a numpy array as its first argument and keyword arguments
that map to the corresponding CLI parameters. All parameters have defaults so
minimal syntax is needed:

    showmatrix(data)
    showslice(data, outfile='slices.pdf')
    showcontour(data, d1=0.01, d2=0.01, colormap='gray')
"""

import os
import sys
import numpy as np
import argparse

_src_dir = os.path.dirname(os.path.abspath(__file__))
if _src_dir not in sys.path:
    sys.path.insert(0, _src_dir)


def _default_args(program):
    """Return an argparse Namespace with all defaults for *program*."""
    from module_getarg import getarg
    parser = argparse.ArgumentParser()
    parser = getarg(parser, program)
    return parser.parse_args([])


# arguments internally consumed via eval() — numeric values are stringified
_EVAL_ARGS = {'slice1', 'slice2', 'slice3'}

# argparse arguments that use nargs='+' and therefore need to be lists
_LIST_ARGS = {
    'outfile', 'infile', 'contours', 'ticks1', 'ticks2', 'ticks3',
    'curve', 'curvestyle', 'curvesize', 'curveselect', 'curvewidth',
    'curvecolor', 'curvefacecolor', 'curveedgecolor', 'curveorder',
    'text', 'textstyle', 'textloc', 'textrotation', 'textsize',
    'textcolor', 'textbox', 'textboxedgecolor', 'textboxfacecolor',
    'textboxalpha', 'textboxpad', 'textboxstyle', 'textorder',
    'arrow', 'arrowfacecolor', 'arrowedgecolor', 'arrowlinestyle',
    'arrowstyle', 'arrowconnect', 'arrowwidth', 'arroworder',
    'polygon', 'polygonfacecolor', 'polygonedgecolor', 'polygonalpha',
    'polygonlinestyle', 'polygonlinewidth', 'polygonorder',
    'circle', 'circlefacecolor', 'circleedgecolor', 'circlealpha',
    'circlelinestyle', 'circlelinewidth', 'circleorder',
    'select', 'linestyle', 'linewidth', 'linecolor', 'linealpha',
    'marker', 'markerevery', 'markersize', 'markerfacecolor',
    'markeredgecolor', 'markeralpha', 'markersizemax', 'markersizemin',
    'plotlabel', 'fillwrt', 'fillcolor', 'fillalpha',
    'wigglecolor', 'wigglewidth', 'wigglestyle', 'tracecolor',
    'contourcolor',
}


def _apply_kwargs(args, kwargs):
    """Overwrite *args* attributes with any user-supplied *kwargs*.

    Scalar strings passed for nargs='+' arguments are automatically wrapped
    in a list so callers can write ``outfile='out.pdf'`` instead of
    ``outfile=['out.pdf']``.
    """
    for key, value in kwargs.items():
        if key in _EVAL_ARGS and isinstance(value, (int, float)):
            value = str(value)
        elif key in _LIST_ARGS and isinstance(value, str):
            value = [value]
        setattr(args, key, value)
    return args


# ---------------------------------------------------------------------------
# showmatrix — 2-D image
# ---------------------------------------------------------------------------

def showmatrix(data, background=None, mask=None, ax=None, **kwargs):
    """Plot a 2-D array as an image.

    Parameters
    ----------
    data : ndarray, shape (n1, n2)
    background : ndarray, optional
    mask : ndarray, optional
    **kwargs
        Any CLI parameter accepted by showmatrix (e.g. ``colormap='gray'``,
        ``legend=True``, ``outfile='out.pdf'``, ``d1=0.01``).

    Returns
    -------
    matplotlib.figure.Figure
    """
    from module_io import inject_data, clear_injected
    import showmatrix as _mod

    args = _default_args('matrix')
    args.n1, args.n2 = data.shape

    inject_data('fore', data)
    if background is not None:
        inject_data('back', background)
        args.background = '__injected__'
    if mask is not None:
        inject_data('mask', mask)
        args.mask = '__injected__'

    _apply_kwargs(args, kwargs)
    args._api_mode = True

    try:
        fig = _mod.run(args, ax=ax)
    finally:
        clear_injected()
    return fig


# ---------------------------------------------------------------------------
# showcontour — 2-D contour plot
# ---------------------------------------------------------------------------

def showcontour(data, background=None, mask=None, ax=None, **kwargs):
    """Plot a 2-D array as contours.

    Parameters
    ----------
    data : ndarray, shape (n1, n2)
    background : ndarray, optional
    mask : ndarray, optional
    **kwargs
        Any CLI parameter accepted by showcontour (e.g. ``contourlevel=0.5``,
        ``overlay=True``, ``outfile='out.pdf'``).

    Returns
    -------
    matplotlib.figure.Figure
    """
    from module_io import inject_data, clear_injected
    import showcontour as _mod

    args = _default_args('contour')
    args.n1, args.n2 = data.shape

    inject_data('fore', data)
    if background is not None:
        inject_data('back', background)
        args.background = '__injected__'
    if mask is not None:
        inject_data('mask', mask)
        args.mask = '__injected__'

    _apply_kwargs(args, kwargs)
    args._api_mode = True

    try:
        fig = _mod.run(args, ax=ax)
    finally:
        clear_injected()
    return fig


# ---------------------------------------------------------------------------
# showslice — 3-D orthogonal slice plot
# ---------------------------------------------------------------------------

def showslice(data, background=None, mask=None, **kwargs):
    """Plot three orthogonal slices of a 3-D array.

    Parameters
    ----------
    data : ndarray, shape (n1, n2, n3)
    background : ndarray, optional
    mask : ndarray, optional
    **kwargs
        Any CLI parameter accepted by showslice (e.g. ``slice1='50'``,
        ``outfile='out.pdf'``, ``render='2d'``).

    Returns
    -------
    matplotlib.figure.Figure
    """
    from module_io import inject_data, clear_injected
    import showslice as _mod

    args = _default_args('slice')
    args.n1, args.n2, args.n3 = data.shape

    inject_data('fore', data)
    if background is not None:
        inject_data('back', background)
        args.background = '__injected__'
    if mask is not None:
        inject_data('mask', mask)
        args.mask = '__injected__'

    _apply_kwargs(args, kwargs)
    args._api_mode = True

    try:
        fig = _mod.run(args)
    finally:
        clear_injected()
    return fig


# ---------------------------------------------------------------------------
# showslicon — 3-D orthogonal slice contour plot
# ---------------------------------------------------------------------------

def showslicon(data, background=None, mask=None, **kwargs):
    """Plot contours on three orthogonal slices of a 3-D array.

    Parameters
    ----------
    data : ndarray, shape (n1, n2, n3)
    background : ndarray, optional
    mask : ndarray, optional
    **kwargs
        Any CLI parameter accepted by showslicon.

    Returns
    -------
    matplotlib.figure.Figure
    """
    from module_io import inject_data, clear_injected
    import showslicon as _mod

    args = _default_args('slicon')
    args.n1, args.n2, args.n3 = data.shape

    inject_data('fore', data)
    if background is not None:
        inject_data('back', background)
        args.background = '__injected__'
    if mask is not None:
        inject_data('mask', mask)
        args.mask = '__injected__'

    _apply_kwargs(args, kwargs)
    args._api_mode = True

    try:
        fig = _mod.run(args)
    finally:
        clear_injected()
    return fig


# ---------------------------------------------------------------------------
# showvolume — 3-D isometric volume plot
# ---------------------------------------------------------------------------

def showvolume(data, ax=None, **kwargs):
    """Plot a 3-D array as an isometric volume.

    Parameters
    ----------
    data : ndarray, shape (n1, n2, n3)
    **kwargs
        Any CLI parameter accepted by showvolume (e.g. ``angle='40,10'``,
        ``outfile='out.pdf'``).

    Returns
    -------
    matplotlib.figure.Figure
    """
    from module_io import inject_data, clear_injected
    import showvolume as _mod

    args = _default_args('volume')
    args.n1, args.n2, args.n3 = data.shape

    inject_data('fore', data)

    _apply_kwargs(args, kwargs)
    args._api_mode = True

    try:
        fig = _mod.run(args, ax=ax)
    finally:
        clear_injected()
    return fig


# ---------------------------------------------------------------------------
# showvolcon — 3-D isometric volume contour plot
# ---------------------------------------------------------------------------

def showvolcon(data, ax=None, **kwargs):
    """Plot isometric volume contours of a 3-D array.

    Parameters
    ----------
    data : ndarray, shape (n1, n2, n3)
    **kwargs
        Any CLI parameter accepted by showvolcon.

    Returns
    -------
    matplotlib.figure.Figure
    """
    from module_io import inject_data, clear_injected
    import showvolcon as _mod

    args = _default_args('volcon')
    args.n1, args.n2, args.n3 = data.shape

    inject_data('fore', data)

    _apply_kwargs(args, kwargs)
    args._api_mode = True

    try:
        fig = _mod.run(args, ax=ax)
    finally:
        clear_injected()
    return fig


# ---------------------------------------------------------------------------
# showwiggle — 2-D wiggle trace plot
# ---------------------------------------------------------------------------

def showwiggle(data, background=None, ax=None, **kwargs):
    """Plot a 2-D array as wiggle traces.

    Parameters
    ----------
    data : ndarray, shape (n1, n2)
        Data to plot as wiggles. Axis 2 (columns) are the default trace
        direction; use ``along=2`` to trace along axis 1.
    background : ndarray, optional
    **kwargs
        Any CLI parameter accepted by showwiggle (e.g. ``every=2``,
        ``fill=1``, ``outfile='out.pdf'``).

    Returns
    -------
    matplotlib.figure.Figure
    """
    import tempfile, shutil
    import showwiggle as _mod

    args = _default_args('wiggle')
    n1, n2 = data.shape
    args.n1 = n1
    args.n2 = n2

    # Apply user kwargs first so datatype is set before we write the file
    _apply_kwargs(args, kwargs)

    # Write data to a temporary binary file
    dt_map = {'float': np.float32, 'double': np.float64, 'int': np.int16}
    dt = dt_map.get(args.datatype, np.float32)
    tmp = tempfile.NamedTemporaryFile(suffix='.bin', delete=False)
    tmp_path = tmp.name
    tmp.close()

    # showwiggle reads in (n2, n1) C-order then transposes, so write accordingly
    data.astype(dt).T.flatten(order='C').tofile(tmp_path)
    args.infile = [tmp_path]

    if background is not None:
        tmp_back = tempfile.NamedTemporaryFile(suffix='.bin', delete=False)
        tmp_back_path = tmp_back.name
        tmp_back.close()
        background.astype(dt).T.flatten(order='C').tofile(tmp_back_path)
        args.background = tmp_back_path
    else:
        tmp_back_path = None

    args._api_mode = True

    try:
        fig = _mod.run(args, ax=ax)
    finally:
        os.unlink(tmp_path)
        if tmp_back_path is not None:
            os.unlink(tmp_back_path)
    return fig


# ---------------------------------------------------------------------------
# showcolorbar — standalone colorbar
# ---------------------------------------------------------------------------

def showcolorbar(ax=None, **kwargs):
    """Plot a standalone colorbar.

    Parameters
    ----------
    **kwargs
        Any CLI parameter accepted by showcolorbar (e.g. ``lloc='right'``,
        ``cmin=0``, ``cmax=1``, ``colormap='jet'``, ``outfile='cb.pdf'``).

    Returns
    -------
    matplotlib.figure.Figure  (or None if showcolorbar.run returns None)
    """
    import showcolorbar as _mod

    args = _default_args('colorbar')
    _apply_kwargs(args, kwargs)
    args._api_mode = True

    return _mod.run(args, ax=ax)


# ---------------------------------------------------------------------------
# showgraph — 1-D / scatter graph
# ---------------------------------------------------------------------------

def showgraph(data, ax=None, **kwargs):
    """Plot 1-D or scatter data.

    Parameters
    ----------
    data : ndarray
        * 1-D array ``(n,)`` → ``ptype=1`` (y values only)
        * 2-D array ``(n, 2)`` → ``ptype=2`` (x, y columns)
        * 2-D array ``(n, 3)`` → ``ptype=3`` (x, y, value scatter)
    **kwargs
        Any CLI parameter accepted by showgraph (e.g. ``linecolor='blue'``,
        ``outfile='graph.pdf'``, ``label1='Time (s)'``).

    Returns
    -------
    matplotlib.figure.Figure
    """
    import tempfile
    import showgraph as _mod

    data = np.asarray(data, dtype=float)
    if data.ndim == 1:
        data = data.reshape(-1, 1)   # column vector for ptype=1

    args = _default_args('graph')

    # infer ptype from shape if not user-supplied
    if 'ptype' not in kwargs:
        ncols = data.shape[1]
        if ncols == 1:
            args.plottype = 1
        elif ncols == 2:
            args.plottype = 2
        else:
            args.plottype = 3

    # showgraph reads ASCII or binary; easiest API path is ASCII
    tmp = tempfile.NamedTemporaryFile(suffix='.txt', delete=False, mode='w')
    tmp_path = tmp.name
    np.savetxt(tmp, data)
    tmp.close()

    args.infile = [tmp_path]
    args.filetype = 'ascii'

    # n1 for graph is a string list (number of rows)
    args.n1 = [str(data.shape[0])]

    _apply_kwargs(args, kwargs)
    args._api_mode = True

    try:
        fig = _mod.run(args, ax=ax)
    finally:
        os.unlink(tmp_path)
    return fig

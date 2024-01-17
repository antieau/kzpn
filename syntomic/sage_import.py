# sage_import.py
import imp
import inspect
import os
import sys

import sage.all


def sage_import(modname, fromlist=None, namespace=None):
    """
    Import a .sage module from the filename <modname>.sage

    Returns the resulting Python module.  If ``fromlist`` is given, returns
    just those members of the module into the global namespace where the
    function was called, or the given namespace.
    """

    filename = modname + ".sage"

    for path in sys.path:
        modpath = os.path.join(path, filename)
        if os.path.isfile(modpath):
            break
    else:
        raise ImportError('no file {} on sys.path'.format(filename))

    with open(modpath) as fobj:
        code = sage.all.preparse(fobj.read())
        mod = imp.new_module(modname)
        mod.__file__ = modpath
        # Fill with all the default Sage globals
        # We could just do a dict.update but we want to exclude dunder
        # and private attributes I guess
        for k, v in sage.all.__dict__.items():
            if not k.startswith('_'):
                mod.__dict__[k] = v

        exec(code,mod.__dict__,mod.__dict__)

    if namespace is None:
        namespace = inspect.currentframe().f_back.f_globals


    if fromlist is not None:
        # First check that each name in fromlist exists before adding
        # any of them to the given namespace.
        for name in fromlist:
            if name not in mod.__dict__:
                raise ImportError('cannot import name {!r} from {}'.format(
                     name, filename))

        for name in fromlist:
            namespace[name] = mod.__dict__[name]
    else:
        namespace[modname] = mod

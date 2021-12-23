from .util import *
from .efoscphotredudef import *
from .efoscfastspecdef import *
from .sofiphotredudef import *
from .efoscspec1Ddef import *
from .efoscspec2Ddef import *
from .sofispec1Ddef import *
from .sofispec2Ddef import *
from .efoscastrodef import *
from .sqlcl import *
from .efosccalibdef import *
from .soficalibdef import *
from . import cosmics

__version__ = "unknown"
try:
    from _version import __version__
except ImportError:
    # We're running in a tree that doesn't have a _version.py, so we don't know what our version is.
    pass

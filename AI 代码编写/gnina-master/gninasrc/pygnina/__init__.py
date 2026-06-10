from .pygnina import *
import sys

from importlib import metadata

try:
  __version__ = metadata.version('pygnina')
except:
  __version__ = 'dev'

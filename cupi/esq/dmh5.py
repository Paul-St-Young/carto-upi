import os
import h5py
import numpy as np

def gen_grid(fp, path):
  """Generate numerical grid based on h5 "Grid" specifications.

  Expect grid specifications: type, npts, start, end.
  H5 file can be generated using pimc++ dmparser.py.

  Args:
    fp (h5py.File): handle to h5file
    path (str): path to grid specifiers e.g. '/Ukj0/Grid'

  Return:
    np.array: 1D array of grid points

  Raises:
    RuntimeError: path name does not end with 'Grid'
    NotImplementedError: unrecognized grid type

  Examples:
    >>> fh5 = 'e-e.sq.h5'
    >>> fp = h5py.File(fh5, 'r')
    >>> gpath = 'Ukj0/Grid'
    >>> myx = gen_grid(fp, gpath)
  """
  name = os.path.basename(fp[path].name)
  if name != 'Grid':
    raise RuntimeError('%s must end with "Grid"' % name)
  gtype = fp['%s/Type' % path].value[0]
  if gtype != 'LINEAR':
    raise NotImplementedError('%s grid type not supported' % gtype)
  rmin = fp['%s/Start' % path].value[0]
  rmax = fp['%s/End' % path].value[0]
  nr = fp['%s/NumGridPoints' % path].value[0]
  return np.linspace(rmin, rmax, nr)

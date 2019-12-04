import os
import h5py
import numpy as np

def get_grid(fh5):
  gpath = 'Ukj0/Grid'
  fp = h5py.File(fh5, 'r')
  grid = gen_grid(fp, gpath)
  fp.close()
  return grid

def get_data(fh5, name='Ukj0'):
  fp = h5py.File(fh5, 'r')
  data = fp['%s/Data' % name][()]
  fp.close()
  return data

def show_sampling(fh5, nsamp=2):
  import matplotlib.pyplot as plt
  grid = get_grid(fh5)
  nr = len(grid)
  sampl = []
  for isamp in range(nsamp):
    name = 'Sampling_%d' % isamp
    samp = get_data(fh5, name)
    sampl.append(samp)
  nr1, nmover, ntau = sampl[0].shape
  assert nr1 == nr

  fig, ax_arr = plt.subplots(nsamp, nmover, sharex=True)
  ax_arr[0, 0].set_ylabel(r'$d\tilde{U}/dR_m$')
  ax_arr[1, 0].set_ylabel(r'$d^2\tilde{U}/dR_m^2$')
  for imover in range(nmover):
    for isamp in range(nsamp):
      ax = ax_arr[isamp, imover]
      ax.plot(grid, sampl[isamp][:, imover, 0])
  for icol in range(nmover):
    ax_arr[-1, icol].set_xlabel('x')
  fig.tight_layout()
  plt.show()

def show_all_ukj(fh5, orders=range(1, 4), itau=0, xmin=0.05, xmax=None):
  import matplotlib.pyplot as plt
  # get raw data
  grid = get_grid(fh5)
  nr = len(grid)
  ukjl = []
  dukjl = []
  for norder in orders:
    name = 'Ukj%d' % (norder-1)
    ukj = get_data(fh5, name)
    assert len(ukj) == nr
    ukjl.append(ukj)
    name = 'dUkjdBeta%d' % (norder-1)
    dukj = get_data(fh5, name)
    assert len(dukj) == nr
    dukjl.append(dukj)
  # decide plotting range from data
  if xmax is None:
    xmax = grid.max()/3.
  sel = (xmin < grid) & (grid < xmax)
  imid = min(1, len(orders)-1)

  ymult = 1.2  # show a big more than the max
  myy = ukjl[imid]
  ymin = ymult*myy[sel].min()
  ymax = ymult*myy[sel].max()

  mydy = dukjl[imid]
  dymin = ymult*mydy[sel].min()
  dymax = ymult*mydy[sel].max()

  nrow = 2  # U and dU/dBeta
  ncol = len(orders)
  fig, ax_arr = plt.subplots(nrow, ncol, sharex=True)
  ax_arr[0, 0].set_xlim(xmin, xmax)  # x axis is shared
  for iax, norder in enumerate(orders):
    # first plot ukj
    ax = ax_arr[0, iax]
    ax.set_title('k = %d' % norder)
    ax.set_ylim(ymin, ymax)
    ukj = ukjl[iax]
    nr, nj, ntau = ukj.shape
    if ntau > 1:
      msg = 'ukj printed at %d temperatures' % ntau
      msg += ' use itau to access higher temperatures'
      print(msg)
    for j in range(nj):
      ax.plot(grid, ukj[:, j, itau])
    # second plot dukj/dbeta
    ax = ax_arr[1, iax]
    ax.set_ylim(dymin, dymax)
    for j in range(nj):
      ax.plot(grid, dukj[:, j, itau], label='j=%d' % j)
  # style the axes
  for icol in range(1, ncol):
    for irow in range(nrow):
      ax_arr[irow, icol].set_yticks([])
  for icol in range(ncol):
    ax_arr[-1, icol].set_xlabel('x')
  ax_arr[0, 0].set_ylabel('ukj')
  ax_arr[1, 0].set_ylabel('dukj/dbeta')
  ax_arr[-1, -1].legend()
  fig.subplots_adjust(wspace=0, hspace=0)
  plt.show()

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
    RuntimeError: grid type is not 'LINEAR' or 'LOG'

  Examples:
    >>> fh5 = 'e-e.sq.h5'
    >>> fp = h5py.File(fh5, 'r')
    >>> gpath = 'Ukj0/Grid'
    >>> myx = gen_grid(fp, gpath)
  """
  name = os.path.basename(fp[path].name)
  if name != 'Grid':
    raise RuntimeError('%s must end with "Grid"' % name)
  gtype = fp['%s/Type' % path][()][0]
  rmin = fp['%s/Start' % path][()][0]
  rmax = fp['%s/End' % path][()][0]
  nr = fp['%s/NumGridPoints' % path][()][0]
  if gtype == 'LINEAR':
    grid = np.linspace(rmin, rmax, nr)
  elif gtype == 'LOG':
    grid = np.exp(np.linspace(np.log(rmin), np.log(rmax), nr))
  else:
    raise RuntimeError('unknown grid type %s' % gtype)
  return grid

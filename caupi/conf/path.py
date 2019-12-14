import numpy as np

def select_slices(path, slices):
  iax = 2
  nslices = path.shape[iax]
  sel = np.zeros(nslices, dtype=bool)
  for i in slices:
    sel[i] = True
  return path[:, :, sel]

def select_particles(path, parts):
  iax = 1
  nparts = path.shape[iax]
  sel = np.zeros(nparts, dtype=bool)
  for i in parts:
    sel[i] = True
  return path[:, sel, :]

def mic_path(path, box):
  ndim, nparts, nslices = path.shape
  path1 = np.zeros([ndim, nparts, nslices])
  for ipart in range(nparts):
    # build path
    cur = path[:, ipart, 0]
    path1[:, ipart, 0] = cur
    for islice in range(1, nslices):
      prev = cur
      dr = path[:, ipart, islice] - prev
      # minimum image convention (mic)
      nint = np.around(dr/box)
      dr -= nint*box
      cur = prev + dr
      path1[:, ipart, islice] = cur
  return path1

def show_path2d(ax, path, box, **kwargs):
  """Show path in 2 spatial dimensions as an elevator in 3D

  Args:
    ax (plt.Axes3D): 3D plot
    path (np.array): path, shape (ndim, nparts, nslices), ndim must be 2
    box (np.array): rectangular box, shape (ndim,)
  Return:
    lines (list): a list of Line3D objects
  Example:
    >>> ndim, nparts, nslices = path.shape
    >>> parts = np.arange(0, nparts, 16)  # every 16th particle
    >>> slices = np.arange(0, nslices, 64)  # every 64th time slice
    >>> path = select_slices(path, slices)
    >>> path = select_slices(path, parts)
    >>> fig, ax = figax3d()
    >>> show_path2d(ax, path, box)
    >>> plt.show()
  """
  marker = kwargs.pop('marker', '.')
  ndim, nparts, nslices = path.shape
  if ndim != 2:
    raise RuntimeError('ndim %d != 2' % ndim)
  path1 = mic_path(path, box)
  lines = []
  for ipart in range(nparts):
    x, y = path1[:, ipart, :]
    z = np.arange(nslices)
    line = ax.plot(x, y, z, marker=marker, **kwargs)
    lines.append(line[0])
  return lines

class Path:
  def __init__(self, path, iperm):
    ndim, nparts, nslices = path.shape
    if len(iperm) != nparts:
      msg = 'wrong length %d of permutation index' % len(iperm)
      msg += ' for %d particles' % nparts
      raise RuntimeError(msg)
    self.path = path
    self.iperm = iperm
    self._ndim = ndim
    self._nparts = nparts
    self._nslices = nslices
  def get_ndim(self):
    return self._ndim
  def get_nparts(self):
    return self._nparts
  def get_nslices(self):
    return self._nslices

def find_cycles(path, check=True):
  # get data from class
  ndim = path.get_ndim()
  nparts = path.get_nparts()
  nslices = path.get_nslices()
  iperm = path.iperm
  # find cycles
  visited = np.zeros(nparts, dtype=bool)
  cycles = []
  for ipart in range(nparts):
    visited[ipart] = True
    this_cycle = [ipart]
    jpart = iperm[ipart]
    for inext in range(nparts):
      if jpart == ipart:
        cycles.append(this_cycle)
        break
      if visited[jpart]:
        continue
      this_cycle.append(jpart)
      visited[jpart] = True
      jpart = path.iperm[jpart]
  if check:
    ncycles = [len(c) for c in cycles]
    nparts1 = sum(ncycles)
    if nparts1 != nparts:
      msg = 'got %d/%d particles after find_cycles' % (nparts1, nparts)
      raise RuntimeError(msg)
  return cycles

def show_cycle(ax, path, cycle, slices, **kwargs):
  markers = ['^', 's', 'd', 'p', 'h', '8', 'P', 's']
  nmarker = len(markers)
  ncycle = len(cycle)
  if nmarker < ncycle:
    raise RuntimeError('not enough markers %d/%d' % (nmarker, ncycle))
  for j, ipart in enumerate(cycle):
    x, y = path[:, ipart, slices]
    z = np.arange(len(x))
    line = ax.plot(x, y, z, **kwargs)
    myc = line[0].get_color()
    # plot the first slice of the connecting particle
    jpart = cycle[(j+1) % ncycle]
    x1, y1 = path[:, jpart, slices[0]]
    myx = np.array([x[-1], x1])
    myy = np.array([y[-1], y1])
    myz = np.array([z[-1], z[-1]])
    alpha = 0.4
    # first connect
    ax.plot(myx, myy, myz, c=myc, alpha=alpha, lw=2)
    # mark connection source and target
    m0 = markers[j]
    m1 = markers[(j+1) % ncycle]
    ax.plot([x[0]], [y[0]], [z[0]],
            ls='', marker=m0, c=myc)
    ax.plot([myx[0]], [myy[0]], [myz[0]],
            ls='', marker=m0, c=myc)
    ax.plot([myx[1]], [myy[1]], [myz[1]],
            ls='', marker=m1, c=myc, alpha=alpha)

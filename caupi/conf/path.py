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

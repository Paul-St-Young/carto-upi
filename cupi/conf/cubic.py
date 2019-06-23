import numpy as np

def pos_in_box(pos, lbox):
  """Positions in [-lbox/2, lbox/2)

  Args:
    pos (np.array): positions in open BC
    lbox (float): cubic box side length
  Return:
    np.array: positions in box centered at 0
  """
  return (pos+lbox/2.) % lbox - lbox/2.

def cubic_pos(nx, ndim=3):
  """Initialize simple cubic lattice at integer points

  Args:
    nx (int): number of points along each dimension
    ndim (int): number of spatial dimensions
  Return:
    np.array: simple cubic lattice positions, shape (nx**3, ndim)
  """
  nl = [range(nx)]*ndim
  pos = np.stack(
    np.meshgrid(*nl, indexing='ij'),
    axis=-1
  ).reshape(-1, ndim)
  return pos

def sc(npart, lbox, ndim=3):
  """Simple cubic positions

  Args:
    npart (int): number of particles
    lbox (float): side length of cubic box
  Return:
    np.array: simple cubic lattice
  """
  nx = int(round(npart**(1./3)))
  if nx**3 != npart:
    msg = 'SC npart %d is not a perfect cube' % npart
    raise RuntimeError(msg)
  pos = lbox*cubic_pos(nx)/nx
  assert len(pos) == npart
  return pos_in_box(pos, lbox)

def bcc(npart, lbox, ndim=3):
  """Body-centered cubic positions

  Args:
    npart (int): number of particles
    lbox (float): side length of cubic box
  Return:
    np.array: body-centered cubic lattice
  """
  nx = int(round((npart/2.)**(1./3)))
  if 2*nx**3 != npart:
    msg = 'BCC npart %d is not 2* a perfect cube' % npart
    raise RuntimeError(msg)
  nvecs = cubic_pos(nx)
  pos = lbox*np.concatenate([
    nvecs/float(nx), (nvecs+0.5)/nx
  ], axis=0)
  assert len(pos) == npart
  return pos_in_box(pos, lbox)

def fcc(npart, lbox, ndim=3):
  """Fody-centered cubic positions

  Args:
    npart (int): number of particles
    lbox (float): side length of cubic box
  Return:
    np.array: face-centered cubic lattice
  """
  nx = int(round((npart/4.)**(1./3)))
  if 4*nx**3 != npart:
    msg = 'FCC npart %d is not 4* a perfect cube' % npart
    raise RuntimeError(msg)
  nvecs = cubic_pos(nx)
  posul = [nvecs/float(nx)]
  for idim in range(ndim):
    shift = 0.5*np.ones(ndim)
    shift[idim] = 0
    posu = (nvecs+shift)/float(nx)
    posul.append(posu)
  pos = lbox*np.concatenate(posul, axis=0)
  assert len(pos) == npart
  return pos_in_box(pos, lbox)

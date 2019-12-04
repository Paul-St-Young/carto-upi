import numpy as np

def extract_path(frs, ndim, nparts, nslices):
  """Extract saved path from restart file

  Args:
    frs (str): restart file name qid.rs
    ndim (int): number of spatial dimensions
    nparts (int): number of particles in total, i.e. all species
    nslices (int): number of time slices
  Return:
    np.array: path of shape (ndim, nparts, nslices)
  Example:
    >>> path = extract_path('h32.rs', 3, 64, 8)  # 32 hydrogens
    >>> pos = path[:, 32:, 0].T  # get reference slice protons
  """
  import struct
  with open(frs, 'rb') as fp:
    content = fp.read()

  # read header character lengths from pickup
  nqdt = 9
  nqtime = 26
  mhead = 4  # what are the first 4 bits of .rs file?
  nhead = nqdt+nqtime  # length of header
  header = content[mhead:mhead+nhead]

  # read path
  ndouble = ndim*nparts*nslices
  istart = mhead+nhead+8
  iend = istart+ndouble*8
  #  check zero paddings around path
  zpad0 = struct.unpack('d', content[istart-8: istart])
  assert np.isclose(zpad0, 0)
  zpad1 = struct.unpack('d', content[iend:iend+8])
  assert np.isclose(zpad1, 0)
  data = struct.unpack('d'*ndouble, content[istart: iend])

  path = np.array(data).reshape([ndim, nparts, nslices], order='F')
  return path

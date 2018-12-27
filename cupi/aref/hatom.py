import numpy as np

def hatom_energy(beta, nmax=4000):
  """Calculate the energy of a hydrogen atom at inverse temperature beta
  assume clamped proton and non-relativistic kinetic energy

  Args:
    beta (float): inverse temperature in hartree atomic units
    nmax (int): number of excited states to sum over
  Return:
    float: total energy
  Example:
    >>> hatom_energy(100.)
    -0.5
  """
  # calculate energy levels
  myn = np.arange(1, nmax)
  en = -1./(2*myn**2)
  # calculate Boltzmann factors
  bfac = np.exp(-beta*en)
  part = bfac.sum()
  return np.dot(bfac, en)/part

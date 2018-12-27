import numpy as np

def freep_energy(beta, lam=0.5):
  """Calculate the energy of a free particle at inverse temperature beta

  Args:
    beta (float): inverse temperature in hartree atomic units
    lam (float, optional): quantumness \hbar^2/(2m), default 0.5 i.e. electron
  Return:
    float: total energy = kinetic energy
  """
  return 3./(4*beta*lam)

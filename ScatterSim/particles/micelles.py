from .base import NanoObject 
import numpy as np
from sasmodels.special import sas_3j1x_x, sas_sinx_x

class SphericalMicelle(NanoObject):
    ''' A spherical micelle defined by radius and SLD.'''

    def __init__(self, pargs={}):
        super().__init__(pargs=pargs)

        self.radius_core = self.pargs.get("radius_core", 4.0)
        self.rho_solvent = self.pargs.get("rho_solvent", 2.0)
        self.rho_core = self.pargs.get("rho_core", 1.0)
        self.rho_corona = self.pargs.get("rho_corona", 1.0)
        self.v_core = self.pargs.get("v_core", 4.0)
        self.v_corona = self.pargs.get("v_corona", 4.0)
        self.rg = self.pargs.get("rg", 1.0)
        self.d_penetration = self.pargs.get("d_penetration", 1.0)
        self.x_solv = self.pargs.get("x_solv", 0.0)

    def form_factor(self, qvec):
        ''' Compute the form factor of a spherical micelle. '''
        qx, qy, qz = qvec
        q = np.sqrt(qx**2 + qy**2 + qz**2)
        phase = self.get_phase(qvec)

        n_aggreg = ((1-self.x_solv)*(4/3)*np.pi*(self.radius_core**3))/self.v_core

        beta_core = self.v_core * (self.rho_core - self.rho_solvent)
        beta_corona = self.v_corona * (self.rho_corona - self.rho_solvent)

        # Self-correlation term of the core
        bes_core = sas_3j1x_x(q*self.radius_core)
        term1 = np.power(n_aggreg*beta_core*bes_core, 2)

        # Self-correlation term of the chains
        qrg2 = np.power(q*self.rg, 2)
        debye_chain = 2.0*(np.vectorize(np.expm1)(-qrg2)+qrg2)/(qrg2**2) 
        debye_chain[qrg2==0.0] = 1.0
        term2 = n_aggreg * beta_corona * beta_corona * debye_chain

        # Interference cross-term between core and chains
        chain_ampl = -np.vectorize(np.expm1)(-qrg2)/qrg2
        chain_ampl[qrg2==0.0] =  1.0 
        bes_corona = sas_sinx_x(q*(self.radius_core + (self.d_penetration * self.rg)))
        term3 = 2.0 * n_aggreg * n_aggreg * beta_core * beta_corona * bes_core * chain_ampl * bes_corona

        # Interference cross-term between chains
        term4 = n_aggreg * (n_aggreg - 1.0)* np.power(beta_corona * chain_ampl * bes_corona, 2)

        # I(q)_micelle : Sum of 4 terms computed above
        Fq = term1 + term2 + term3 + term4

        return Fq*phase
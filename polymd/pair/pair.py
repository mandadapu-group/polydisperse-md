# Copyright (c) 2009-2022 The Regents of the University of Michigan
# This file is part of the HOOMD-blue project, released under the BSD 3-Clause License.

#It's possible that I only need to bring in hoomd's pair python script for this 
import copy
import warnings

import hoomd
from hoomd.md import _md
from hoomd.md import force

from hoomd.polymd import _polymd
import hoomd.md.pair as md_pair

from hoomd.data.parameterdicts import ParameterDict, TypeParameterDict
from hoomd.data.typeparam import TypeParameter
import numpy as np
from hoomd.data.typeconverter import OnlyFrom, nonnegative_real

class Polydisperse(md_pair.Pair):
    r"""Polydisperse pair force.
    TO DO: a Better Write-Up
    """
    _cpp_class_name = "PotentialPairPolydisperse"

    def __init__(self,
                 nlist,
                 default_r_cut=None,
                 default_r_on=0.,
                 mode='none',
                 tail_correction=False):
        super().__init__(nlist, default_r_cut, default_r_on, mode)
        params = TypeParameter(
            'params', 'particle_types',
            TypeParameterDict(v0=float,rcut=float,eps=float, m_expnt=int, n_expnt=int, c0=float, c1=float, c2=float, len_keys=8))
       
        #Compute the constants here and do it once!
        c0 = 0.125*(-(2 + m_expnt)*(4 + m_expnt)/rcut**m_expnt + (2 + n_expnt)*(4 + n_expnt)/rcut**n_expnt)*v0;
        c1 = 0.25/rcut**(2 + m_expnt + n_expnt)*(-n_expnt*(4 + n_expnt)*rcut**m_expnt*v0 + m_expnt*(4 + m_expnt)*rcut**n_expnt*v0);
        c2 = 0.125/rcut**(4 + m_expnt + n_expnt)*(n_expnt*(2 + n_expnt)*rcut**m_expnt*v0 - m_expnt*(2 + m_expnt)*rcut**n_expnt*v0);
        params.c0 = c0
        params.c1 = c1
        params.c2 = c2
        
        self._add_typeparam(params)
        self._param_dict.update(
            ParameterDict(tail_correction=bool(tail_correction)))

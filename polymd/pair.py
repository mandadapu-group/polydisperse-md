# Copyright (c) 2009-2019 The Regents of the University of Michigan
# This file is part of the HOOMD-blue project, released under the BSD 3-Clause License.

R""" Pair potentials.

Generally, pair forces are short range and are summed over all non-bonded particles
within a certain cutoff radius of each particle. Any number of pair forces
can be defined in a single simulation. The net force on each particle due to
all types of pair forces is summed.

Pair forces require that parameters be set for each unique type pair. Coefficients
are set through the aid of the :py:class:`coeff` class. To set these coefficients, specify
a pair force and save it in a variable::

    my_force = pair.some_pair_force(arguments...)

Then the coefficients can be set using the saved variable::

    my_force.pair_coeff.set('A', 'A', epsilon=1.0, sigma=1.0)
    my_force.pair_coeff.set('A', 'B', epsilon=1.0, sigma=2.0)
    my_force.pair_coeff.set('B', 'B', epsilon=2.0, sigma=1.0)

This example set the parameters *epsilon* and *sigma*
(which are used in :py:class:`lj`). Different pair forces require that different
coefficients are set. Check the documentation of each to see the definition
of the coefficients.
"""
#It's possible that I only need to bring in hoomd's pair python script for this 
from hoomd.md import force;
from hoomd.md import nlist as nl # to avoid naming conflicts
import hoomd;

import math;
import sys;

#from collections import OrderedDict
from hoomd.polymd import _polymd
from hoomd.md import _md
from hoomd import _hoomd
import hoomd.md.pair as md_pair

class lj_plugin(md_pair.pair):
    R""" Lennard-Jones pair potential.

    Args:
        r_cut (float): Default cutoff radius (in distance units).
        nlist (:py:mod:`hoomd.md.nlist`): Neighbor list
        name (str): Name of the force instance.

    :py:class:`lj` specifies that a Lennard-Jones pair potential should be applied between every
    non-excluded particle pair in the simulation.

    .. math::
        :nowrap:

        \begin{eqnarray*}
        V_{\mathrm{LJ}}(r)  = & 4 \varepsilon \left[ \left( \frac{\sigma}{r} \right)^{12} -
                          \alpha \left( \frac{\sigma}{r} \right)^{6} \right] & r < r_{\mathrm{cut}} \\
                            = & 0 & r \ge r_{\mathrm{cut}} \\
        \end{eqnarray*}

    See :py:class:`pair` for details on how forces are calculated and the available energy shifting and smoothing modes.
    Use :py:meth:`pair_coeff.set <coeff.set>` to set potential coefficients.

    The following coefficients must be set per unique pair of particle types:

    - :math:`\varepsilon` - *epsilon* (in energy units)
    - :math:`\sigma` - *sigma* (in distance units)
    - :math:`\alpha` - *alpha* (unitless) - *optional*: defaults to 1.0
    - :math:`r_{\mathrm{cut}}` - *r_cut* (in distance units)
      - *optional*: defaults to the global r_cut specified in the pair command
    - :math:`r_{\mathrm{on}}`- *r_on* (in distance units)
      - *optional*: defaults to the global r_cut specified in the pair command

    Example::

        nl = nlist.cell()
        lj = pair.lj(r_cut=3.0, nlist=nl)
        lj.pair_coeff.set('A', 'A', epsilon=1.0, sigma=1.0)
        lj.pair_coeff.set('A', 'B', epsilon=2.0, sigma=1.0, alpha=0.5, r_cut=3.0, r_on=2.0);
        lj.pair_coeff.set('B', 'B', epsilon=1.0, sigma=1.0, r_cut=2**(1.0/6.0), r_on=2.0);
        lj.pair_coeff.set(['A', 'B'], ['C', 'D'], epsilon=1.5, sigma=2.0)

    """
    def __init__(self, r_cut, nlist, name=None):
        hoomd.util.print_status_line();

        # tell the base class how we operate

        # initialize the base class
        md_pair.pair.__init__(self, r_cut, nlist, name);

        # create the c++ mirror class
        if not hoomd.context.exec_conf.isCUDAEnabled():
            self.cpp_force = _polymd.PotentialPairLJPlugin(hoomd.context.current.system_definition, self.nlist.cpp_nlist, self.name);
            self.cpp_class = _polymd.PotentialPairLJPlugin;
        else:
            self.nlist.cpp_nlist.setStorageMode(_md.NeighborList.storageMode.full);
            self.cpp_force = _polymd.PotentialPairLJPluginGPU(hoomd.context.current.system_definition, self.nlist.cpp_nlist, self.name);
            self.cpp_class = _polymd.PotentialPairLJPluginGPU;

        hoomd.context.current.system.addCompute(self.cpp_force, self.force_name);

        # setup the coefficient options
        self.required_coeffs = ['epsilon', 'sigma', 'alpha'];
        self.pair_coeff.set_default_coeff('alpha', 1.0);

    def process_coeff(self, coeff):
        epsilon = coeff['epsilon'];
        sigma = coeff['sigma'];
        alpha = coeff['alpha'];

        lj1 = 4.0 * epsilon * math.pow(sigma, 12.0);
        #lj2 = alpha * 4.0 * epsilon #* math.pow(sigma, 6.0);
        #return _hoomd.make_scalar2(lj1, lj2);
        #lj1 = 4.0 * epsilon #* math.pow(sigma, 12.0);
        lj2 = alpha * 4.0 * epsilon * math.pow(sigma, 6.0);
        #return lj1;#_hoomd.make_scalar2(lj1);
        return _hoomd.make_scalar2(lj1,lj2);

class force_shifted_lj_plugin(md_pair.pair):
    R""" Lennard-Jones pair potential.

    Args:
        r_cut (float): Default cutoff radius (in distance units).
        nlist (:py:mod:`hoomd.md.nlist`): Neighbor list
        name (str): Name of the force instance.

    :py:class:`lj` specifies that a Lennard-Jones pair potential should be applied between every
    non-excluded particle pair in the simulation.

    .. math::
        :nowrap:

        \begin{eqnarray*}
        V_{\mathrm{LJ}}(r)  = & 4 \varepsilon \left[ \left( \frac{\sigma}{r} \right)^{12} -
                          \alpha \left( \frac{\sigma}{r} \right)^{6} \right] & r < r_{\mathrm{cut}} \\
                            = & 0 & r \ge r_{\mathrm{cut}} \\
        \end{eqnarray*}

    See :py:class:`pair` for details on how forces are calculated and the available energy shifting and smoothing modes.
    Use :py:meth:`pair_coeff.set <coeff.set>` to set potential coefficients.

    The following coefficients must be set per unique pair of particle types:

    - :math:`\varepsilon` - *epsilon* (in energy units)
    - :math:`\sigma` - *sigma* (in distance units)
    - :math:`\alpha` - *alpha* (unitless) - *optional*: defaults to 1.0
    - :math:`r_{\mathrm{cut}}` - *r_cut* (in distance units)
      - *optional*: defaults to the global r_cut specified in the pair command
    - :math:`r_{\mathrm{on}}`- *r_on* (in distance units)
      - *optional*: defaults to the global r_cut specified in the pair command

    Example::

        nl = nlist.cell()
        lj = pair.lj(r_cut=3.0, nlist=nl)
        lj.pair_coeff.set('A', 'A', epsilon=1.0, sigma=1.0)
        lj.pair_coeff.set('A', 'B', epsilon=2.0, sigma=1.0, alpha=0.5, r_cut=3.0, r_on=2.0);
        lj.pair_coeff.set('B', 'B', epsilon=1.0, sigma=1.0, r_cut=2**(1.0/6.0), r_on=2.0);
        lj.pair_coeff.set(['A', 'B'], ['C', 'D'], epsilon=1.5, sigma=2.0)

    """
    def __init__(self, r_cut, nlist, name=None):
        hoomd.util.print_status_line();

        # tell the base class how we operate

        # initialize the base class
        md_pair.pair.__init__(self, r_cut, nlist, name);

        # create the c++ mirror class
        if not hoomd.context.exec_conf.isCUDAEnabled():
            self.cpp_force = _polymd.PotentialPairForceShiftedLJPlugin(hoomd.context.current.system_definition, self.nlist.cpp_nlist, self.name);
            self.cpp_class = _polymd.PotentialPairForceShiftedLJPlugin;
        else:
            self.nlist.cpp_nlist.setStorageMode(_md.NeighborList.storageMode.full);
            self.cpp_force = _polymd.PotentialPairForceShiftedLJPluginGPU(hoomd.context.current.system_definition, self.nlist.cpp_nlist, self.name);
            self.cpp_class = _polymd.PotentialPairForceShiftedLJPluginGPU;

        hoomd.context.current.system.addCompute(self.cpp_force, self.force_name);

        # setup the coefficient options
        self.required_coeffs = ['epsilon', 'sigma', 'alpha'];
        self.pair_coeff.set_default_coeff('alpha', 1.0);

    def process_coeff(self, coeff):
        epsilon = coeff['epsilon'];
        sigma = coeff['sigma'];
        alpha = coeff['alpha'];

        lj1 = 4.0 * epsilon * math.pow(sigma, 12.0);
        #lj2 = alpha * 4.0 * epsilon #* math.pow(sigma, 6.0);
        #return _hoomd.make_scalar2(lj1, lj2);
        #lj1 = 4.0 * epsilon #* math.pow(sigma, 12.0);
        lj2 = alpha * 4.0 * epsilon * math.pow(sigma, 6.0);
        #return lj1;#_hoomd.make_scalar2(lj1);
        return _hoomd.make_scalar2(lj1,lj2);

class polydisperse(md_pair.pair):
    R""" Polydisperse's custom pair potential.

    """
    def __init__(self, r_cut, nlist, model,name=None, d_max = None):
        hoomd.util.print_status_line();
        
        # initialize the base class
        md_pair.pair.__init__(self, r_cut, nlist, name);
        
        # update the neighbor list
        if d_max is None :
            sysdef = hoomd.context.current.system_definition;
            d_max = sysdef.getParticleData().getMaxDiameter()
            hoomd.context.msg.notice(2, "Notice: polydisperse set d_max=" + str(d_max) + "\n");

        # SLJ requires diameter shifting to be on
        self.nlist.cpp_nlist.setDiameterShift(True);
        self.nlist.cpp_nlist.setMaximumDiameter(d_max);
        
        # create the c++ mirror class
        if (model == "polydisperse12"):
            if not hoomd.context.exec_conf.isCUDAEnabled():
                self.cpp_force = _polymd.PotentialPairPolydisperse(hoomd.context.current.system_definition, self.nlist.cpp_nlist, self.name);
                self.cpp_class = _polymd.PotentialPairPolydisperse;
            else:
                self.nlist.cpp_nlist.setStorageMode(_md.NeighborList.storageMode.full);
                self.cpp_force = _polymd.PotentialPairPolydisperseGPU(hoomd.context.current.system_definition, self.nlist.cpp_nlist, self.name);
                self.cpp_class = _polymd.PotentialPairPolydisperseGPU;
        elif (model == "lennardjones"):
            if not hoomd.context.exec_conf.isCUDAEnabled():
                self.cpp_force = _polymd.PotentialPairPolydisperseLJ(hoomd.context.current.system_definition, self.nlist.cpp_nlist, self.name);
                self.cpp_class = _polymd.PotentialPairPolydisperseLJ;
            else:
                self.nlist.cpp_nlist.setStorageMode(_md.NeighborList.storageMode.full);
                self.cpp_force = _polymd.PotentialPairPolydisperseLJGPU(hoomd.context.current.system_definition, self.nlist.cpp_nlist, self.name);
                self.cpp_class = _polymd.PotentialPairPolydisperseLJGPU;
        elif (model == "polydisperse18"):
            if not hoomd.context.exec_conf.isCUDAEnabled():
                self.cpp_force = _polymd.PotentialPairPolydisperse18(hoomd.context.current.system_definition, self.nlist.cpp_nlist, self.name);
                self.cpp_class = _polymd.PotentialPairPolydisperse18;
            else:
                self.nlist.cpp_nlist.setStorageMode(_md.NeighborList.storageMode.full);
                self.cpp_force = _polymd.PotentialPairPolydisperse18GPU(hoomd.context.current.system_definition, self.nlist.cpp_nlist, self.name);
                self.cpp_class = _polymd.PotentialPairPolydisperse18GPU;
        elif (model == "polydisperse10"):
            if not hoomd.context.exec_conf.isCUDAEnabled():
                self.cpp_force = _polymd.PotentialPairPolydisperse10(hoomd.context.current.system_definition, self.nlist.cpp_nlist, self.name);
                self.cpp_class = _polymd.PotentialPairPolydisperse10;
            else:
                self.nlist.cpp_nlist.setStorageMode(_md.NeighborList.storageMode.full);
                self.cpp_force = _polymd.PotentialPairPolydisperse10GPU(hoomd.context.current.system_definition, self.nlist.cpp_nlist, self.name);
                self.cpp_class = _polymd.PotentialPairPolydisperse10GPU;
        elif (model == "polydisperse106"):
            if not hoomd.context.exec_conf.isCUDAEnabled():
                self.cpp_force = _polymd.PotentialPairPolydisperseLJ106(hoomd.context.current.system_definition, self.nlist.cpp_nlist, self.name);
                self.cpp_class = _polymd.PotentialPairPolydisperseLJ106;
            else:
                self.nlist.cpp_nlist.setStorageMode(_md.NeighborList.storageMode.full);
                self.cpp_force = _polymd.PotentialPairPolydisperseLJ106GPU(hoomd.context.current.system_definition, self.nlist.cpp_nlist, self.name);
                self.cpp_class = _polymd.PotentialPairPolydisperseLJ106GPU;
        hoomd.context.current.system.addCompute(self.cpp_force, self.force_name);

        # setup the coefficient options
        if (model == "polydisperse12"):
            self.required_coeffs = ['v0', 'eps', 'scaledr_cut'];
            self.pair_coeff.set_default_coeff('v0', 1.0);
            self.pair_coeff.set_default_coeff('eps', 0.2);
            self.pair_coeff.set_default_coeff('scaledr_cut', 1.25);
        elif (model == "polydisperse18"):
            self.required_coeffs = ['v0', 'eps', 'scaledr_cut'];
            self.pair_coeff.set_default_coeff('v0', 1.0);
            self.pair_coeff.set_default_coeff('eps', 0.0);
            self.pair_coeff.set_default_coeff('scaledr_cut', 1.25);
        elif (model == "polydisperse10"):
            self.required_coeffs = ['v0', 'eps', 'scaledr_cut'];
            self.pair_coeff.set_default_coeff('v0', 1.0);
            self.pair_coeff.set_default_coeff('eps', 0.0416667);
            self.pair_coeff.set_default_coeff('scaledr_cut', 1.48);
        elif (model == "polydisperse106"):
            self.required_coeffs = ['v0', 'eps', 'scaledr_cut'];
            self.pair_coeff.set_default_coeff('v0', 1.0);
            self.pair_coeff.set_default_coeff('eps', 0.1);
            self.pair_coeff.set_default_coeff('scaledr_cut', 2.5);
        elif (model == "lennardjones"):
            self.required_coeffs = ['v0', 'eps', 'scaledr_cut'];
            self.pair_coeff.set_default_coeff('v0', 1.0);
            self.pair_coeff.set_default_coeff('eps', 0.2);
            self.pair_coeff.set_default_coeff('scaledr_cut', 2.5);
    def process_coeff(self, coeff):
        v0 = coeff['v0'];
        eps = coeff['eps'];
        scaledr_cut = coeff['scaledr_cut'];
        
        return _hoomd.make_scalar3(v0,eps,scaledr_cut);

class polydisperseyukawa(md_pair.pair):
    R""" Polydisperse's yukawa

    Args:
        r_cut (float): Default cutoff radius (in distance units).
        nlist (:py:mod:`hoomd.md.nlist`): Neighbor list
        name (str): Name of the force instance.

    :py:class:`lj` specifies that a Lennard-Jones pair potential should be applied between every
    non-excluded particle pair in the simulation.

    .. math::
        :nowrap:

        \begin{eqnarray*}
        V_{\mathrm{LJ}}(r)  = & 4 \varepsilon \left[ \left( \frac{\sigma}{r} \right)^{12} -
                          \alpha \left( \frac{\sigma}{r} \right)^{6} \right] & r < r_{\mathrm{cut}} \\
                            = & 0 & r \ge r_{\mathrm{cut}} \\
        \end{eqnarray*}

    See :py:class:`pair` for details on how forces are calculated and the available energy shifting and smoothing modes.
    Use :py:meth:`pair_coeff.set <coeff.set>` to set potential coefficients.

    The following coefficients must be set per unique pair of particle types:

    - :math:`\varepsilon` - *epsilon* (in energy units)
    - :math:`\sigma` - *sigma* (in distance units)
    - :math:`\alpha` - *alpha* (unitless) - *optional*: defaults to 1.0
    - :math:`r_{\mathrm{cut}}` - *r_cut* (in distance units)
      - *optional*: defaults to the global r_cut specified in the pair command
    - :math:`r_{\mathrm{on}}`- *r_on* (in distance units)
      - *optional*: defaults to the global r_cut specified in the pair command

    Example::

        nl = nlist.cell()
        lj = pair.lj(r_cut=3.0, nlist=nl)
        lj.pair_coeff.set('A', 'A', epsilon=1.0, sigma=1.0)
        lj.pair_coeff.set('A', 'B', epsilon=2.0, sigma=1.0, alpha=0.5, r_cut=3.0, r_on=2.0);
        lj.pair_coeff.set('B', 'B', epsilon=1.0, sigma=1.0, r_cut=2**(1.0/6.0), r_on=2.0);
        lj.pair_coeff.set(['A', 'B'], ['C', 'D'], epsilon=1.5, sigma=2.0)

    """
    def __init__(self, r_cut, nlist, model):
        hoomd.util.print_status_line();

        # tell the base class how we operate

        # initialize the base class
        md_pair.pair.__init__(self, r_cut, nlist, "polydisperse_yukawa");
        
        # create the c++ mirror class
        if not hoomd.context.exec_conf.isCUDAEnabled():
            self.cpp_force = _polymd.PotentialPairPolydisperseYukawa(hoomd.context.current.system_definition, self.nlist.cpp_nlist, self.name);
            self.cpp_class = _polymd.PotentialPairPolydisperseYukawa;
        else:
            self.nlist.cpp_nlist.setStorageMode(_md.NeighborList.storageMode.full);
            self.cpp_force = _polymd.PotentialPairPolydisperseYukawaGPU(hoomd.context.current.system_definition, self.nlist.cpp_nlist, self.name);
            self.cpp_class = _polymd.PotentialPairPolydisperseYukawaGPU;
        hoomd.context.current.system.addCompute(self.cpp_force, self.force_name);

        # setup the coefficient options
        self.required_coeffs = ['v0', 'eps', 'scaledr_cut', 'kappa'];
        self.pair_coeff.set_default_coeff('v0', 10.0);
        self.pair_coeff.set_default_coeff('eps', 0.0);
        self.pair_coeff.set_default_coeff('scaledr_cut', 3.0);
        self.pair_coeff.set_default_coeff('kappa', 3.0);
    
    def process_coeff(self, coeff):
        v0 = coeff['v0'];
        eps = coeff['eps'];
        scaledr_cut = coeff['scaledr_cut'];
        kappa = coeff['kappa'];
        
        return _hoomd.make_scalar4(v0,eps,scaledr_cut,kappa);

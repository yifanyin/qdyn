# QDYN <img src="docs/img/qdyn_logo_small.jpeg" alt="QDYN logo" align="right" />

## A Quasi-DYNamic earthquake simulator

[![Build Status](https://travis-ci.com/ydluo/qdyn.svg?branch=master)](https://travis-ci.com/ydluo/qdyn) [![Documentation Status](https://readthedocs.org/projects/ansicolortags/badge/?version=latest)](https://ydluo.github.io/qdyn/) [![GPLv3 license](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://perso.crans.org/besson/LICENSE.html)

*QDYN* is a boundary element software to simulate earthquake cycles (seismic and aseismic slip on tectonic faults) under the quasi-dynamic approximation (quasi-static elasticity combined with radiation damping) on faults governed by rate-and-state friction and embedded in elastic media.

*QDYN* includes various forms of rate-and-state friction and state evolution laws, and handles non-planar fault geometry in 3D and 2D media, as well as spring-block simulations. Loading is controlled by remote displacement, steady creep or oscillatory load. In 3D it handles free surface effects in a half-space, including normal stress coupling. The medium surrounding the fault is linear, isotropic and elastic, and may be uniform or (in 2D) contain a damaged layer.

*QDYN* implements adaptive time stepping, shared-memory parallelization, and can deal with multi-scale earthquake cycle simulations with fine details in both time and space. It is equipped with user-friendly MATLAB and Python interfaces and graphical output utilities.

To get started with QDYN, see [the documentation](https://ydluo.github.io/qdyn/).


------

## Features

- rate-and-state friction, with velocity cut-offs, aging and slip laws

- microphysically based frictional model (*Chen-Niemeijer-Spiers* model)

- heterogeneous frictional properties

- slow and fast, aseismic and seismic slip transients

- dynamic weakening (thermal pressurization)

- non-planar faults (currently limited to variable dip, rectangular elements)

- 3D, 2D and 1D (spring-block)

- tectonic and transient loads

- normal stress coupling

- faults surrounded by damaged zones

- MATLAB and Python wrappers, and graphic output display utilities

- parallelized for shared memory systems (OpenMP)

- parallelized for distributed memory systems (MPI)


## LSODA solver

This fork contains the implementation of LSODA solver from the [ODEPACK](https://computing.llnl.gov/projects/odepack), Lawrence Livermore National Laboratory.

The solver allows faster time-stepping when integrating the rate-and-state friction ODEs in larger relative tolerance and not loosing the overall model behaviors. The LSODA solver only works for solving rate-and-state friction law.

 The CNS friction law and other evolutions like thermal pressurization is not tested in this fork.

--------------------------------

## Downloads, documentation and support

Previous (stable) versions of QDYN can be downloaded from the [release page](https://github.com/ydluo/qdyn/releases). Development versions are available as separate [branches](https://github.com/ydluo/qdyn/branches) following the naming convention `release/x.x.x`.

To install QDYN please follow the *Getting started* section in [the documentation](http://ydluo.github.io/qdyn/).

Questions, feedback or suggestions can be submitted via our [issue tracking system](https://github.com/ydluo/qdyn/issues).



-------------------------

## Introduction Poster

![](https://lh4.googleusercontent.com/-OjKBE5_Ipf8/T9wk2GtVRXI/AAAAAAAAABg/a1diUWu7tFU/s763/Poster_QDYN.jpg)

[Download Poster](http://code.google.com/p/qdyn/downloads/detail?name=Poster_QDYN.pdf) 

-------------------------


## Featured Simulations

### [Simulation_Tohoku](https://github.com/ydluo/qdyn/wiki/Simulation_Tohoku)
![](https://lh5.googleusercontent.com/-JPaTpBXo5eA/USdSArzQ0QI/AAAAAAAAKew/9wnVu30Lhf4/s900/Tohoku_cycle_logo.gif)

### [Simulation_Cascadia_Tremor](https://github.com/ydluo/qdyn/wiki/Simulation_Cascadia_Tremor)
![](https://lh5.googleusercontent.com/-a_2MRxcUgf8/T-v2JCjmxBI/AAAAAAAAAB8/NlQTwfra4fY/s900/Tremor_3D_Cascadia.gif)

------------------------
## Developers

*(listed alphabetically)*

[Jean-Paul Ampuero](http://www.seismolab.caltech.edu/ampuero_jp.html) (IRD/UCA, Géoazur, France; Caltech Seismolab, USA)

[Martijn van den Ende](https://www.linkedin.com/in/martijnvandenende) (Université Côte d'Azur, Géoazur, France)

[Percy Galvez](https://smi.kaust.edu.sa/Pages/People-Galvez.aspx) (KAUST, Saudi Arabia; AECOM, Switzerland)

[Benjamin Idini](http://www.seismolab.caltech.edu/idini_b.html) (Caltech Seismolab, USA)

[Yingdi Luo](https://science.jpl.nasa.gov/people/YLuo/) (NASA JPL, USA)

-------------------------

## Suggested References

#### For all uses of the QDYN software

Luo, Y., Ampuero, J. P., Galvez,  P., van den Ende, M., & Idini, B. (2017). 
QDYN: a Quasi-DYNamic earthquake simulator (v1.1) [Data set]. Zenodo. doi:10.5281/zenodo.322459  
 [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.322459.svg)](https://doi.org/10.5281/zenodo.322459)

#### For simulations on heterogeneous faults

Luo, Y., & Ampuero, J. P. (2018). 
Stability of faults with heterogeneous friction properties and effective normal stress. 
Tectonophysics, 733, 257-272, doi:[10.1016/j.tecto.2017.11.006](https://doi.org/10.1016/j.tecto.2017.11.006)

Luo, Y., Ampuero, J. P., Miyakoshi, K., & Irikura, K. (2017). Surface rupture effects on earthquake moment-area scaling relations. Pure and Applied Geophysics, 174(9), 3331-3342, doi:[10.1007/s00024-017-1467-4](https://link.springer.com/article/10.1007/s00024-017-1467-4).
In Topical Volume on *"Best Practices in Physics-based Fault Rupture Models for Seismic Hazard Assessment of Nuclear Installations"*.  
[PDF](https://rdcu.be/oOL9)

Luo, Y., & Ampuero, J. P. (2012), Simulation of Complex Tremor Migration Patterns, AGU Fall Meeting 2012 Abstract S44B-02

Luo, Y., & Ampuero, J. P. (2011), Numerical Simulation of Tremor Migration Triggered by Slow Slip and Rapid Tremor Reversals, AGU Fall Meeting 2011, abstract S33C-02

#### For the microphysically based (CNS) simulations

van den Ende, M. P. A., Chen, J., Ampuero, J. P., & Niemeijer, A. R. (2018).
A comparison between rate-and-state friction and microphysical models, based on numerical simulations of fault slip.
Tectonophysics, 733, 273-295, doi:[10.1016/j.tecto.2017.11.040](https://doi.org/10.1016/j.tecto.2017.11.040)

#### For simulations on faults surrounded by damaged zones

Idini, B., & Ampuero, J. P. (2017).
Rupture complexity promoted by damaged fault zones in earthquake cycle models.
AGU Fall Meeting 2017, abstract T41C-0632.  
Poster [PDF](https://www.essoar.org/doi/abs/10.1002/essoar.10500080.1), doi:[10.1002/essoar.10500080.1](https://dx.doi.org/10.1002/essoar.10500080.1)

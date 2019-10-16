# User’s Manual for Walkabout Version 1.0

## 1. Introduction

Walkabout Version 1.0 (LA-CC-11-033) performs random walk particle tracking
simulations of solute transport based on groundwater flow solutions that use
fully unstructured control volume grids. Walkabout 1.0 is designed to work
within the FEHM (Zyvoloski, 2007) code system, and accepts groundwater flow
solutions from FEHM and computational mesh descriptions from LaGrit (Los Alamos
Grid Toolbox, http://lagrit.lanl.gov, 2011). Walkabout also provides input to
the PLUMECALC (Robinson, et al. 2011) particle tracking post processor. Although
designed to work within the FEHM system, other control-volume flow codes and
mesh generators could, with appropriate reformatting of the output, provide the
flow fields and mesh descriptions.

A typical workflow for Walkabout within the FEHM system would use LaGrit to
generate unstructured grids. FEHM then provides a discretized representation of
the steady-state flow field to Walkabout. Given this discrete solution,
Walkabout then reconstructs a groundwater flow field, and performs the random
walk particle tracking calculation. Output is provided in a form compatible with
the PLUMECALC software. PLUMECALC may be used to efficiently post-process the
particle tracking results from Walkabout to add effects of retention/retardation
and arbitrary source histories. An option exists to also record particle
positions versus time, thus allowing other post-processing codes or
visualization systems to be used.

Other features and limitations of Walkabout are

1.  Walkabout works on fully unstructured tetrahedral meshes in three spatial
dimensions. Two-dimensional meshes and meshes other than tetrahedral meshes are
not supported in Version 1.0.

2.  A control-volume solution for steady groundwater flow is required. Finite-
element solutions are not supported. Transient flow is not supported in Version
1.0.

3.  All particles are launched at time zero. The PLUMECALC software may be used
to postprocess the resulting particle tracks to obtain concentration for an
arbitrary source history.

4.  Particles are moved through the system without decay or retardation. The
PLUMECALC system may be used to postprocess the particle tracks to obtain
concentration with decay and matrix diffusion or other retardation/retention
processes.

5.  The Burnett and Frind (1987) model for dispersion coefficient is presumed.

6.  Full heterogeneity in porosity, liquid density, liquid saturation index, and
dispersivity is supported.
 
This document summarizes the technical basis for the code, describes the input
and output formats, documents verification testing, and provides example test
cases that exercise the full range of capabilities.
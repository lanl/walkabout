# Walkabout Version 1 (LA-CC-11-033) (Open Source) 


Performs random walk particle tracking simulations of solute transport based on groundwater flow solutions that use fully unstructured control volume grids. Walkabout 1.0 is designed to work within the FEHM (Zyvoloski, 2007) code system, and accepts groundwater flow solutions from FEHM and computational mesh descriptions from LaGrit (Los Alamos Grid Toolbox).

A typical workflow for Walkabout within the FEHM system would use LaGrit to generate unstructured grids. FEHM then provides a discretized representation of the steady-state flow field to Walkabout. Given this discrete solution, Walkabout then reconstructs a groundwater flow field, and performs the random walk particle tracking calculation. Output is provided in a form compatible with the PLUMECALC software. PLUMECALC may be used to efficiently post-process the particle tracking results from Walkabout to add effects of retention/retardation and arbitrary source histories. An option exists to also record particle positions versus time, thus allowing other post-processing codes or visualization systems to be used.

[External Contributer Agreement (NON LANL Employees)](https://www.clahub.com/agreements/lanl/walkabout)
   
[Copyright text for Walkabout](./COPYRIGHT.md)
    
[Install directions for Walkabout](./INSTALL.md)

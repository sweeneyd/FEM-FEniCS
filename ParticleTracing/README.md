# FEM-FEniCS/ParticleTracing

## Project Description
This mini-project contains scripts to simulate fluid flow through a pipe using the FEniCS finite element framework in Docker. The pipe was made in FreeCAD (`cylinder_flow.step`), with GMSH used for first-order tetrahedral meshing (`cylinder_flow.geo`). FEniCS code is written in Python and run through a temporary Docker instance (`fenics.sh`) and is used to first call convert the mesh file output from GMSH to a .xml file, then run the fluid model script (`pipe_flow.py`), and finally the particle tracing script (`particle_flow.py`).

## File Summary
- `cylinder_flow.step`: STEP file created in FreeCAD of a 1 mm radius x 10 mm height cylinder
- `cylinder_flow.geo`: GMSH file containing meshing configuration
- `fenics.sh`: Bash script to facilitate running a temporary FEniCS instance in Docker
- `flow_simulation.sh`: Bash script facilitating meshing, mesh conversion, flow simulation, and particle tracing
- `pipe_flow.py`: 3D FEM simulation of pressure-drivein incompressible pipe flow using Taylor-Hood elements
- `particle_flow.py`: Particle tracing algorithm

## Output Files
- `cylinder_flow.msh`: Mesh file created by GMSH
- `cylinder_flow.xml`: Mesh file converted to .xml from .msh using dolfin-convert
- `cylinder_flow_physical_region.xml`: Mesh file with subdomain markers
- `cylinder_flow_facet_region.xml`: Mesh file with boundary markers
- `velocity.h5`: File containing velocity data created during FEM simulation
- `pressure.h5`: File containing pressure data created during FEM simulation
- `particles.csv`: File containing particle identifier, position (x,y,z), and normal velocity vector components

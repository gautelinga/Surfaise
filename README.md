# Surfaise
A Python package for solving PDEs on parametrized surfaces using FEniCS.

To install:
```
sudo python3 setup.py install
```
or
```
python3 setup.py install --user
```

Demonstrations, tests, and examples to appear...

### Features 
* Differential geometry:
For a given mapping from reference element to a curved 2D surface embedded in 3D space, calculates the metric tensor, curvature tensor and invariant quantities symbolically using SymPy. Supported shapes: Ellipsoid, Sphere, Cylinder, Gaussian bump, random bumpy surface, Torus, etc. Feel free to add shapes as a pull request. (Note that ellipsoid and sphere require extra care at the poles.)
* Timeseries (neat exporting):
```
from surfaise.common.io import Timeseries
```
* Postprocessing (preliminary):
```
surfaise-postprocess folder=path/to/folder-with-Timeseries-folder-inside/0/ method=plot
```
* Visualization (maps the solution on the reference element to the surface):
```
surfaise-visualize path/to/folder
```
This creates a file `visualize.xdmf` that you can open in e.g. ParaView.

### Dependencies
* FEniCS (Dolfin, UFL, ...)
* SymPy
* NumPy
* json
* cloudpickle

### Developed and maintained by:
* Bjarke Frost Nielsen
* Gaute Linga

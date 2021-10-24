# openFDTD 

An open source program to simulate electromagnetic
signals and interactions with materials. 

The simulations employ the FDTD method, on a staggered
grid as formulated by Kane S. Yee.

## Installation

In order to run this program, clone the repository and
from the top directory run:

```bash
cmake -S . -B build
make -C build
```

The binary file `openFDTD` is generated in the project's
`bin` directory.

## Gallery
|<img src="./gallery/1D_FDTD_demo.gif" width="480" height="400"/>|
|:--:|
| *1D FDTD with TFSF, material, and 1st order Mur ABC* |

## Planned Features
* One dimension
  - [x] Simple propagation
  - [x] Simple Mur boundary conditions
  - [x] Total Field - Scattered Field formulation
  - [x] Materials in grid
  - [x] Automatic calculation of required parameters
  - [x] Specral analysis
  - [ ] Lossy materials
  - [ ] More sophisticated boundary conditions
* Two dimension
* Three dimension
  - [ ] Implement a mesher for 3D objects (ex: `.stl` file to cubic grid)
* General
  - [ ] Create a GUI

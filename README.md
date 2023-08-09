This repository archives the model source code, input data, and experiment setups needed to run the forward simulations in Huth et al (2023). "Simulating the processes controlling ice-shelf rift paths using damage mechanics". Journal of Glaciology.

An installation of Elmer FEM with Elmer/Ice is required. The simulations in the paper were run using Elmer v9.0.

The source code for the SSA and damage model is given in `Model/PROG`. It is based on the MPM-SSA-Damage model: Huth et al (2021). "A Generalized Interpolation Material Point Method for Shallow Ice Shelves. 2: Anisotropic Nonlocal Damage Mechanics and Rift Propagation". JAMES. https://doi.org/10.1029/2020MS002292. This model is modified here to run using the Lagrangian finite element method, and to include the rift-flank boundary scheme, as detailed in Huth et al., 2023.

The mesh for the simulations is contained within ./Larsen_C_Grid

The input data (e.g. ice geometry, bed elevation, ice temperature, and the inverted fields from Figure 5) for all forward simulations is provided as an Elmer restart file: `./Larsen_C_Grid/m_init.result`.

***
* Compile the model following the README in Model/PROG
* uncompress the mesh file: `(cd Larsen_C_Grid/ ; tar -xzvf m_Init.tar.gz ; cd ..)`
* Run the simulations (summarized in Figure 8) using the scripts provided in `./Model/Simulations`, which contains a README with instructions.
***

Please contact me with any questions!

Alex Huth
Alexander.Huth@noaa.gov

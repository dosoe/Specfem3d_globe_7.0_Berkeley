# SPECFEM 3D Globe with Berkeley Subroutines
The software package SPECFEM3D_GLOBE simulates three-dimensional global and regional seismic wave propagation and performs full waveform imaging (FWI) or adjoint tomography based upon the spectral-element method (SEM). Effects due to lateral variations in compressional-wave speed, shear-wave speed, density, a 3D crustal model, ellipticity, topography and bathymetry, the oceans, rotation, and self-gravitation are included. The package can accommodate full 21-parameter anisotropy [Chen and Tromp, 2007] as well as lateral variations in attenuation [Savage et al., 2010].

#### Compilation
We have tested the compilation of the SPECFEM3D_globe on the two supercomputing clusters - Stampede2, and NERSC CORI.

##### To compile on NERSC CORI (intel compiler):
```
make -f Makefile.intel.CORI clean
make -f Makefile.intel.CORI
```


##### To compile on Stampede2 (intel compiler):
```
make -f Makefile.intel.XSEDE clean
make -f Makefile.intel.XSEDE
```

##### To compile on Perlmutter (gnu compiler)
1. Use gnu compiler: PrgEnv-gnu/8.3.3
    <!-- ./configure FC=ftn CC=cc MPIFC=ftn -->
    <!-- module swap PrgEnv-gnu/8.3.3 PrgEnv-cray/8.3.3 -->
1. Some codes that need fixes:
    1. src/specfem3D/mirror.F90 # line 1069
    1. src/specfem3D/setup_sources_receivers.f90
    1. src/specfem3D/iterate_time.F90
    1. src/specfem3D/get_cmt.f90
    1. src/specfem3D/ucb_dfour.f90 #implicit types are not defined

1. Compilation
    ```
    make -f Makefile.gnu.PERLM clean
    make -f Makefile.gnu.PERLM
    ```

#### Instructions to run

- Set the run params in `DATA/Par_file` and `setup/constants.h`

In `DATA/Par_file`:
- For PREM model `MODEL = PREM_A3d_1Dcrust`
- For SEMUCB model `MODEL = SEMUCB_A3d`
- Set filter passband for SEM source-time function:
    ```
    # Source frequency content (i.e., heaviside function)
    SOURCE_T1                       = 400.d0
    SOURCE_T2                       = 250.d0
    SOURCE_T3                       =  53.d0
    SOURCE_T4                       =  40.d0
    ```
- See [manual](https://specfem3d-globe.readthedocs.io/en/latest/03_running_the_mesher/) for the description of other parameters in the Par_file.

#### Instructions of hybrid simulation to output the mirror_src.dat 
- [x] generate a mirror.txt file where the displacement will be calculated by the specfem3d, and copy it to DATA/.
- [x] set the save_mirror in `setup/constants.h` to true


In `setup/constants.h`:
- Set `R_EARTH = 6368000.d0` for PREM with oceans
- Set `A3d_folder = "DATA/PREM_aniso_noocean_A3d.zeros/"` for PREM_A3d_1Dcrust
- Set `A3d_folder = "DATA/SEMUCB_A3d/"` for SEMUCB_A3d

- Main params:
    ```
    logical, parameter :: SAVE_MIRRORS = .false. 
    A3d_folder = "DATA/PREM_aniso_noocean_A3d.zeros/"
    double precision, parameter :: R_EARTH = 6371000.d0 
    ```

- Setting up the number of nodes for execution
    ```
    # compute total number of nodes needed
    NPROC_XI=`grep ^NPROC_XI DATA/Par_file | cut -c 34- `
    NPROC_ETA=`grep ^NPROC_ETA DATA/Par_file | cut -c 34- `
    NCHUNKS=`grep ^NCHUNKS DATA/Par_file | cut -c 34- `

    # total number of nodes is the product of the values read
    numnodes=$(( $NCHUNKS * $NPROC_XI * $NPROC_ETA ))
    
##### Checklist
- [x] Edit `DATA/Par_file`
- [x] Confirm `setup/constants.h`
- [x] Create `DATABASES_MPI` in the specfem home dir
- [x] Compile the script: `make clean; make` (system dependent Makefile)
- [x] Input the number of required processes in the slurm script. This will depend on the input parameters in the Par_file.

## Instructions to add a new model with format similar to the SEMUCB_WM1 
We first need to create a model directory (similar to `SEMUCB_A3d`) and add the information of this model into these files:
1. DATA/Par_file
1. setup/constants.h
1. src/shared/get_model_parameters.f90

#### References
1. [SPECFEM3D_Globe](https://github.com/SPECFEM/specfem3d_globe): original codes
1. [User Manual for SPECFEM3D_Globe](https://specfem3d-globe.readthedocs.io/en/latest/)
1. S. W. French, B. A. Romanowicz, Whole-mantle radially anisotropic shear velocity structure from spectral-element waveform tomography, Geophysical Journal International, Volume 199, Issue 3, December 2014, Pages 1303–1327, https://doi.org/10.1093/gji/ggu334

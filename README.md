
# Shift Current _(shift-current)_

> Calculates the shift current of a periodic system from electronic wavefunctions.

Original Authors: Steve M. Young, Fan Zheng

## Install

Edit the `make.inc` file with the appropriate paths and executables. Then run `make`.

See compile script and modify for your system. Use of intel compiler is strongly recommended.
The `iotk` module from Quantum Espresso must be presesent.

## Usage

First, generate wavefunctions from quantum espresso or wannier90. 
Next, run shiftcurrent code. See examples/MgO for an MgO quantum espresso example. See examples/ascii_input_example for reading in wavefunctions in ascii format, for a toy model.

If using quantum espresso, use `outdir='./'` for the quantum espresso output directory, and use `wf_collect= .true.` Use

```
KPOINTS automatic
n1 n2 n3 s1 s2 s3
```

Quantum Espresso must be compiled with the `__OLD_XML` option (present for versions <= 6.2.1)

Shift current calculation may be run on any number of processors. Parallelization is done over kpoints.   

### Input file description

Must be named "input_sc". This file takes the structure:

```
&input  
    dirname = ...  
    cutoff = ...  
    ...  
/
```

Detailed explanation of the fields are provided below.

#### dirname

type: char  
default: none  

QE prefix

#### cutoff

type: real  
default: 3.0  

Photon energy cutoff, Hartree

#### block_tol

type: real  
default: 0.001 

Eigenstates with energies within block_tol are placed in a single block. block_tol has units of Hartree.


#### broadwidth

type: real  
default: 0.003

Broadening of spectrum, Hartree.

#### sg_calc

type: logical  
default: .false.

Set to true to do second harmonic generation calculation.     


#### data_format

type: char  
default: "auto"

"auto" : automatically detects Quantum Espresso or Wannier90 data  
"ascii": reads in user supplied ascii data files, which must conform to the following standards:  

Full kpoint grid nk1 x nk2 x nk3 must be given, symmetries not implemented  
gvecs.dat: g-vectors in crystal coords. Each row is gx gy gz, which are integers  
kcrys.dat: k-points in crystal coords. Each row is kx ky kz, which are floats  

data-file: Fortran namelist style input file with these fields  
&systemdata  
a1 = a1x, a1y, a1z ! lattice vectors, bohr   
a2 = a2x, a2y, a2z  
a3 = a3x, a3y, a3z  
nelec = int  ! must be even, 2 electrons per band  
nbnd = int   
kptdim = nk1, nk2, nk3 !int  
kptoff = sk1, sk2, sk3 ! 1 means offset by half step, else 0  
nspin = 1 ! only nspin = 1 implemented for ascii input for now  
/

work/k$i/ where $i is an integer, contains 

work/k$i/ens, band energies in eV (one per row)  
work/k$i/evecs, wavefunctions in the format  
real(wfn1(g1)) imag(wfn1(g1)) real(wfn2(g1)) imag(wfn2(g1)) ...  
real(wfn1(g2)) imag(wfn1(g2)) real(wfn2(g2)) imag(wfn2(g2)) ...  
...


#### kpts_type

type: char    
default: "mp"  

"mp" : Monkhorst-pack grid.  

"nonmp" :Non Monkhorst-pack grid. Currently implemented only for nxnxn grid. Must be in the order of k1 fastest, k3 slowest.

"path": k-path, will just print dipole matrix elements in the dump files along the path


### Format of output files

`sc.ijk.dat` gives the shift current response tensor as a function of light frequency, defined by

J_k (omega) = sum_ij sigma_ijk (omega) * E_i * E_j

Alternatively, the current density can be written as a function of light intensity instead of electric field

J = (sigma/epsilon) I

where epsilon is the dielectric function and I = epsilon E^2 is the light intensity.

column 1 : light frequency [eV]  
column 2 : sigma response tensor [A V^-2]  
column 3 : sigma/epsilon [A W^-1]  

In spinful systems, the format of `sc.ijk.dat`

column 1 : light frequency [eV]  
column 2 : sigma response tensor (spin down) [A V^-2]  
column 3 : sigma response tensor (spin up) [A V^-2]  
column 4 : sigma response tensor (total) [A V^-2]  
column 5 : sigma/epsilon [A W^-1]  

When running second harmonic generation calculations, `sg.ijk.dat` gives the nonlinear susceptibility tensor as a function of light frequency. It prints

column 1 : light frequency [eV]  
column 2 : Re[sigma/(2 i omega)] [C V^-2]  
blank line   
column 1 : light frequency [eV]  
column 2 : Im[sigma/(2 i omega)] [C V^-2]  

The code prints detailed information in the dump files. The format of these files is

kpoints in fractional coordinates  
kpoints in 2pi/a coordinates    
vband1 vband2 cband1 cband2 vblk cblk  
occupation difference, transition energy  
shift vector integrated over constant energy surface (3 components)  
transition intensity, <c| p_i |v><v| p_j |c> (9 components, column order)  
shift current (column order)  

Each transition is calculated between a block of valence bands and a block of conduction bands.  
The index of the valence block is vblk, and it comprises bands [vband1,vband2].  
The index of the conduction block is cblk, and it comprises bands [cband1,cband2].  
The transition energy should have the hartree unit. Others should have SI unit. Check sc_tools.F90 where the dump files are output. See examples of reading dump files in shift_current_postprocessing/.

## Known bugs/issues
Sometimes, when there are many almost-degenerate bands, setting block_tol too high can result in apparently divergent shift current response (1e5 A/V^2). The error comes when calculating the shift vector matrix for bands in that block. The solution is to lower block_tol
# Momentum_matrix_from_steve
# Momentum_matrix_from_steve
# Momentum_matrix_from_steve

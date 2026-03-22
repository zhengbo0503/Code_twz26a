# Code_twz26

## 1. Description
Version: January 29, 2026  
Description: MATLAB codes for paper "*Computing accurate singular values using a mixed-precision one-sided Jacobi algorithm*" by Zhengbo Zhou, Françoise Tisseur, and Marcus Webb.   
This repository has been tested in the following configuration:
```
OS: macOS 26.0.1 25A362 arm64 
CPU: Apple M3 Pro 
GPU: Apple M3 Pro 
Memory: 36864MiB 
MATLAB version: 25.2.0.3042426 (R2025b) Update 1
```

## 2. Requirements 
To run the code, you need the following setup:

### (1) [Anymatrix](https://github.com/north-numerical-computing/anymatrix)
We are using version 1.4. 

### (2) [Advanpix Multiprecision Computing Toolbox](https://www.advanpix.com/)
You inevitably need some packages to simulate the quadruple precision. 
For our project, we chose to use **Advanpix (Version 1.4)**. 
If you use other packages, such as `vpa()` in the MATLAB Symbolic Math Toolbox,
you need to change all occurance of `mp()` in this repo to the desired command.

### (3) MEX facilities for LAPACK SVD routines
In order to compile files 
 - `sgesvj_mex.mexmaca64`,
 - `dgesvj_mex.mexmaca64`,
 - `sgejsv_mex.mexmaca64`, and 
 - `dgejsv_mex.mexmaca64`,

you need to run the the following build files in the corresponding folder.
#### (i) Example compilation of `sgesvj`
> **Warning**: This example is only for macOS system. and for other systems, you may need to change the compilation commands in `build_sgesvj_mex.m`.

Before compiling, you need the following two installations:
- [OpenBLAS](https://github.com/OpenMathLib/OpenBLAS):
We are using `OpenBLAS 0.3.29.dev`.
- [gfortran](https://fortran-lang.org/learn/os_setup/install_gfortran/): 
We are using `GNU Fortran 15.1.0 (Homebrew GCC 15.1.0)`.

Then, open the folder `get_sgesvj`, there are three files within the folder
- `sgesvj_mex.c`: the main c file for compiling the MEX file.
- `mex.h`: dependency for `sgesvj_mex.c`.
- `build_sgesvj_mex.m`: the setup file for `sgesvj_mex`.

Before run `build_sgesvj_mex.m` in MATLAB command window, you should specify
the settings in the file:
```matlab
openblas_lib = '/opt/homebrew/opt/openblas/lib';
gfortran_lib = '/opt/homebrew/Cellar/gcc/15.1.0/lib/gcc/15'
```
Once you update this part, you can run `build_sgesvj_mex.m` in MATLAB command window.
You will get `sgesvj_mex.mexmaca64` in the folder `get_sgesvj`. Move it to the main folder
`Code_twz26`, then you are good to go. 

## 3. Codes description
- The main algorithm is
   [`mposj.m`](https://github.com/zhengbo0503/Code_twz26/blob/main/mposj.m)
   which implements our algorithm with (single, double, quadruple) setting. 
- Another implementation of our algorithm using (single, single, double)
   setting is
   [`mposj_ssd`](https://github.com/zhengbo0503/Code_twz26/blob/main/mposj_ssd.m). 
- `test1.m`: (*Figure 1*) accuracy tests for random matrices with varying condition numbers. 
- `test2.m`: (*Figure 2*) accuracy tests for random matrices with varying number of columns.
- `test3.m`: (*Figure 3*) accuracy tests for matrices from Anymatrix. 
- `test4.m`: (*Figure 4*) timing tests for random matrices with varying number of columns.
- `test5.m`: (*Figure 5*) timing tests for random matrices using `mposj_ssd`.

## 4. Quick Example 
To run a quick example of our algorithm, you can run the following code in MATLAB command window:
```matlab
clc; 
m = 100; n = 80; 
A = randn(m,n); 
[U,S,V] = mposj(A); 
fprintf('Backward error: %e\n', norm(A - U*S*V', 'fro') / norm(A, 'fro'));
fprintf('Orthogonality of U: %e\n', norm(U'*U - eye(n), 'fro'));
fprintf('Orthogonality of V: %e\n', norm(V'*V - eye(n), 'fro'));
```

## 5. Warning
(Updated March 22, 2026)

Today, when I tried to use my code again, MATLAB started crashing
regardless of how I modified `dgejsv_mex.c`, see [crash file](matlab_crash_report.txt).
After several hours of debugging, I identified the root cause as a conflict between two `libomp` runtimes being loaded into the same process.

Starting from version `0.3.31`, installing [OpenBLAS from Homebrew](https://formulae.brew.sh/formula/openblas) introduces `libomp` as a dependency. This was not the case for version `0.3.30`; see [here](https://web.archive.org/web/20251130191823/https://formulae.brew.sh/formula/openblas).
However, MATLAB already ships with its own `libomp` for the MEX runtime. On my laptop, for example, it uses
```text
/Applications/MATLAB_R2025b.app/bin/maca64/libomp.dylib
```
As a result, when the MEX file is executed, two different `libomp` runtimes are loaded into the same MATLAB process, which leads to segmentation faults.

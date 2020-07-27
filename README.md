PAAKAT - An HPC in-situ analysis tool
=====================================


Motivation
----------

One of the main goals of this tool is related to its capability of handling data arising from the C++ VTK API.



# Table of contents
A. [Compilation](#Compilation)

1. [PAAKAT](#subparagraph1)   
2.  NEK5K V1901 

[comment]: #  (  )


------------
## A. Compilation <a name="Compilation"></a>


### Download 

```
git clone https://github.com/jmake/INSITU_ANALYSIS

```

This project contains tree directories 
```
PAAKAT, 
NEK5K1901, and 		
TARs. 
```


Compilation is divided in two main steps : 

A.1. **PAAKAT**

[PAAKAT][PAAKAT] is a fast and scalable ... 



A.2. **Nek5000**

[Nek5000][ALYA] is a fast and scalable ... 


[ALYA]: https://nek5000.mcs.anl.gov/ 
[PAAKAT]: https://PAAKAT.PAAKAT 

## A.1. PAAKAT 

PAAKAT is a colecction of tools which allows working with the vtk library. 

It includes some useful examples 
(*INSITU_ANALYSIS/PAAKAT/EXAMPLES*)
related to how to perform parallel operations 
between vtk-filters. 
Some of which have been successfully
tested
in thousands of processors. 
Examples related to code instrumentation are also included (*INSITU_ANALYSIS/PAAKAT/PATCHES*). 

Additionally, 
PAAKAT allows 
the execution and interaction of 
multiple-instances of differents (or same) programs 
by using the Multiple Program Multiple Data (MPMD) MPI syntaxis
(*INSITU_ANALYSIS/PAAKAT/LIBPLEPP*). 


### Directory structure

```
INSITU_ANALYSIS/PAAKAT
├── EXAMPLES
│   └── PY
│       └── MPMD02_1
│       └── SIMPLEST
│   └── CXX 
│       └── SIMPLEST
├── LIBPLEPP
│   ├── CMakeLists.txt
│   ├── Tests
│   ├── include
│   ├── src
│   └── wrappers
├── PATCHES
│   └── NEK5K
│       └── V1901
│           └── PY
└── SRCs
    ├── E01_SIMPLEST
    ├── E02_SIMPLEST
    ├── getIncludeAndLibs.py
    ├── insitu_py.hpp
    └── vtktools.hpp

```

### \$INSITU\_ANALYSIS/PAAKAT/EXAMPLES/CXX/SIMPLEST  

This is the simplest in-situ test. 

```
## SOURCES
ln -s ../../../SRCs/E01_SIMPLEST/* .
ln -s ../../../SRCs/vtktools.hpp

## COMPILATION 
FC=mpif90
CC=mpicc
CXX=mpicxx
ROOT_CATALYST=WHERE/ARE/PARAVIEW/SOURCES?

make -f makefile.dummy compilelib CXX=$CXX CC=$CC ROOT_CATALYST=$ROOT_CATALYST

```

**simplest.mac** 
contains a concrete example of how this example could be compiled in macOS Mojave 


Compilation achieved in Beskow system (**simplest.beskow** )

```
>> cmake . -DCMAKE_CXX_COMPILER=$(CXX) -DCMAKE_C_COMPILER=$(CC) -DParaView_DIR=$(ROOT_CATALYST) 
... 
-- Cray Programming Environment 2.6.1 C
-- Check for working C compiler: /opt/cray/pe/craype/2.6.1/bin/cc
...  
-- Cray Programming Environment 2.6.1 CXX
-- Check for working CXX compiler: /opt/cray/pe/craype/2.6.1/bin/CC
... 
-- Configuring done
-- Generating done
-- Build files have been written to: $INSITU_ANALYSIS/PAAKAT/EXAMPLES/CXX/SIMPLEST

>> make 
Scanning dependencies of target test_cxx.x
[ 50%] Building CXX object CMakeFiles/test_cxx.x.dir/test.cxx.o
[100%] Linking CXX executable test_cxx.x
[100%] Built target test_cxx.x

```


### \$INSITU\_ANALYSIS/PAAKAT/EXAMPLES/PY/SIMPLEST  

The main difference with the previous example is related 
to the module **vtkPVPythonCatalyst**. 
This module is necessary when 
**vtkCPPythonScriptPipeline** is used 
which is the case in the most of the *catalyst* examples. 


```
## SOURCES
ln -s ../../../SRCs/E01_SIMPLEST/* .
ln -s ../../../SRCs/vtktools.hpp

## COMPILATION 
FC=mpif90
CC=mpicc
CXX=mpicxx
ROOT_CATALYST=WHERE/ARE/PARAVIEW/SOURCES?

make -f makefile.dummy compilelib CXX=$CXX CC=$CC ROOT_CATALYST=$ROOT_CATALYST

```

**simplest.mac** 
contains a concrete example of how this example could be compiled in macOS Mojave 




## 1.A. Nek5000 


### Download 

From 

```
wget https://github.com/Nek5000/Nek5000/releases/download/v19.0/Nek5000-19.0.tar.gz

```
Additionally, 
Nek5000-V19.0.1 (*Nek5000-19.0.tar.gz*) could be found in TARs. 

In any case, it is assumed that sources are found 
in *INSITU_ANALYSIS/NEK5K1901/SRC*.


### Directory structure


```
INSITU_ANALYSIS/NEK5K1901/SRC
├── 3rd_party
├── bin
├── core
├── examples
├── run
├── short_tests
└── tools
```

### INSITU\_ANALYSIS/NEK5K1901/EXAMPLES/OCYL01_3  

This example corresponds to the one found in 
*INSITU_ANALYSIS/NEK5K1901/SRC/examples/ocyl*. 
Which has been included in *OCYL00_1*. 


It runs multiple-instances of Nek5000 using 
the *Multiple Program Multiple Data* (MPMD) MPI syntaxis : 

```
mpirun -np n1 ./program1 : -np n1 ./program1 : ... : -np nN ./programN 

```

Each of these programs (*program1, program2, ..., programN*) 
reads a *pipe.py* file 
which corresponds to the *catalyst pipeline*. 
In this particular case a line (*\*_line%t.pvtp*) is saved each time step. 
Additinally, *\*_grid%t.pvtu* 
the complete grid is saved each five time steps. 

MPI communication between python pipelines 
must be performed by using *vtk*. 
This example includes some possible options 
(*vtkMultiProcessController, vtkProcessModule, etc*)
which could be used for this purpose. 


#### Directory structure 


```

NEK5K1901/
├── EXAMPLES
│   ├── OCYL00_1
│   └── OCYL01_3
└── SRC 

INSITU_ANALYSIS/PAAKAT
├── EXAMPLES
│   └── PY
│       └── MPMD02_1
└── LIBPLEPP 

```

#### PAAKAT Compilation 

This example depends on 
*$INSITU\_ANALYSIS/PAAKAT/EXAMPLES/PY/MPMD02_1*
which is compiled in a similar way as  
*$INSITU\_ANALYSIS/PAAKAT/EXAMPLES/PY/SIMPLEST*


```
PLEPP_PATH=$PAAKAT_PATH/LIBPLEPP
CASE_PATH=$PAAKAT_PATH/EXAMPLES/PY/MPMD02_1
ROOT_CATALYST=

make --no-print-directory -C ${CASE_PATH} -f makefile.dummy PLEPP_ROOT=${PLEPP_PATH} ROOT_CATALYST=${ROOT_CATALYST} compilelib
make --no-print-directory -C ${CASE_PATH} -f makefile.dummy getincludesandlibs
```

Output files 
*libPLEPP.a,libs.nek5k,includes.nek5k,* and *mod_plepp.mod*
must be copied here. 



#### NEK Patch  

In order to run this example, 
some modifications are needed, 
which can be applied by doing  

```
cd $INSITU\_ANALYSIS/NEK5K1901/SRC/
cp -r $INSITU\_ANALYSIS/PAAKAT/PATCHES/NEK5K/V1901/PY/core/* . 
```

#### NEK compilation   

```
NEK_PATH=$INSITU\_ANALYSIS/NEK5K1901/SRC
CASEA_NAME=ocyl

export CATALYST_INCS=$(<$./includes.nek5k); 
export CATALYST_LIBS=$(<$./libs.nek5k); 
export USR+="${PAAKAT_FFILE}.o"
export LDFLAGS="${CATALYST_LIBS}"

export CASE="$CASEA_NAME"
export FC="ftn -J./ "
export FFLAGS+="-DPLEPP -Wno-argument-mismatch"
export CC="cc"
export CFLAGS="-D_Float128=__float128"
export NEK_SOURCE_ROOT="${NEK_PATH}"

bash $NEK_PATH/bin/makenek ${CASEA_NAME}

```


#### NEK execution    


```
echo $CASEA_NAME  >   SESSION.NAME
echo $INSITU_ANALYSIS/NEK5K1901/EXAMPLES/OCYL01_3 >>  SESSION.NAME

ln -s pipe.py Y.py 
ln -s pipe.py X.py 
ln -s pipe.py Z.py 

>> mpirun \
  -np 5 ./nek5000 --namei Y : \
  -np 6 ./nek5000 --namei X : \
  -np 7 ./nek5000 --namei Z

```


```


/----------------------------------------------------------\\
|      _   __ ______ __ __  ______  ____   ____   ____     |
|     / | / // ____// //_/ / ____/ / __ \\/ __ \\/ __ \\   |
|    /  |/ // __/  / ,<   /___ \\ / / / // / / // / / /    |
|   / /|  // /___ / /| | ____/ / / /_/ // /_/ // /_/ /     |
|  /_/ |_//_____//_/ |_|/_____/  \\___/ \\___/ \\___/      |
|                                                          |
|----------------------------------------------------------|
|                                                          |
| COPYRIGHT (c) 2008-2019 UCHICAGO ARGONNE, LLC            |
| Version:  19.0                                           |
| Web:      https://nek5000.mcs.anl.gov                    |
|                                                          |
\\----------------------------------------------------------/

 Number of MPI ranks :           5
 ...
 Number of MPI ranks :           6
 ... 
 Number of MPI ranks :           7
 ... 
[CreateVtkMPIController] 5 Processes Initialized!!
[TestProcessController] size.rank:5.0
[Initialize] nRanks:5 pyfile:'Y.py' 
 ... 
[CreateVtkMPIController] 6 Processes Initialized!!
[TestProcessController] size.rank:6.0
[Initialize] nRanks:5 pyfile:'X.py' 
 ... 
[CreateVtkMPIController] 7 Processes Initialized!!
[TestProcessController] size.rank:7.0
[Initialize] nRanks:5 pyfile:'Z.py' 
 ... 
[Coprocess] timeStep:10, request:1 
end of time-step loop
...
[Finalize]
run successful: dying ...
... 
```

**ocyl01_3.beskow** 
contains a concrete example of how this example could be compiled and executed in Beskow system. 
 
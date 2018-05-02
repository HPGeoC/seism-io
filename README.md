# SEISM-IO

## 1. Introduction
The SEISM-IO library is designed to reduce the amount of optimization efforts and lower the barrier of parallel I/O implementation. Currently, SEISM-IO library only works with structured meshes, because unstructured meshes require much more complexed partitioning. Compared with other generalized I/O libraries, our SEISM-IO library has many specialized functions which aim at improving the programming efficiency of seismic applications, such as grid partition and buffering output. We also develop an easy-to-use application programming interface (API) for both C and Fortran language, which integrates different initialization, open, read, write and finalize processes in underlying MPI-IO, PHDF5, PnetCDF and ADIOS I/O libraries. By calling this light-weight interface instead of calling the SEISM-IO directly, we hide complex low-level operation in the SEISM-IO framework and high-level operations in the underlying libraries from user so they can focus on their scientific difficulties. Although the target application we focused on are mostly seismic applications, the SEISM-IO library can be used by many HPC applications based on structured meshes.

## 2. Installation

The SEISM-IO library requires multiple 3rd party software packages, We strongly recommend our user use the software packages listed as followings.

### 2.1 Pre-Install

1. szip-2.1.1
```shell
wget https://support.hdfgroup.org/ftp/lib-external/szip/2.1.1/src/szip-2.1.1.tar.gz
tar -zxvf szip-2.1.1.tar.gz
cd szip-2.1.1
./configure --prefix=/path/to/szip-2.1.1
make install
```

2. zlib-1.2.11
```shell
wget https://zlib.net/zlib-1.2.11.tar.gz
tar -zxvf zlib-1.2.11.tar.gz
cd zlib-1.2.11
./configure --prefix=/path/to/zlib-1.2.11
make install
```

3. phdf5-1.8.19
```shell
wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8/hdf5-1.8.19/src/hdf5-1.8.19.tar.gz
tar -zxvf hdf5-1.8.19.tar.gz
cd hdf5-1.8.19
CC=mpicc FC=mpif90 ./configure --enable-parallel --prefix=/path/to/phdf5-1.8.19 --enable-fortran
make install
```

4. pnetcdf-1.7.0
```shell
wget http://cucis.ece.northwestern.edu/projects/PnetCDF/Release/parallel-netcdf-1.7.0.tar.gz
tar -zxvf parallel-netcdf-1.7.0.tar.gz
cd parallel-netcdf-1.7.0
./configure --prefix=/path/to/pnetcdf-1.7.0
make
make install prefix=/path/to/pnetcdf-1.7.0
```

5. mxml-2.11
```shell
wget https://github.com/michaelrsweet/mxml/releases/download/v2.11/mxml-2.11.tar.gz
mkdir mxml-2.11
tar -zxvf mxml-2.11.tar.gz -C mxml-2.11
cd mxml-2.11
./configure --prefix=/path/to/mxml-2.11
make install
```

6. adios-1.9.0
```shell
wget http://users.nccs.gov/~pnorbert/adios-1.9.0.tar.gz
tar -zxvf adios-1.9.0.tar.gz
cd adios-1.9.0
CC=mpicc FC=mpif90 ./configure --prefix=/path/to/adios-1.9.0 --with-mxml=/path/to/mxml-2.11
make
make install

note: if "./src/core/adios_internals.c" complains about "isnan" or "isfinite", replace line 4970
isnan ( data [size] )    => isnan ( (float)data [size] )
isfinite ( data [size] ) => isfinite ( (float)data [size] )
```

After installed all software packages above, add all libraries' path into environment as followings,
```
export PATH=/path/to/phdf5-1.8.19/bin:$PATH
export LD_LIBRARY_PATH=/path/to/phdf5-1.8.19/lib:$LD_LIBRARY_PATH

export PATH=/path/to/pnetcdf-1.7.0/bin:$PATH
export LD_LIBRARY_PATH=/path/to/pnetcdf-1.7.0/lib:$LD_LIBRARY_PATH

export PATH=/path/to/adios-1.9.0/bin:$PATH
export LD_LIBRARY_PATH=/path/to/adios-1.9.0/lib:$LD_LIBRARY_PATH

export PATH=/path/to/mxml-2.11/bin:$PATH
export LD_LIBRARY_PATH=/path/to/mxml-2.11/lib:$LD_LIBRARY_PATH
```

### 2.2 Install
1. general Install
```shell
git clone https://github.com/HPGeoC/seism-io.git
cd seism-io
autoreconf -fi
CC=mpicc FC=mpif90 ./configure --prefix=/path/to/seismio \
--with-phdf5=/path/to/phdf5-1.8.19 --with-pnetcdf=/path/to/pnetcdf-1.7.0 \
--with-adios=/path/to/adios-1.9.0 --with-mxml=/path/to/mxml-2.11
make
make install
```
2. Install on Comet
```shell
git clone https://github.com/HPGeoC/seism-io.git
cd seism-io
autoreconf -fi
CC=mpicc FC=mpif90 ./configure --prefix=/path/to/seismio \
--with-phdf5=/opt/hdf5/intel/mvapich2_ib --with-pnetcdf=/opt/netcdf/4.3.2/intel/mvapich2_ib \
--with-adios=/path/to/adios --with-mxml=/path/to/local
make
make install
```

3. Install on BlueWaters
```shell
git clone https://github.com/HPGeoC/seism-io.git
cd seism-io
autoreconf -fi
CC=cc FC=ftn ./configure --prefix=/path/to/seismio \
--with-phdf5=/opt/cray/hdf5-parallel/1.8.16/gnu/4.9 --with-pnetcdf=/opt/cray/parallel-netcdf/1.7.0/gnu/4.9 \
--with-adios=/path/to/adios --with-mxml=/sw/xe/mxml/2.7/cnl4.2_gnu4.8.2
make
make install
```

4. Install on Titan
```shell
git clone https://github.com/HPGeoC/seism-io.git
cd seism-io
autoreconf -fi
CC=cc FC=ftn ./configure --prefix=/path/to/seismio \
--with-phdf5=/opt/cray/hdf5-parallel/1.8.16/gnu/4.9 \
--with-pnetcdf=/opt/cray/parallel-netcdf/1.7.0/gnu/4.9 \
--with-adios=/ccs/home/dmu/project/Tools --with-mxml=/sw/xk6/mxml/2.9/cle5.2_gcc4.8.2
make
make install
```

## 3. Run Example
```shell
cd examples/c
[modify libraries path in the Makefile]
make
./clean.sh
mpirun -np 8 ./testInterface

OR

cd examples/f
[modify libraries path in the Makefile]
make
./clean.sh
mpirun -np 8 ./testInterface
```

## 4. License
SEISM-IO is licensed under BSD-2

## 5. Contact
If you have any question, please contact yfcui@sdsc.edu or visit our homepage hpgeoc.ucsd.edu.

**Yifeng Cui**

2017.01

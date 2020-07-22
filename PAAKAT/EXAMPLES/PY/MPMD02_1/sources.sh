find . -type l -delete

ln -s  ../../../LIBPLEPP/Tests/DUMMY/dummy.cxx
ln -s  ../../../LIBPLEPP/Tests/DUMMY/dummy.f90
ln -s  ../../../SRCs/getIncludeAndLibs.py
ln -s  ../../../SRCs/insitu_py.hpp
ln -s  ../../../LIBPLEPP/Tests/DUMMY/makefile.dummy
ln -s  ../../../SRCs/vtktools.hpp
ln -s  ../../../SRCs/E01_SIMPLEST/wrapper.cxx

make -f makefile.dummy clean 


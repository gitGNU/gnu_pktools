//SWIG interface for pktools
%include <std_string.i>
%include <std_vector.i>
%include <std_iostream.i>
%include <std_shared_ptr.i>
%shared_ptr(ImgRaster)

%module pktools
%{
  #include <memory>
  #include "config.h"
  #include "imageclasses/ImgRaster.h"
  #include "imageclasses/ImgCollection.h"
  #include "apps/AppFactory.h"
  #include "algorithms/Filter2d.h"
  %}

%inline %{
extern "C"
{
 void *__dso_handle = 0;
}
 %}

//Parse the header file
%include "imageclasses/ImgRaster.h"
%include "imageclasses/ImgCollection.h"
%include "apps/AppFactory.h"
%include "algorithms/Filter2d.h"
// Instantiate some templates
%template(getGDALDataType_char) getGDALDataType<char>;
%template(getGDALDataType_uchar) getGDALDataType<unsigned char>;
%template(getGDALDataType_ushort) getGDALDataType<unsigned short>;
%template(getGDALDataType_short) getGDALDataType<short>;
%template(getGDALDataType_uint) getGDALDataType<unsigned int>;
%template(getGDALDataType_int) getGDALDataType<int>;
%template(getGDALDataType_ulong) getGDALDataType<unsigned long>;
%template(getGDALDataType_long) getGDALDataType<long>;
%template(getGDALDataType_float) getGDALDataType<float>;
%template(getGDALDataType_double) getGDALDataType<double>;

%template(readData_char) ImgRaster::readData<char>;
%template(readData_uchar) ImgRaster::readData<unsigned char>;
%template(readData_ushort) ImgRaster::readData<unsigned short>;
%template(readData_short) ImgRaster::readData<short>;
%template(readData_uint) ImgRaster::readData<unsigned int>;
%template(readData_int) ImgRaster::readData<int>;
%template(readData_ulong) ImgRaster::readData<unsigned long>;
%template(readData_long) ImgRaster::readData<long>;
%template(readData_float) ImgRaster::readData<float>;
%template(readData_double) ImgRaster::readData<double>;

%template(readDataBlock_char) ImgRaster::readDataBlock<char>;
%template(readDataBlock_uchar) ImgRaster::readDataBlock<unsigned char>;
%template(readDataBlock_ushort) ImgRaster::readDataBlock<unsigned short>;
%template(readDataBlock_short) ImgRaster::readDataBlock<short>;
%template(readDataBlock_uint) ImgRaster::readDataBlock<unsigned int>;
%template(readDataBlock_int) ImgRaster::readDataBlock<int>;
%template(readDataBlock_ulong) ImgRaster::readDataBlock<unsigned long>;
%template(readDataBlock_long) ImgRaster::readDataBlock<long>;
%template(readDataBlock_float) ImgRaster::readDataBlock<float>;
%template(readDataBlock_double) ImgRaster::readDataBlock<double>;

%template(writeData_char) ImgRaster::writeData<char>;
%template(writeData_uchar) ImgRaster::writeData<unsigned char>;
%template(writeData_ushort) ImgRaster::writeData<unsigned short>;
%template(writeData_short) ImgRaster::writeData<short>;
%template(writeData_uint) ImgRaster::writeData<unsigned int>;
%template(writeData_int) ImgRaster::writeData<int>;
%template(writeData_ulong) ImgRaster::writeData<unsigned long>;
%template(writeData_long) ImgRaster::writeData<long>;
%template(writeData_float) ImgRaster::writeData<float>;
%template(writeData_double) ImgRaster::writeData<double>;

%template(writeDataBlock_char) ImgRaster::writeDataBlock<char>;
%template(writeDataBlock_uchar) ImgRaster::writeDataBlock<unsigned char>;
%template(writeDataBlock_ushort) ImgRaster::writeDataBlock<unsigned short>;
%template(writeDataBlock_short) ImgRaster::writeDataBlock<short>;
%template(writeDataBlock_uint) ImgRaster::writeDataBlock<unsigned int>;
%template(writeDataBlock_int) ImgRaster::writeDataBlock<int>;
%template(writeDataBlock_ulong) ImgRaster::writeDataBlock<unsigned long>;
%template(writeDataBlock_long) ImgRaster::writeDataBlock<long>;
%template(writeDataBlock_float) ImgRaster::writeDataBlock<float>;
%template(writeDataBlock_double) ImgRaster::writeDataBlock<double>;

%ignore ImgRaster::operator=;
// Instantiate templates used by example
/* %rename(__assignImgRaster__) ImgRaster::operator=; */
/*   %template(ImgVector) std::vector<ImgRaster>; */
/*   %template(IntVector) std::vector<int>; */
//swig -c++ -I../.. -I.. -python -o pktools_wrap.cc pktools.i
//g++ -fPIC -std=c++11 -I../.. -I.. -I/usr/include/python2.7 -c pktools_wrap.cc $(python-config --cflags) -o pktools_wrap.o
//g++ -shared -v -nostartfiles -L/usr/local/lib pktools_wrap.o -limageClasses -lalgorithms -lgsl -ldl -lgdal $(python-config --ldflags) -o _pktools.so


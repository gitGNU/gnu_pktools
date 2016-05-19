//SWIG interface for ImgRaster
%module ImgRaster
%{
  #include "ImgRasterGdal.h"
  %}

//Parse the header file
%include "ImgRasterGdal.h"

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

%template(readData_char) ImgRasterGdal::readData<char>;
%template(readData_uchar) ImgRasterGdal::readData<unsigned char>;
%template(readData_ushort) ImgRasterGdal::readData<unsigned short>;
%template(readData_short) ImgRasterGdal::readData<short>;
%template(readData_uint) ImgRasterGdal::readData<unsigned int>;
%template(readData_int) ImgRasterGdal::readData<int>;
%template(readData_ulong) ImgRasterGdal::readData<unsigned long>;
%template(readData_long) ImgRasterGdal::readData<long>;
%template(readData_float) ImgRasterGdal::readData<float>;
%template(readData_double) ImgRasterGdal::readData<double>;

%template(readDataBlock_char) ImgRasterGdal::readDataBlock<char>;
%template(readDataBlock_uchar) ImgRasterGdal::readDataBlock<unsigned char>;
%template(readDataBlock_ushort) ImgRasterGdal::readDataBlock<unsigned short>;
%template(readDataBlock_short) ImgRasterGdal::readDataBlock<short>;
%template(readDataBlock_uint) ImgRasterGdal::readDataBlock<unsigned int>;
%template(readDataBlock_int) ImgRasterGdal::readDataBlock<int>;
%template(readDataBlock_ulong) ImgRasterGdal::readDataBlock<unsigned long>;
%template(readDataBlock_long) ImgRasterGdal::readDataBlock<long>;
%template(readDataBlock_float) ImgRasterGdal::readDataBlock<float>;
%template(readDataBlock_double) ImgRasterGdal::readDataBlock<double>;

%template(writeData_char) ImgRasterGdal::writeData<char>;
%template(writeData_uchar) ImgRasterGdal::writeData<unsigned char>;
%template(writeData_ushort) ImgRasterGdal::writeData<unsigned short>;
%template(writeData_short) ImgRasterGdal::writeData<short>;
%template(writeData_uint) ImgRasterGdal::writeData<unsigned int>;
%template(writeData_int) ImgRasterGdal::writeData<int>;
%template(writeData_ulong) ImgRasterGdal::writeData<unsigned long>;
%template(writeData_long) ImgRasterGdal::writeData<long>;
%template(writeData_float) ImgRasterGdal::writeData<float>;
%template(writeData_double) ImgRasterGdal::writeData<double>;

%template(writeDataBlock_char) ImgRasterGdal::writeDataBlock<char>;
%template(writeDataBlock_uchar) ImgRasterGdal::writeDataBlock<unsigned char>;
%template(writeDataBlock_ushort) ImgRasterGdal::writeDataBlock<unsigned short>;
%template(writeDataBlock_short) ImgRasterGdal::writeDataBlock<short>;
%template(writeDataBlock_uint) ImgRasterGdal::writeDataBlock<unsigned int>;
%template(writeDataBlock_int) ImgRasterGdal::writeDataBlock<int>;
%template(writeDataBlock_ulong) ImgRasterGdal::writeDataBlock<unsigned long>;
%template(writeDataBlock_long) ImgRasterGdal::writeDataBlock<long>;
%template(writeDataBlock_float) ImgRasterGdal::writeDataBlock<float>;
%template(writeDataBlock_double) ImgRasterGdal::writeDataBlock<double>;

/* swig -c++ -I../imageclasses -python -o ImgRaster_wrap.cc ImgRaster.i */
/* add following lines to ImgRaster_wrap.cc */
/* extern "C" */
/* { */
/* void *__dso_handle = 0; */
/* } */
/* g++ -fPIC -I.. -I../imageclasses -I/usr/include/python2.7 -c ImgRaster_wrap.cc $(python-config --cflags) -o ImgRaster_wrap.o */
/* g++ -v -nostartfiles -L/usr/local/lib ImgRaster_wrap.o -limageClasses -ldl -lgdal $(python-config --ldflags) -limageClasses -o _ImgRaster.so */

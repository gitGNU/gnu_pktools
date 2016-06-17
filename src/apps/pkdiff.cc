/**********************************************************************
pkdiff.cc: program to compare two raster image files
Copyright (C) 2008-2014 Pieter Kempeneers

This file is part of pktools

pktools is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

pktools is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with pktools.  If not, see <http://www.gnu.org/licenses/>.
***********************************************************************/
#include <assert.h>
#include "imageclasses/ImgRaster.h"
#include "imageclasses/ImgReaderOgr.h"
#include "imageclasses/ImgWriterOgr.h"
#include "base/Optionpk.h"
#include "algorithms/ConfusionMatrix.h"

/******************************************************************************/
/*! \page pkdiff pkdiff
 program to compare two raster image files
## SYNOPSIS

<code>
  Usage: pkdiff -i input -ref reference
</code>

<code>

  Options: [-ln layer] [-b band] [-cm] [-lr attribute] [-c name -r value]* [-nodata value]* [-m mask] [-msknodata value]*

  Advanced options:
       [-o output] [-f OGR format] [-lc attribute] [-bnd value [-hom] [-circ]] [-ct colortable] [-co NAME=VALUE]* 

</code>

\section pkdiff_description Description

The utility pkdiff compares two datasets. The reference can either be a raster or a vector, but the input must be a raster dataset.  
In case the reference is a raster dataset, a pixel by pixel comparison is performed. With no further options, the utility reports if the rasters are identical or different. If required, an output raster dataset can be written with a qualitative information per pixel: 0 (input=reference), 1 (input>reference) or 2 (input<reference). 
If, however, the reference is a vector dataset, it must consist of point features. Polygon features are automatically converted to the centroid points before analyzing. 

A typical use of the utility is to assess the accuracy of an input raster land cover map, based on a reference vector dataset. The reference dataset must contain an attribute (label) for each class. A confusion matrix is produced if the option -cm|--confusion is set. Here too, an output dataset can be written, which will be a vector dataset in this case. It contains the reference feature points with the extracted data value of the raster input dataset as a new attribute.
\section pkdiff_options Options
 - use either `-short` or `--long` options (both `--long=value` and `--long value` are supported)
 - short option `-h` shows basic options only, long option `--help` shows all options
|short|long|type|default|description|
|-----|----|----|-------|-----------|
 | i      | input                | std::string |       |Input raster dataset. | 
 | ref    | reference            | std::string |       |Reference (raster or vector) dataset | 
 | ln     | ln                   | std::string |       |Layer name(s) in sample. Leave empty to select all (for vector reference datasets only) | 
 | b      | band                 | short | 0     |Input (reference) raster band. Optionally, you can define different bands for input and reference bands respectively: -b 1 -b 0. | 
 | rmse   | rmse                 | bool | false |Report root mean squared error | 
 | reg    | reg                  | bool | false |Report linear regression (Input = c0+c1*Reference) | 
 | cm     | confusion            | bool | false |Create confusion matrix (to std out) | 
 | lr     | lref                 | std::string | label |Attribute name of the reference label (for vector reference datasets only) | 
 | c      | class                | std::string |       |List of class names. | 
 | r      | reclass              | short |       |List of class values (use same order as in classname option). | 
 | nodata | nodata               | double |       |No data value(s) in input or reference dataset are ignored | 
 | m      | mask                 | std::string |       |Use the first band of the specified file as a validity mask. Nodata values can be set with the option msknodata. | 
 | msknodata | msknodata            | double | 0     |Mask value(s) where image is invalid. Use negative value for valid data (example: use -t -1: if only -1 is valid value) | 
 | o      | output               | std::string |       |Output dataset (optional) | 
 | f      | f                    | std::string | SQLite |OGR format for output vector (for vector reference datasets only) | 
 | of     | oformat              | std::string | GTiff |Output image format (see also gdal_translate).| 
 | lc     | lclass               | std::string | class |Attribute name of the classified label (for vector reference datasets only) | 
 | cmf    | cmf                  | std::string | ascii |Format for confusion matrix (ascii or latex) | 
 | cmo    | cmo                  | std::string |       |Output file for confusion matrix | 
 | se95   | se95                 | bool | false |Report standard error for 95 confidence interval | 
 | bnd    | boundary             | short | 1     |Boundary for selecting the sample (for vector reference datasets only) | 
 | hom    | homogeneous          | bool | false |Only take regions with homogeneous boundary into account (for reference datasets only) | 
 | circ   | circular             | bool | false |Use circular boundary (for vector reference datasets only) | 
 | ct     | ct                   | std::string |       |Color table in ASCII format having 5 columns: id R G B ALFA (0: transparent, 255: solid). | 
 | co     | co                   | std::string |       |Creation option for output file. Multiple options can be specified. | 
 |        | commission           | short | 2     |Value for commission errors: input label < reference label | 
 | mem    | mem                  | unsigned long int | 0 |Buffer size (in MB) to read image data blocks in memory | 

Usage: pkdiff -i input -ref reference


Examples
========
Some examples how to use pkdiff can be found \ref examples_pkdiff "here"
**/

using namespace std;

int main(int argc, char *argv[])
{
  Optionpk<string> input_opt("i", "input", "Input raster dataset.");
  Optionpk<string> reference_opt("ref", "reference", "Reference (raster or vector) dataset");
  Optionpk<string> layer_opt("ln", "ln", "Layer name(s) in sample. Leave empty to select all (for vector reference datasets only)");
  Optionpk<string> mask_opt("m", "mask", "Use the first band of the specified file as a validity mask. Nodata values can be set with the option msknodata.");
  Optionpk<double> msknodata_opt("msknodata", "msknodata", "Mask value(s) where image is invalid. Use negative value for valid data (example: use -t -1: if only -1 is valid value)", 0);
  Optionpk<double> nodata_opt("nodata", "nodata", "No data value(s) in input or reference dataset are ignored");
  Optionpk<unsigned int> band_opt("b", "band", "Input (reference) raster band. Optionally, you can define different bands for input and reference bands respectively: -b 1 -b 0.", 0);
  Optionpk<bool> rmse_opt("rmse", "rmse", "Report root mean squared error", false);
  Optionpk<bool> regression_opt("reg", "reg", "Report linear regression (Input = c0+c1*Reference)", false);
  Optionpk<bool> confusion_opt("cm", "confusion", "Create confusion matrix (to std out)", false);
  Optionpk<string> cmformat_opt("cmf","cmf","Format for confusion matrix (ascii or latex)","ascii");
  Optionpk<string> cmoutput_opt("cmo","cmo","Output file for confusion matrix");
  Optionpk<bool> se95_opt("se95","se95","Report standard error for 95 confidence interval",false);
  Optionpk<string> labelref_opt("lr", "lref", "Attribute name of the reference label (for vector reference datasets only)", "label");
  Optionpk<string> classname_opt("c", "class", "List of class names."); 
  Optionpk<short> classvalue_opt("r", "reclass", "List of class values (use same order as in classname option)."); 
  Optionpk<string> output_opt("o", "output", "Output dataset (optional)");
  Optionpk<string> ogrformat_opt("f", "f", "OGR format for output vector (for vector reference datasets only)","SQLite");
  Optionpk<string>  oformat_opt("of", "oformat", "Output image format (see also gdal_translate).","GTiff");
  Optionpk<string> labelclass_opt("lc", "lclass", "Attribute name of the classified label (for vector reference datasets only)", "class");
  Optionpk<short> boundary_opt("bnd", "boundary", "Boundary for selecting the sample (for vector reference datasets only)", 1,1);
  Optionpk<bool> homogeneous_opt("hom", "homogeneous", "Only take regions with homogeneous boundary into account (for reference datasets only)", false,1);
  Optionpk<bool> disc_opt("circ", "circular", "Use circular boundary (for vector reference datasets only)", false,1);
  Optionpk<string> colorTable_opt("ct", "ct", "Color table in ASCII format having 5 columns: id R G B ALFA (0: transparent, 255: solid).");
  Optionpk<string> option_opt("co", "co", "Creation option for output file. Multiple options can be specified.");
  Optionpk<short> valueE_opt("\0", "correct", "Value for correct pixels", 0,2);
  Optionpk<short> valueO_opt("\0", "omission", "Value for omission errors: input label > reference label", 1,2);
  Optionpk<short> valueC_opt("\0", "commission", "Value for commission errors: input label < reference label", 2,1);
  Optionpk<unsigned long int>  memory_opt("mem", "mem", "Buffer size (in MB) to read image data blocks in memory",0,1);
  Optionpk<short> verbose_opt("v", "verbose", "Verbose level", 0,2);

  output_opt.setHide(1);
  ogrformat_opt.setHide(1);
  oformat_opt.setHide(1);
  labelclass_opt.setHide(1);
  boundary_opt.setHide(1);
  homogeneous_opt.setHide(1);
  disc_opt.setHide(1);
  colorTable_opt.setHide(1);
  option_opt.setHide(1);
  memory_opt.setHide(1);

  bool doProcess;//stop process when program was invoked with help option (-h --help)
  try{
    doProcess=input_opt.retrieveOption(argc,argv);
    reference_opt.retrieveOption(argc,argv);
    layer_opt.retrieveOption(argc,argv);
    band_opt.retrieveOption(argc,argv);
    rmse_opt.retrieveOption(argc,argv);
    regression_opt.retrieveOption(argc,argv);
    confusion_opt.retrieveOption(argc,argv);
    labelref_opt.retrieveOption(argc,argv);
    classname_opt.retrieveOption(argc,argv);
    classvalue_opt.retrieveOption(argc,argv);
    nodata_opt.retrieveOption(argc,argv);
    mask_opt.retrieveOption(argc,argv);
    msknodata_opt.retrieveOption(argc,argv);
    output_opt.retrieveOption(argc,argv);
    ogrformat_opt.retrieveOption(argc,argv);
    labelclass_opt.retrieveOption(argc,argv);
    cmformat_opt.retrieveOption(argc,argv);
    cmoutput_opt.retrieveOption(argc,argv);
    se95_opt.retrieveOption(argc,argv);
    boundary_opt.retrieveOption(argc,argv);
    homogeneous_opt.retrieveOption(argc,argv);
    disc_opt.retrieveOption(argc,argv);
    colorTable_opt.retrieveOption(argc,argv);
    option_opt.retrieveOption(argc,argv);
    // class_opt.retrieveOption(argc,argv);
    valueE_opt.retrieveOption(argc,argv);
    valueO_opt.retrieveOption(argc,argv);
    valueC_opt.retrieveOption(argc,argv);
    memory_opt.retrieveOption(argc,argv);
    verbose_opt.retrieveOption(argc,argv);
  }
  catch(string predefinedString){
    std::cout << predefinedString << std::endl;
    exit(0);
  }
  if(!doProcess){
    cout << endl;
    cout << "Usage: pkdiff -i input -ref reference" << endl;
    cout << endl;
    std::cout << "short option -h shows basic options only, use long option --help to show all options" << std::endl;
    exit(0);//help was invoked, stop processing
  }

  ImgRaster inputReader;
  ImgRaster maskReader;

  if(verbose_opt[0]){
    cout << "flag(s) set to";
    for(int iflag=0;iflag<nodata_opt.size();++iflag)
      cout << " " << nodata_opt[iflag];
    cout << endl;
  }

  if(input_opt.empty()){
    std::cerr << "No input file provided (use option -i). Use --help for help information" << std::endl;
    exit(0);
  }
  if(reference_opt.empty()){
    std::cerr << "No reference file provided (use option -ref). Use --help for help information" << std::endl;
    exit(0);
  }

  //band_opt[0] is for input
  //band_opt[1] is for reference
  if(band_opt.size()<2)
    band_opt.push_back(band_opt[0]);

  if(mask_opt.size())
    while(mask_opt.size()<input_opt.size())
      mask_opt.push_back(mask_opt[0]);
  vector<short> inputRange;
  vector<short> referenceRange;
  confusionmatrix::ConfusionMatrix cm;
  int nclass=0;
  map<string,short> classValueMap;
  vector<std::string> nameVector(255);//the inverse of the classValueMap
  vector<string> classNames;

  unsigned int ntotalValidation=0;
  unsigned int nflagged=0;
  Vector2d<unsigned int> resultClass;
  vector<float> user;
  vector<float> producer;
  vector<unsigned int> nvalidation;

  if(confusion_opt[0]){
    // if(class_opt.size()>1)
    //   inputRange=class_opt;
    // if(classvalue_opt.size()>1)
    //   inputRange=classvalue_opt;
    // else{
      try{
        if(verbose_opt[0])
          cout << "opening input image file " << input_opt[0] << endl;
        inputReader.open(input_opt[0],memory_opt[0]);//,imagicX_opt[0],imagicY_opt[0]);
      }
      catch(string error){
        cerr << error << endl;
        exit(1);
      }
      inputReader.getRange(inputRange,band_opt[0]);
      inputReader.close();
    // }
  
    for(int iflag=0;iflag<nodata_opt.size();++iflag){
      vector<short>::iterator fit;
      fit=find(inputRange.begin(),inputRange.end(),static_cast<short>(nodata_opt[iflag]));
      if(fit!=inputRange.end())
        inputRange.erase(fit);
    }
    nclass=inputRange.size();
    if(verbose_opt[0]){
      cout << "nclass (inputRange.size()): " << nclass << endl;
      cout << "input range: " << endl;
    }
    if(classname_opt.size()){
      assert(classname_opt.size()==classvalue_opt.size());
      for(int iclass=0;iclass<classname_opt.size();++iclass){
        classValueMap[classname_opt[iclass]]=classvalue_opt[iclass];
        assert(classvalue_opt[iclass]<nameVector.size());
        nameVector[classvalue_opt[iclass]]=classname_opt[iclass];
      }
    }
    // nclass=classValueMap.size();
    for(int rc=0;rc<inputRange.size();++rc){
      classNames.push_back(type2string(inputRange[rc]));
      if(verbose_opt[0])
        cout << inputRange[rc] << endl;
    }
    cm.setClassNames(classNames);
    if(verbose_opt[0]){
      cout << "class names: " << endl;
      for(int iclass=0;iclass<cm.nClasses();++iclass)
        cout << iclass << " " << cm.getClass(iclass) << endl;
    }
    resultClass.resize(nclass,nclass);
    user.resize(nclass);
    producer.resize(nclass);
    nvalidation.resize(nclass);
    //initialize
    for(int rc=0;rc<nclass;++rc){
      for(int ic=0;ic<nclass;++ic)
        resultClass[rc][ic]=0;
      nvalidation[rc]=0;
    }
  }
  
  bool isDifferent=false;
  bool refIsRaster=false;

  ImgReaderOgr referenceReaderOgr;
  try{
    referenceReaderOgr.open(reference_opt[0]);
    referenceReaderOgr.close();
  }
  catch(string errorString){
    //todo: sampleIsRaster will not work from GDAL 2.0!!?? (unification of driver for raster and vector datasets)
    refIsRaster=true;
  }
  const char* pszMessage;
  void* pProgressArg=NULL;
  GDALProgressFunc pfnProgress=GDALTermProgress;
  float progress=0;
  // if(reference_opt[0].find(".shp")!=string::npos){
  if(!refIsRaster){
    for(int iinput=0;iinput<input_opt.size();++iinput){
      if(verbose_opt[0])
	cout << "Processing input " << input_opt[iinput] << endl;
      if(output_opt.size())
        assert(reference_opt.size()==output_opt.size());
      for(int iref=0;iref<reference_opt.size();++iref){
	cout << "reference " << reference_opt[iref] << endl;
        // assert(reference_opt[iref].find(".shp")!=string::npos);
        try{
          inputReader.open(input_opt[iinput],memory_opt[0]);//,imagicX_opt[0],imagicY_opt[0]);
          if(mask_opt.size()){
            maskReader.open(mask_opt[iinput],memory_opt[0]);
            assert(inputReader.nrOfCol()==maskReader.nrOfCol());
            assert(inputReader.nrOfRow()==maskReader.nrOfRow());
          }
          referenceReaderOgr.open(reference_opt[iref]);
        }
        catch(string error){
          cerr << error << endl;
          exit(1);
        }
        if(confusion_opt[0])
          referenceRange=inputRange;

        ImgWriterOgr ogrWriter;
	if(output_opt.size()){
	  try{
	    ogrWriter.open(output_opt[iref],ogrformat_opt[0]);
	  }
	  catch(string error){
	    cerr << error << endl;
	    exit(1);
	  }
	}
	int nlayer=referenceReaderOgr.getDataSource()->GetLayerCount();
	for(int ilayer=0;ilayer<nlayer;++ilayer){
	  progress=0;
	  OGRLayer *readLayer=referenceReaderOgr.getLayer(ilayer);
	  //	  readLayer = referenceReaderOgr.getDataSource()->GetLayer(ilayer);
	  string currentLayername=readLayer->GetName();
	  if(layer_opt.size())
	    if(find(layer_opt.begin(),layer_opt.end(),currentLayername)==layer_opt.end())
	      continue;
	  if(!verbose_opt[0])
	    pfnProgress(progress,pszMessage,pProgressArg);
	  else
	    cout << "processing layer " << readLayer->GetName() << endl;

	  readLayer->ResetReading();
	  OGRLayer *writeLayer;
	  if(output_opt.size()){
	    if(verbose_opt[0])
	      cout << "creating output vector file " << output_opt[0] << endl;
	    // assert(output_opt[0].find(".shp")!=string::npos);
	    char     **papszOptions=NULL;
	    if(verbose_opt[0])
	      cout << "creating layer: " << readLayer->GetName() << endl;
	    // if(ogrWriter.createLayer(layername, referenceReaderOgr.getProjection(ilayer), referenceReaderOgr.getGeometryType(ilayer), papszOptions)==NULL)
	    writeLayer=ogrWriter.createLayer(readLayer->GetName(), referenceReaderOgr.getProjection(ilayer), wkbPoint, papszOptions);
	    assert(writeLayer);
	    if(verbose_opt[0]){
	      cout << "created layer" << endl;
	      cout << "copy fields from " << reference_opt[iref] << endl;
	    }
	    ogrWriter.copyFields(referenceReaderOgr,ilayer,ilayer);
	    //create extra field for classified label
	    short theDim=boundary_opt[0];
	    for(int windowJ=-theDim/2;windowJ<(theDim+1)/2;++windowJ){
	      for(int windowI=-theDim/2;windowI<(theDim+1)/2;++windowI){
		if(disc_opt[0]&&(windowI*windowI+windowJ*windowJ>(theDim/2)*(theDim/2)))
		  continue;
		ostringstream fs;
		if(theDim>1)
		  fs << labelclass_opt[0] << "_" << windowJ << "_" << windowI;
		else
		  fs << labelclass_opt[0];
		if(verbose_opt[0])
		  cout << "creating field " << fs.str() << endl;
		ogrWriter.createField(fs.str(),OFTInteger,ilayer);
	      }
	    }
	  }
	  OGRFeature *readFeature;
	  OGRFeature *writeFeature;
	  int isample=0;
	  unsigned int nfeatureInLayer=readLayer->GetFeatureCount();
	  unsigned int ifeature=0;
	  while( (readFeature = readLayer->GetNextFeature()) != NULL ){
	    if(verbose_opt[0])
	      cout << "sample " << ++isample << endl;
	    //get x and y from readFeature
	    double x,y;
	    OGRGeometry *poGeometry;
	    OGRPoint centroidPoint;
	    OGRPoint *poPoint;
	    poGeometry = readFeature->GetGeometryRef();
	    // assert( poGeometry != NULL && wkbFlatten(poGeometry->getGeometryType()) == wkbPoint );
	    if(poGeometry==NULL)
	      continue;
	    else if(wkbFlatten(poGeometry->getGeometryType()) == wkbMultiPolygon){
	      OGRMultiPolygon readPolygon = *((OGRMultiPolygon *) poGeometry);
	      readPolygon = *((OGRMultiPolygon *) poGeometry);
	      readPolygon.Centroid(&centroidPoint);
	      poPoint=&centroidPoint;
	    }
	    else if(wkbFlatten(poGeometry->getGeometryType()) == wkbPolygon){
	      OGRPolygon readPolygon=*((OGRPolygon *) poGeometry);
	      readPolygon.Centroid(&centroidPoint);
	      poPoint=&centroidPoint;
	    }
	    else if(wkbFlatten(poGeometry->getGeometryType()) == wkbPoint )
	      poPoint = (OGRPoint *) poGeometry;
	    else{
	      std::cerr << "Warning: skipping feature (not of type point or polygon)" << std::endl;
	      continue;
	    }
	    x=poPoint->getX();
	    y=poPoint->getY();
	    double inputValue;
	    vector<double> inputValues;
	    bool isHomogeneous=true;
	    short maskValue;
	    short outputValue;
	    //read referenceValue from feature
	    unsigned short referenceValue;
	    string referenceClassName;
	    if(classValueMap.size()){
	      referenceClassName=readFeature->GetFieldAsString(readFeature->GetFieldIndex(labelref_opt[0].c_str()));
	      referenceValue=classValueMap[referenceClassName];
	    }
	    else
	      referenceValue=readFeature->GetFieldAsInteger(readFeature->GetFieldIndex(labelref_opt[0].c_str()));
	    if(verbose_opt[0])
	      cout << "reference value: " << referenceValue << endl;
	    
	    bool pixelFlagged=false;
	    bool maskFlagged=false;
	    for(int iflag=0;iflag<nodata_opt.size();++iflag){
	      if(referenceValue==nodata_opt[iflag])
		pixelFlagged=true;
	    }
	    if(pixelFlagged)
	      continue;
	    double i_centre,j_centre;
	    //input reader is georeferenced!
	    inputReader.geo2image(x,y,i_centre,j_centre);
	    //       else{
	    //         i_centre=x;
	    //         j_centre=y;
	    //       }
	    //nearest neighbour
	    j_centre=static_cast<unsigned int>(j_centre);
	    i_centre=static_cast<unsigned int>(i_centre);
	    //check if j_centre is out of bounds
	    if(static_cast<unsigned int>(j_centre)<0||static_cast<unsigned int>(j_centre)>=inputReader.nrOfRow())
	      continue;
	    //check if i_centre is out of bounds
	    if(static_cast<unsigned int>(i_centre)<0||static_cast<unsigned int>(i_centre)>=inputReader.nrOfCol())
	      continue;

	    if(output_opt.size()){
	      writeFeature = OGRFeature::CreateFeature(writeLayer->GetLayerDefn());
	      assert(readFeature);
	      int nfield=readFeature->GetFieldCount();
	      writeFeature->SetGeometry(poPoint);
	      if(verbose_opt[0])
	      	cout << "copying fields from " << reference_opt[0] << endl;
	      assert(readFeature);
	      assert(writeFeature);
	      vector<int> panMap(nfield);
	      vector< int>::iterator panit=panMap.begin();
	      for(unsigned int ifield=0;ifield<nfield;++ifield)
	      	panMap[ifield]=ifield;
	      writeFeature->SetFieldsFrom(readFeature,&(panMap[0]));
	      // if(writeFeature->SetFrom(readFeature)!= OGRERR_NONE)
	      // 	cerr << "writing feature failed" << endl;
	      // if(verbose_opt[0])
	      // 	cout << "feature written" << endl;
	    }
	    bool windowAllFlagged=true;
	    bool windowHasFlag=false;
	    short theDim=boundary_opt[0];
	    for(int windowJ=-theDim/2;windowJ<(theDim+1)/2;++windowJ){
	      for(int windowI=-theDim/2;windowI<(theDim+1)/2;++windowI){
		if(disc_opt[0]&&(windowI*windowI+windowJ*windowJ>(theDim/2)*(theDim/2)))
		  continue;
		int j=j_centre+windowJ;
		//check if j is out of bounds
		if(static_cast<unsigned int>(j)<0||static_cast<unsigned int>(j)>=inputReader.nrOfRow())
		  continue;
		int i=i_centre+windowI;
		//check if i is out of bounds
		if(static_cast<unsigned int>(i)<0||static_cast<unsigned int>(i)>=inputReader.nrOfCol())
		  continue;
		if(verbose_opt[0])
		  cout << setprecision(12) << "reading image value at x,y " << x << "," << y << " (" << i << "," << j << "), ";
		inputReader.readData(inputValue,i,j,band_opt[0]);
		inputValues.push_back(inputValue);
		if(inputValues.back()!=*(inputValues.begin()))
		  isHomogeneous=false;
		if(verbose_opt[0])
		  cout << "input value: " << inputValue << endl;
		pixelFlagged=false;
		for(int iflag=0;iflag<nodata_opt.size();++iflag){
		  if(inputValue==nodata_opt[iflag]){
		    pixelFlagged=true;
		    break;
		  }
		}
		maskFlagged=false;//(msknodata_opt[ivalue]>=0)?false:true;
		if(mask_opt.size()){
		  maskReader.readData(maskValue,i,j,0);
		  for(int ivalue=0;ivalue<msknodata_opt.size();++ivalue){
		    if(msknodata_opt[ivalue]>=0){//values set in msknodata_opt are invalid
		      if(maskValue==msknodata_opt[ivalue]){
			maskFlagged=true;
			break;
		      }
		    }
		    else{//only values set in msknodata_opt are valid
		      if(maskValue!=-msknodata_opt[ivalue])
			maskFlagged=true;
		      else{
			maskFlagged=false;
			break;
		      }
		    }
		  }
		}
		pixelFlagged=pixelFlagged||maskFlagged;
		if(pixelFlagged)
		  windowHasFlag=true;
		else
		  windowAllFlagged=false;//at least one good pixel in neighborhood
	      }
	    }
	    //at this point we know the values for the entire window

	    if(homogeneous_opt[0]){//only centre pixel
	      //flag if not all pixels are homogeneous or if at least one pixel flagged
          
	      if(!windowHasFlag&&isHomogeneous){
		if(output_opt.size())
		  writeFeature->SetField(labelclass_opt[0].c_str(),static_cast<int>(inputValue));
		if(confusion_opt[0]){
		  ++ntotalValidation;
		  if(classValueMap.size()){
		    assert(inputValue<nameVector.size());
		    string className=nameVector[static_cast<unsigned short>(inputValue)];
		    cm.incrementResult(type2string<short>(classValueMap[referenceClassName]),type2string<short>(classValueMap[className]),1);
		  }
		  else{
		    int rc=distance(referenceRange.begin(),find(referenceRange.begin(),referenceRange.end(),static_cast<unsigned short>(referenceValue)));
		    int ic=distance(inputRange.begin(),find(inputRange.begin(),inputRange.end(),static_cast<unsigned short>(inputValue)));
		    assert(rc<nclass);
		    assert(ic<nclass);
		    ++nvalidation[rc];
		    ++resultClass[rc][ic];
		    if(verbose_opt[0]>1)
		      cout << "increment: " << rc << " " << referenceRange[rc] << " " << ic << " " << inputRange[ic] << endl;
		    cm.incrementResult(cm.getClass(rc),cm.getClass(ic),1);
		  }
		}
		if(inputValue==referenceValue){//correct
		  outputValue=valueE_opt[0];
		  if(nodata_opt.size()){
		    if(valueE_opt[0]==nodata_opt[0])
		      outputValue=inputValue;
		  }
		}
		else if(inputValue>referenceValue)//1=forest,2=non-forest
		  outputValue=valueO_opt[0];//omission error
		else
		  outputValue=valueC_opt[0];//commission error
	      }
	    }
	    else{
	      for(int windowJ=-theDim/2;windowJ<(theDim+1)/2;++windowJ){
		for(int windowI=-theDim/2;windowI<(theDim+1)/2;++windowI){
		  if(disc_opt[0]&&(windowI*windowI+windowJ*windowJ>(theDim/2)*(theDim/2)))
		    continue;
		  int j=j_centre+windowJ;
		  //check if j is out of bounds
		  if(static_cast<unsigned int>(j)<0||static_cast<unsigned int>(j)>=inputReader.nrOfRow())
		    continue;
		  int i=i_centre+windowI;
		  //check if i is out of bounds
		  if(static_cast<unsigned int>(i)<0||static_cast<unsigned int>(i)>=inputReader.nrOfCol())
		    continue;
		  if(!windowAllFlagged){
		    ostringstream fs;
		    if(theDim>1)
		      fs << labelclass_opt[0] << "_" << windowJ << "_" << windowI;
		    else
		      fs << labelclass_opt[0];
		    if(output_opt.size())
		      writeFeature->SetField(fs.str().c_str(),inputValue);
		    if(!windowJ&&!windowI){//centre pixel
		      if(confusion_opt[0]){
			++ntotalValidation;
			if(classValueMap.size()){
			  assert(inputValue<nameVector.size());
			  string className=nameVector[static_cast<unsigned short>(inputValue)];
			  cm.incrementResult(type2string<short>(classValueMap[referenceClassName]),type2string<short>(classValueMap[className]),1);
			}
			else{
			  int rc=distance(referenceRange.begin(),find(referenceRange.begin(),referenceRange.end(),static_cast<unsigned short>(referenceValue)));
			  int ic=distance(inputRange.begin(),find(inputRange.begin(),inputRange.end(),static_cast<unsigned short>(inputValue)));
			  if(rc>=nclass)
			    continue;
			  if(ic>=nclass)
			    continue;
			  // assert(rc<nclass);
			  // assert(ic<nclass);
			  ++nvalidation[rc];
			  ++resultClass[rc][ic];
			  if(verbose_opt[0]>1)
			    cout << "increment: " << rc << " " << referenceRange[rc] << " " << ic << " " << inputRange[ic] << endl;
			  cm.incrementResult(cm.getClass(rc),cm.getClass(ic),1);
			}
		      }
		      if(inputValue==referenceValue){//correct
			outputValue=valueE_opt[0];
			if(nodata_opt.size()){
			  if(valueE_opt[0]==nodata_opt[0])
			    outputValue=inputValue;
			}
		      }
		      else if(inputValue>referenceValue)//1=forest,2=non-forest
			outputValue=valueO_opt[0];//omission error
		      else
			outputValue=valueC_opt[0];//commission error
		    }
		  }
		}
	      }
	    }
	    if(output_opt.size()){
	      if(!windowAllFlagged){
		if(verbose_opt[0])
		  cout << "creating feature" << endl;
		if(writeLayer->CreateFeature( writeFeature ) != OGRERR_NONE ){
		  string errorString="Failed to create feature in OGR vector file";
		  throw(errorString);
		}
	      }
	      OGRFeature::DestroyFeature( writeFeature );
	    }
	    ++ifeature;
	    progress=static_cast<float>(ifeature+1)/nfeatureInLayer;
	    pfnProgress(progress,pszMessage,pProgressArg);
	  }//next feature
	}//next layer
        if(output_opt.size())
          ogrWriter.close();
        referenceReaderOgr.close();
        inputReader.close();
        if(mask_opt.size())
          maskReader.close();
      }//next reference
    }//next input
    pfnProgress(1.0,pszMessage,pProgressArg);
  }//reference is OGR vector
  else{//reference is GDAL raster
    ImgRaster gdalWriter;
    try{
      inputReader.open(input_opt[0],memory_opt[0]);
      if(mask_opt.size())
        maskReader.open(mask_opt[0],memory_opt[0]);
      if(output_opt.size()){
        if(verbose_opt[0])
          cout << "opening output image " << output_opt[0] << endl;
        if(option_opt.findSubstring("INTERLEAVE=")==option_opt.end()){
          string theInterleave="INTERLEAVE=";
          theInterleave+=inputReader.getInterleave();
          option_opt.push_back(theInterleave);
        }
        gdalWriter.open(output_opt[0],inputReader.nrOfCol(),inputReader.nrOfRow(),1,inputReader.getDataType(),oformat_opt[0],memory_opt[0],option_opt);
	if(nodata_opt.size())
	  gdalWriter.GDALSetNoDataValue(nodata_opt[0]);
	gdalWriter.copyGeoTransform(inputReader);
        if(colorTable_opt.size())
          gdalWriter.setColorTable(colorTable_opt[0]);
        else if(inputReader.getColorTable()!=NULL){
          if(verbose_opt[0])
            cout << "set colortable from input image" << endl;
          gdalWriter.setColorTable(inputReader.getColorTable());
        }
      }
      else if(verbose_opt[0])
        cout << "no output image defined" << endl;
        
    }
    catch(string error){
      cout << error << endl;
      exit(2);
    }
    //todo: support different data types!
    vector<double> lineInput(inputReader.nrOfCol());
    vector<double> lineMask(maskReader.nrOfCol());
    vector<double> lineOutput;
    vector<double> bufferInput;//for regression
    vector<double> bufferReference;//for regression
    if(output_opt.size())
      lineOutput.resize(inputReader.nrOfCol());

    unsigned int irow=0;
    unsigned int icol=0;
    double oldreferencerow=-1;
    double oldmaskrow=-1;
    ImgRaster referenceReaderGdal;
    try{
      referenceReaderGdal.open(reference_opt[0]);//,rmagicX_opt[0],rmagicY_opt[0]);
    }
    catch(string error){
      cerr << error << endl;
      exit(1);
    }
    if(inputReader.isGeoRef()){
      assert(referenceReaderGdal.isGeoRef());
      if(inputReader.getProjection()!=referenceReaderGdal.getProjection())
        cerr << "Warning: projection of input image and reference image are different" << endl;
    }
    vector<double> lineReference(referenceReaderGdal.nrOfCol());
    if(confusion_opt[0]){
      referenceReaderGdal.getRange(referenceRange,band_opt[1]);
      for(int iflag=0;iflag<nodata_opt.size();++iflag){
        vector<short>::iterator fit;
        fit=find(referenceRange.begin(),referenceRange.end(),static_cast<unsigned short>(nodata_opt[iflag]));
        if(fit!=referenceRange.end())
          referenceRange.erase(fit);
      }
      if(verbose_opt[0]){
        cout << "reference range: " << endl;
        for(int rc=0;rc<referenceRange.size();++rc)
          cout << referenceRange[rc] << endl;
      }
      if(referenceRange.size()!=inputRange.size()){
        if(confusion_opt[0]||output_opt.size()){
          cout << "reference range is not equal to input range!" << endl;
          cout << "Kappa: " << 0 << endl;    
          cout << "total weighted: " << 0 << endl;
        }
        else
          cout << "reference range is not equal to input range!" << endl;
        cout << input_opt[0] << " and " << reference_opt[0] << " are different" << endl;
        exit(1);
      }
    }
    double rmse=0;
    // for(irow=0;irow<inputReader.nrOfRow()&&!isDifferent;++irow){
    for(irow=0;irow<inputReader.nrOfRow();++irow){
      //read line in lineInput, lineReference and lineMask
      inputReader.readData(lineInput,irow,band_opt[0]);
      double x,y;//geo coordinates
      double ireference,jreference;//image coordinates in reference image
      double imask,jmask;//image coordinates in mask image
      for(icol=0;icol<inputReader.nrOfCol();++icol){
        //find col in reference
        inputReader.image2geo(icol,irow,x,y);
        referenceReaderGdal.geo2image(x,y,ireference,jreference);
        if(ireference<0||ireference>=referenceReaderGdal.nrOfCol()){
	  if(rmse_opt[0]||regression_opt[0])
	    continue;
	  else{
	    cerr << ireference << " out of reference range!" << endl;
	    cerr << x << " " << y << " " << icol << " " << irow << endl;
	    cerr << x << " " << y << " " << ireference << " " << jreference << endl;
	    exit(1);
	  }
        }
        if(jreference!=oldreferencerow){
          if(jreference<0||jreference>=referenceReaderGdal.nrOfRow()){
	    if(rmse_opt[0]||regression_opt[0])
	      continue;
	    else{
	      cerr << jreference << " out of reference range!" << endl;
	      cerr << x << " " << y << " " << icol << " " << irow << endl;
	      cerr << x << " " << y << " " << ireference << " " << jreference << endl;
	      exit(1);
	    }
          }
          else{
            referenceReaderGdal.readData(lineReference,static_cast<unsigned int>(jreference),band_opt[1]);
            oldreferencerow=jreference;
          }
        }
        bool flagged=false;
        for(int iflag=0;iflag<nodata_opt.size();++iflag){
          if((lineInput[icol]==nodata_opt[iflag])||(lineReference[ireference]==nodata_opt[iflag])){
            if(output_opt.size())
              lineOutput[icol]=nodata_opt[iflag];
            flagged=true;
            break;
          }
        }
        if(mask_opt.size()){
	  maskReader.geo2image(x,y,imask,jmask);
	  if(jmask>=0&&jmask<maskReader.nrOfRow()){
	    if(jmask!=oldmaskrow)
	      maskReader.readData(lineMask,jmask);
	    for(int ivalue=0;ivalue<msknodata_opt.size();++ivalue){
	      if(lineMask[icol]==msknodata_opt[ivalue]){
		flagged=true;
		break;
	      }
	    }
          }
        }
        if(!flagged){
	  if(rmse_opt[0]){//divide by image size to prevent overflow. At the end we need to take care about flagged pixels by normalizing...
	    rmse+=static_cast<double>(lineInput[icol]-lineReference[ireference])*(lineInput[icol]-lineReference[ireference])/inputReader.nrOfCol()/inputReader.nrOfRow();
	  }
	  else if(regression_opt[0]){
	    bufferInput.push_back(lineInput[icol]);
	    bufferReference.push_back(lineReference[ireference]);
	  }

          if(confusion_opt[0]){
            ++ntotalValidation;
            int rc=distance(referenceRange.begin(),find(referenceRange.begin(),referenceRange.end(),lineReference[ireference]));
            int ic=distance(inputRange.begin(),find(inputRange.begin(),inputRange.end(),lineInput[icol]));
            assert(rc<nclass);
            assert(ic<nclass);
            ++nvalidation[rc];
            ++resultClass[rc][ic];
            if(verbose_opt[0]>1)
              cout << "increment: " << rc << " " << referenceRange[rc] << " " << ic << " " << inputRange[ic] << endl;
            cm.incrementResult(cm.getClass(rc),cm.getClass(ic),1);
          }
          if(lineInput[icol]==lineReference[ireference]){//correct
            if(output_opt.size()){
	      lineOutput[icol]=valueE_opt[0];
	      if(nodata_opt.size()){
		if(valueE_opt[0]==nodata_opt[0])
		  lineOutput[icol]=lineInput[icol];
	      }
            }
          }
          else{//error
            if(output_opt.empty()&&!confusion_opt[0]&&!rmse_opt[0]&&!regression_opt[0]){
              isDifferent=true;
              break;
            }
            if(output_opt.size()){
              if(lineInput[icol]>lineReference[ireference])
		lineOutput[icol]=valueO_opt[0];//omission error
	      else
		lineOutput[icol]=valueC_opt[0];//commission error
            }
          }
	}
        else{
          ++nflagged;
          if(output_opt.size()){
	    if(nodata_opt.size())
	      lineOutput[icol]=nodata_opt[0];
	    else //should never occur?
	      lineOutput[icol]=0;
	  }	      
        }
      }
      if(output_opt.size()){
        try{
          gdalWriter.writeData(lineOutput,irow);
        }
        catch(string errorstring){
          cerr << "lineOutput.size(): " << lineOutput.size() << endl;
          cerr << "gdalWriter.nrOfCol(): " << gdalWriter.nrOfCol() << endl;
          cerr << errorstring << endl;
          exit(1);
        }
      }
      else if(isDifferent&&!confusion_opt[0]&&!rmse_opt[0]&&!regression_opt[0]){//we can break off here, files are different...
        if(!verbose_opt[0])
          pfnProgress(1.0,pszMessage,pProgressArg);
        break;
      }
      progress=static_cast<float>(irow+1.0)/inputReader.nrOfRow();
      if(!verbose_opt[0])
        pfnProgress(progress,pszMessage,pProgressArg);
    }
    if(output_opt.size())
      gdalWriter.close();
    else if(!confusion_opt[0]){
      if(rmse_opt[0]){
	double normalization=1.0*inputReader.nrOfCol()*inputReader.nrOfRow()/(inputReader.nrOfCol()*inputReader.nrOfRow()-nflagged);
	if(verbose_opt[0]){
	  cout << "normalization: " << normalization << endl;
	  cout << "rmse before sqrt and normalization: " << rmse << endl;
	}
	cout << "--rmse " << sqrt(rmse/normalization) << endl;
      }
      else if(regression_opt[0]){
	double err=0;
	double c0=0;
	double c1=1;
	statfactory::StatFactory stat;
	if(bufferInput.size()&&bufferReference.size()){
	  err=stat.linear_regression_err(bufferInput,bufferReference,c0,c1);
	}
	if(verbose_opt[0]){
	  cout << "bufferInput.size(): " << bufferInput.size() << endl;
	  cout << "bufferReference.size(): " << bufferReference.size() << endl;
	  double theMin=0;
	  double theMax=0;
	  stat.minmax(bufferInput,bufferInput.begin(),bufferInput.end(),theMin,theMax);
	  cout << "min,  max input: " << theMin << ", " << theMax << endl;
	  theMin=0;
	  theMax=0;
	  stat.minmax(bufferReference,bufferReference.begin(),bufferReference.end(),theMin,theMax);
	  cout << "min,  max reference: " << theMin << ", " << theMax << endl;
	}
	cout << "--c0 " << c0 << "--c1 " << c1 << " --rmse: " << err << endl;
	
      }
      else if(isDifferent)
        cout << input_opt[0] << " and " << reference_opt[0] << " are different" << endl;
      else
        cout << input_opt[0] << " and " << reference_opt[0] << " are identical" << endl;
    }
    referenceReaderGdal.close();
    inputReader.close();
    if(mask_opt.size())
      maskReader.close();
  }//raster dataset

  if(confusion_opt[0]){
    cm.setFormat(cmformat_opt[0]);
    cm.reportSE95(se95_opt[0]);
    ofstream outputFile;
    if(cmoutput_opt.size()){
      outputFile.open(cmoutput_opt[0].c_str(),ios::out);
      outputFile << cm << endl;
    }
    else
      cout << cm << endl;
    // cout << "class #samples userAcc prodAcc" << endl;
    // double se95_ua=0;
    // double se95_pa=0;
    // double se95_oa=0;
    // double dua=0;
    // double dpa=0;
    // double doa=0;
    // for(int iclass=0;iclass<cm.nClasses();++iclass){
    //   dua=cm.ua_pct(classNames[iclass],&se95_ua);
    //   dpa=cm.pa_pct(classNames[iclass],&se95_pa);
    //   cout << cm.getClass(iclass) << " " << cm.nReference(cm.getClass(iclass)) << " " << dua << " (" << se95_ua << ")" << " " << dpa << " (" << se95_pa << ")" << endl;
    // }
    // doa=cm.oa(&se95_oa);
    // cout << "Kappa: " << cm.kappa() << endl;
    // cout << "Overall Accuracy: " << 100*doa << " (" << 100*se95_oa << ")"  << endl;
  }
}

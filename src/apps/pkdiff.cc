/**********************************************************************
pkdiff.cc: program to compare two raster image files
Copyright (C) 2008-2012 Pieter Kempeneers

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
#include "imageclasses/ImgReaderGdal.h"
#include "imageclasses/ImgWriterGdal.h"
#include "imageclasses/ImgReaderOgr.h"
#include "imageclasses/ImgWriterOgr.h"
#include "Optionpk.h"
#include "algorithms/ConfusionMatrix.h"

int main(int argc, char *argv[])
{
  Optionpk<string> input_opt("i", "input", "Input image file.");
  Optionpk<string> reference_opt("ref", "reference", "Reference image file");
  Optionpk<string> output_opt("o", "output", "Output image file. Default is empty: no output image, only report difference or identical.");
  Optionpk<string> mask_opt("m", "mask", "Mask image file. A single mask is supported only, but several mask values can be used. See also msknodata option. (default is empty)");
  Optionpk<int> masknodata_opt("msknodata", "msknodata", "Mask value(s) where image is invalid. Use negative value for valid data (example: use -t -1: if only -1 is valid value)", 0);
  Optionpk<string> colorTable_opt("ct", "ct", "color table (file with 5 columns: id R G B ALFA (0: transparent, 255: solid)");
  Optionpk<short> valueE_opt("\0", "correct", "Value for correct pixels (0)", 0,1);
  Optionpk<short> valueO_opt("\0", "omission", "Value for omission errors: input label > reference label (default value is 1)", 1,1);
  Optionpk<short> valueC_opt("\0", "commission", "Value for commission errors: input label < reference label (default value is 2)", 2,1);
  Optionpk<short> nodata_opt("nodata", "nodata", "No value flag(s)");
  Optionpk<short> band_opt("b", "band", "Band to extract (0)", 0);
  Optionpk<bool> confusion_opt("cm", "confusion", "create confusion matrix (to std out) (default value is 0)", false);
  Optionpk<string> labelref_opt("lr", "lref", "name of the reference label in case reference is shape file(default is label)", "label");
  Optionpk<string> labelclass_opt("lc", "lclass", "name of the classified label in case output is shape file (default is class)", "class");
  Optionpk<short> boundary_opt("bnd", "boundary", "boundary for selecting the sample (default: 1)", 1,1);
  Optionpk<bool> disc_opt("\0", "circular", "use circular disc kernel boundary)", false,1);
  Optionpk<bool> homogeneous_opt("hom", "homogeneous", "only take homogeneous regions into account", false,1);
  Optionpk<string> option_opt("co", "co", "Creation option for output file. Multiple options can be specified.");
  Optionpk<string> classname_opt("c", "class", "list of class names."); 
  Optionpk<short> classvalue_opt("r", "reclass", "list of class values (use same order as in classname opt."); 
  Optionpk<short> verbose_opt("v", "verbose", "verbose (default value is 0)", 0);

  bool doProcess;//stop process when program was invoked with help option (-h --help)
  try{
    doProcess=input_opt.retrieveOption(argc,argv);
    output_opt.retrieveOption(argc,argv);
    option_opt.retrieveOption(argc,argv);
    reference_opt.retrieveOption(argc,argv);
    mask_opt.retrieveOption(argc,argv);
    colorTable_opt.retrieveOption(argc,argv);
    valueE_opt.retrieveOption(argc,argv);
    valueO_opt.retrieveOption(argc,argv);
    valueC_opt.retrieveOption(argc,argv);
    nodata_opt.retrieveOption(argc,argv);
    masknodata_opt.retrieveOption(argc,argv);
    band_opt.retrieveOption(argc,argv);
    confusion_opt.retrieveOption(argc,argv);
    labelref_opt.retrieveOption(argc,argv);
    labelclass_opt.retrieveOption(argc,argv);
    // class_opt.retrieveOption(argc,argv);
    boundary_opt.retrieveOption(argc,argv);
    disc_opt.retrieveOption(argc,argv);
    homogeneous_opt.retrieveOption(argc,argv);
    classname_opt.retrieveOption(argc,argv);
    classvalue_opt.retrieveOption(argc,argv);
    verbose_opt.retrieveOption(argc,argv);
  }
  catch(string predefinedString){
    std::cout << predefinedString << std::endl;
    exit(0);
  }
  if(!doProcess){
    std::cout << "short option -h shows basic options only, use long option --help to show all options" << std::endl;
    exit(0);//help was invoked, stop processing
  }

  ImgReaderGdal inputReader;
  ImgReaderGdal maskReader;

  if(verbose_opt[0]){
    cout << "flag(s) set to";
    for(int iflag=0;iflag<nodata_opt.size();++iflag)
      cout << " " << nodata_opt[iflag];
    cout << endl;
  }

  assert(input_opt.size());
  assert(reference_opt.size());
  if(mask_opt.size())
    while(mask_opt.size()<input_opt.size())
      mask_opt.push_back(mask_opt[0]);
  vector<short> inputRange;
  vector<short> referenceRange;
  ConfusionMatrix cm;
  int nclass=0;
  map<string,short> classValueMap;
  vector<std::string> nameVector(255);//the inverse of the classValueMap
  vector<string> classNames;
  if(confusion_opt[0]){
    // if(class_opt.size()>1)
    //   inputRange=class_opt;
    // if(classvalue_opt.size()>1)
    //   inputRange=classvalue_opt;
    // else{
      try{
        if(verbose_opt[0])
          cout << "opening input image file " << input_opt[0] << endl;
        inputReader.open(input_opt[0]);//,imagicX_opt[0],imagicY_opt[0]);
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
      fit=find(inputRange.begin(),inputRange.end(),nodata_opt[iflag]);
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
  }

  unsigned int ntotalValidation=0;
  unsigned int nflagged=0;
  Vector2d<int> resultClass(nclass,nclass);
  vector<float> user(nclass);
  vector<float> producer(nclass);
  vector<unsigned int> nvalidation(nclass);

  if(confusion_opt[0]){
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
  const char* pszMessage;
  void* pProgressArg=NULL;
  GDALProgressFunc pfnProgress=GDALTermProgress;
  double progress=0;
  if(!verbose_opt[0])
    pfnProgress(progress,pszMessage,pProgressArg);
  if(reference_opt[0].find(".shp")!=string::npos){
    for(int iinput=0;iinput<input_opt.size();++iinput){
      if(output_opt.size())
        assert(reference_opt.size()==output_opt.size());
      for(int iref=0;iref<reference_opt.size();++iref){
        if(verbose_opt[0])
          cout << "reference is " << reference_opt[iref] << endl;
        assert(reference_opt[iref].find(".shp")!=string::npos);
        ImgReaderOgr referenceReader;
        try{
          inputReader.open(input_opt[iinput]);//,imagicX_opt[0],imagicY_opt[0]);
          if(mask_opt.size()){
            maskReader.open(mask_opt[iinput]);
            assert(inputReader.nrOfCol()==maskReader.nrOfCol());
            assert(inputReader.nrOfRow()==maskReader.nrOfRow());
          }
          referenceReader.open(reference_opt[iref]);
        }
        catch(string error){
          cerr << error << endl;
          exit(1);
        }
        if(confusion_opt[0])
          referenceRange=inputRange;

        ImgWriterOgr ogrWriter;
        OGRLayer *writeLayer;
        if(output_opt.size()){
          if(verbose_opt[0])
            cout << "creating output vector file " << output_opt[0] << endl;
          assert(output_opt[0].find(".shp")!=string::npos);
          try{
            ogrWriter.open(output_opt[iref]);
          }
          catch(string error){
            cerr << error << endl;
            exit(1);
          }
          char     **papszOptions=NULL;
          string layername=output_opt[0].substr(0,output_opt[0].find(".shp"));
          if(verbose_opt[0])
            cout << "creating layer: " << layername << endl;
          if(ogrWriter.createLayer(layername, "EPSG:3035", wkbPoint, papszOptions)==NULL)
            cout << "Error: create layer failed!" << endl;
          else if(verbose_opt[0])
            cout << "created layer" << endl;
          if(verbose_opt[0])
            cout << "copy fields from " << reference_opt[iref] << endl;
          ogrWriter.copyFields(referenceReader);
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
              ogrWriter.createField(fs.str(),OFTInteger);
            }
          }
          writeLayer=ogrWriter.getDataSource()->GetLayer(0);
        }
        OGRLayer  *readLayer;
        readLayer = referenceReader.getDataSource()->GetLayer(0);
        readLayer->ResetReading();
        OGRFeature *readFeature;
        int isample=0;
        while( (readFeature = readLayer->GetNextFeature()) != NULL ){
          if(verbose_opt[0])
            cout << "sample " << ++isample << endl;
          //get x and y from readFeature
          double x,y;
          OGRGeometry *poGeometry;
          OGRPolygon readPolygon;
          OGRPoint centroidPoint;
          OGRPoint *poPoint;
          poGeometry = readFeature->GetGeometryRef();
          // assert( poGeometry != NULL && wkbFlatten(poGeometry->getGeometryType()) == wkbPoint );
          if(poGeometry==NULL)
            continue;
          else if(wkbFlatten(poGeometry->getGeometryType()) == wkbPolygon){
            readPolygon = *((OGRPolygon *) poGeometry);
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
          short inputValue;
          vector<short> inputValues;
          bool isHomogeneous=true;
          short maskValue;
          short outputValue;
          //read referenceValue from feature
          short referenceValue;
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
          j_centre=static_cast<int>(j_centre);
          i_centre=static_cast<int>(i_centre);
          //check if j_centre is out of bounds
          if(static_cast<int>(j_centre)<0||static_cast<int>(j_centre)>=inputReader.nrOfRow())
            continue;
          //check if i_centre is out of bounds
          if(static_cast<int>(i_centre)<0||static_cast<int>(i_centre)>=inputReader.nrOfCol())
            continue;
          OGRFeature *writeFeature;
          if(output_opt.size()){
            writeFeature = OGRFeature::CreateFeature(writeLayer->GetLayerDefn());
            if(verbose_opt[0])
              cout << "copying fields from " << reference_opt[0] << endl;
            if(writeFeature->SetFrom(readFeature)!= OGRERR_NONE)
              cerr << "writing feature failed" << endl;
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
              if(static_cast<int>(j)<0||static_cast<int>(j)>=inputReader.nrOfRow())
                continue;
              int i=i_centre+windowI;
              //check if i is out of bounds
              if(static_cast<int>(i)<0||static_cast<int>(i)>=inputReader.nrOfCol())
                continue;
              if(verbose_opt[0])
                cout << setprecision(12) << "reading image value at x,y " << x << "," << y << " (" << i << "," << j << "), ";
              inputReader.readData(inputValue,GDT_Int16,i,j,band_opt[0]);
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
              maskFlagged=false;//(masknodata_opt[ivalue]>=0)?false:true;
              if(mask_opt.size()){
                maskReader.readData(maskValue,GDT_Int16,i,j,band_opt[0]);
                for(int ivalue=0;ivalue<masknodata_opt.size();++ivalue){
                  if(masknodata_opt[ivalue]>=0){//values set in masknodata_opt are invalid
                    if(maskValue==masknodata_opt[ivalue]){
                      maskFlagged=true;
                      break;
                    }
                  }
                  else{//only values set in masknodata_opt are valid
                    if(maskValue!=-masknodata_opt[ivalue])
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
            int j=j_centre;
            int i=i_centre;
            //flag if not all pixels are homogeneous or if at least one pixel flagged
          
            if(!windowHasFlag&&isHomogeneous){
              if(output_opt.size())
                writeFeature->SetField(labelclass_opt[0].c_str(),static_cast<int>(inputValue));
              if(confusion_opt[0]){
                ++ntotalValidation;
                if(classValueMap.size()){
                  assert(inputValue<nameVector.size());
                  string className=nameVector[inputValue];
                  cm.incrementResult(type2string<short>(classValueMap[referenceClassName]),type2string<short>(classValueMap[className]),1);
                }
                else{
                  int rc=distance(referenceRange.begin(),find(referenceRange.begin(),referenceRange.end(),referenceValue));
                  int ic=distance(inputRange.begin(),find(inputRange.begin(),inputRange.end(),inputValue));
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
                if(static_cast<int>(j)<0||static_cast<int>(j)>=inputReader.nrOfRow())
                  continue;
                int i=i_centre+windowI;
                //check if i is out of bounds
                if(static_cast<int>(i)<0||static_cast<int>(i)>=inputReader.nrOfCol())
                  continue;
                if(!windowAllFlagged){
                  ostringstream fs;
                  if(theDim>1)
                    fs << labelclass_opt[0] << "_" << windowJ << "_" << windowI;
                  else
                    fs << labelclass_opt[0];
                  if(output_opt.size())
                    writeFeature->SetField(fs.str().c_str(),static_cast<int>(inputValue));
                  if(!windowJ&&!windowI){//centre pixel
                    if(confusion_opt[0]){
                      ++ntotalValidation;
                      if(classValueMap.size()){
                        assert(inputValue<nameVector.size());
                        string className=nameVector[inputValue];
                        cm.incrementResult(type2string<short>(classValueMap[referenceClassName]),type2string<short>(classValueMap[className]),1);
                      }
                      else{
                        int rc=distance(referenceRange.begin(),find(referenceRange.begin(),referenceRange.end(),referenceValue));
                        int ic=distance(inputRange.begin(),find(inputRange.begin(),inputRange.end(),inputValue));
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
                string errorString="Failed to create feature in shapefile";
                throw(errorString);
              }
            }
            OGRFeature::DestroyFeature( writeFeature );
          }
        }
        if(output_opt.size())
          ogrWriter.close();
        referenceReader.close();
        inputReader.close();
        if(mask_opt.size())
          maskReader.close();
      }
    }
  }//reference is shape file
  else{
    ImgWriterGdal imgWriter;
    try{
      inputReader.open(input_opt[0]);//,imagicX_opt[0],imagicY_opt[0]);
      if(mask_opt.size())
        maskReader.open(mask_opt[0]);
      if(output_opt.size()){
        if(verbose_opt[0])
          cout << "opening output image " << output_opt[0] << endl;
        if(option_opt.findSubstring("INTERLEAVE=")==option_opt.end()){
          string theInterleave="INTERLEAVE=";
          theInterleave+=inputReader.getInterleave();
          option_opt.push_back(theInterleave);
        }
        imgWriter.open(output_opt[0],inputReader.nrOfCol(),inputReader.nrOfRow(),1,inputReader.getDataType(),inputReader.getImageType(),option_opt);
	if(nodata_opt.size())
	  imgWriter.GDALSetNoDataValue(nodata_opt[0]);
        if(inputReader.isGeoRef()){
          imgWriter.copyGeoTransform(inputReader);
        }
        if(colorTable_opt.size())
          imgWriter.setColorTable(colorTable_opt[0]);
        else if(inputReader.getColorTable()!=NULL){
          if(verbose_opt[0])
            cout << "set colortable from input image" << endl;
          imgWriter.setColorTable(inputReader.getColorTable());
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
    vector<short> lineInput(inputReader.nrOfCol());
    vector<short> lineMask(maskReader.nrOfCol());
    vector<short> lineOutput;
    if(output_opt.size())
      lineOutput.resize(inputReader.nrOfCol());

    int irow=0;
    int icol=0;
    double oldreferencerow=-1;
    ImgReaderGdal referenceReader;
    try{
      referenceReader.open(reference_opt[0]);//,rmagicX_opt[0],rmagicY_opt[0]);
    }
    catch(string error){
      cerr << error << endl;
      exit(1);
    }
    if(inputReader.isGeoRef()){
      assert(referenceReader.isGeoRef());
      if(inputReader.getProjection()!=referenceReader.getProjection())
        cout << "projection of input image and reference image are different!" << endl;
    }
    vector<short> lineReference(referenceReader.nrOfCol());
    if(confusion_opt[0]){
      referenceReader.getRange(referenceRange,band_opt[0]);
      for(int iflag=0;iflag<nodata_opt.size();++iflag){
        vector<short>::iterator fit;
        fit=find(referenceRange.begin(),referenceRange.end(),nodata_opt[iflag]);
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
    for(irow=0;irow<inputReader.nrOfRow()&&!isDifferent;++irow){
      //read line in lineInput, lineReference and lineMask
      inputReader.readData(lineInput,GDT_Int16,irow,band_opt[0]);
      if(mask_opt.size())
        maskReader.readData(lineMask,GDT_Int16,irow);
      double x,y;//geo coordinates
      double ireference,jreference;//image coordinates in reference image
      for(icol=0;icol<inputReader.nrOfCol();++icol){
        //find col in reference
        inputReader.image2geo(icol,irow,x,y);
        referenceReader.geo2image(x,y,ireference,jreference);
        if(ireference<0||ireference>=referenceReader.nrOfCol()){
          cerr << ireference << " out of reference range!" << endl;
          cerr << x << " " << y << " " << icol << " " << irow << endl;
          cerr << x << " " << y << " " << ireference << " " << jreference << endl;
          exit(1);
        }
        if(jreference!=oldreferencerow){
          if(jreference<0||jreference>=referenceReader.nrOfRow()){
            cerr << jreference << " out of reference range!" << endl;
            cerr << x << " " << y << " " << icol << " " << irow << endl;
            cerr << x << " " << y << " " << ireference << " " << jreference << endl;
            exit(1);
          }
          else{
            referenceReader.readData(lineReference,GDT_Int16,static_cast<int>(jreference),band_opt[0]);
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
          for(int ivalue=0;ivalue<masknodata_opt.size();++ivalue){
            if(lineMask[icol]==masknodata_opt[ivalue]){
              flagged=true;
              break;
            }
          }
        }
        if(!flagged){
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
            if(output_opt.empty()&&!confusion_opt[0]){
              isDifferent=true;
              break;
            }
            if(output_opt.size()){
              // if(lineInput[icol]<20){//forest
              //   if(lineReference[icol]>=20)//gain
              //     lineOutput[icol]=lineInput[icol]*10+1;//GAIN is 111,121,131
              //   else//forest type changed: mixed
              //     lineOutput[icol]=130;//MIXED FOREST
              // }
              // else if(lineReference[icol]<20){//loss
              //   lineOutput[icol]=20*10+lineReference[icol];//LOSS is 211 212 213
              // }
              // else//no forest
              //   lineOutput[icol]=20*10;//NON FOREST is 200
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
          imgWriter.writeData(lineOutput,GDT_Int16,irow);
        }
        catch(string errorstring){
          cerr << "lineOutput.size(): " << lineOutput.size() << endl;
          cerr << "imgWriter.nrOfCol(): " << imgWriter.nrOfCol() << endl;
          cerr << errorstring << endl;
          exit(1);
        }
      }
      else if(isDifferent&&!confusion_opt[0]){//we can break off here, files are different...
        if(!verbose_opt[0])
          pfnProgress(1.0,pszMessage,pProgressArg);
        break;
      }
      progress=static_cast<float>(irow+1.0)/inputReader.nrOfRow();
      if(!verbose_opt[0])
        pfnProgress(progress,pszMessage,pProgressArg);
    }
    if(output_opt.size())
      imgWriter.close();
    else if(!confusion_opt[0]){
      if(isDifferent)
        cout << input_opt[0] << " and " << reference_opt[0] << " are different" << endl;
      else
        cout << input_opt[0] << " and " << reference_opt[0] << " are identical" << endl;
    }
    referenceReader.close();
    inputReader.close();
    if(mask_opt.size())
      maskReader.close();
  }

  if(confusion_opt[0]){
    // double totalResult=0;
    // cout << " ";
    // for(int ic=0;ic<inputRange.size();++ic)
    //   cout << inputRange[ic] << " ";
    // cout << endl;
    // unsigned int ntotal=0;
    // vector<unsigned int> ntotalclass(referenceRange.size());
    // for(int rc=0;rc<referenceRange.size();++rc){
    //   ntotalclass[rc]=0;
    //   cout << referenceRange[rc] << " ";
    //   //initialize
    //   for(int ic=0;ic<inputRange.size();++ic){
    //     unsigned int result=0;
    //     user[ic]=0;
    //     producer[ic]=0;	
    //     for(int k=0;k<nclass;++k){
    //       user[ic]+=resultClass[k][ic];
    //       producer[ic]+=resultClass[ic][k];
    //     }	  
    //     result=resultClass[rc][ic];
    //     ntotal+=result;
    //     ntotalclass[rc]+=result;
    //     if(ic==rc){
    //       totalResult+=result;
    //     }
    //     cout << result << " ";
    //   }
    //   cout << endl;
    // }
    // if(verbose_opt[0]){
    //   cout << "totalResult: " << totalResult << endl;
    //   cout << "ntotalValidation: " << ntotalValidation << endl;
    //   cout << "nflagged: " << nflagged << endl;
    //   cout << "ntotal: " << ntotal << endl;
    // }
    // totalResult*=100.0/ntotal;
    // int nclass0=0;//number of classes without any reference
    // for(int rc=0;rc<referenceRange.size();++rc)
    //   if(!nvalidation[rc])
    //     ++nclass0;
    // double pChance=0;
    // double pCorrect=0;
    // double totalEntries=0;
    // cout << "class #samples userAcc prodAcc" << endl;    
    // for(int rc=0;rc<referenceRange.size();++rc){
    //   totalEntries+=user[rc];
    //   pChance+=user[rc]*producer[rc];
    //   pCorrect+=resultClass[rc][rc];
    //   cout << referenceRange[rc] << " " << nvalidation[rc] << " ";
    //   cout << 100.0*resultClass[rc][rc]/user[rc] << " "//user accuracy
    //        << 100.0*resultClass[rc][rc]/producer[rc] << " " << endl;//producer accuracy
    // }
    // pCorrect/=totalEntries;
    // pChance/=totalEntries*totalEntries;    
    // double kappa=pCorrect-pChance;
    // kappa/=1-pChance;
    // cout << "Kappa: " << kappa << endl;    
    // cout << "total weighted: " << static_cast<int>(0.5+totalResult) << endl;
    assert(cm.nReference());
    cout << cm << endl;
    cout << "class #samples userAcc prodAcc" << endl;
    double se95_ua=0;
    double se95_pa=0;
    double se95_oa=0;
    double dua=0;
    double dpa=0;
    double doa=0;
    for(int iclass=0;iclass<cm.nClasses();++iclass){
      dua=cm.ua_pct(classNames[iclass],&se95_ua);
      dpa=cm.pa_pct(classNames[iclass],&se95_pa);
      cout << cm.getClass(iclass) << " " << cm.nReference(cm.getClass(iclass)) << " " << dua << " (" << se95_ua << ")" << " " << dpa << " (" << se95_pa << ")" << endl;
    }
    doa=cm.oa(&se95_oa);
    cout << "Kappa: " << cm.kappa() << endl;
    cout << "Overall Accuracy: " << 100*doa << " (" << 100*se95_oa << ")"  << endl;
  }
}

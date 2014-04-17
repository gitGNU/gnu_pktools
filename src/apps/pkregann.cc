/**********************************************************************
pkregann.cc: regression with artificial neural network (multi-layer perceptron)
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
#include <vector>
#include <fstream>
#include "base/Optionpk.h"
#include "fileclasses/FileReaderAscii.h"
#include "floatfann.h"
#include "algorithms/myfann_cpp.h"
using namespace std;

int main(int argc, char *argv[])
{
  //--------------------------- command line options ------------------------------------
  Optionpk<string> input_opt("i", "input", "input ASCII file"); 
  Optionpk<string> output_opt("o", "output", "output ASCII file for result"); 
  Optionpk<int> inputCols_opt("ic", "inputCols", "input columns (e.g., for three dimensional input data in first three columns use: -ic 0 -ic 1 -ic 2"); 
  Optionpk<int> outputCols_opt("oc", "outputCols", "output columns (e.g., for two dimensional output in columns 3 and 4 (starting from 0) use: -oc 3 -oc 4"); 
  Optionpk<string> training_opt("t", "training", "training ASCII file (each row represents one sampling unit. Input features should be provided as columns, followed by output)"); 
  Optionpk<double> from_opt("from", "from", "start from this row in training file (start from 0)",0); 
  Optionpk<double> to_opt("to", "to", "read until this row in training file (start from 0 or set leave 0 as default to read until end of file)", 0); 
  Optionpk<short> band_opt("b", "band", "band index (starting from 0, either use band option or use start to end)");
  Optionpk<double> offset_opt("\0", "offset", "offset value for each spectral band input features: refl[band]=(DN[band]-offset[band])/scale[band]", 0.0);
  Optionpk<double> scale_opt("\0", "scale", "scale value for each spectral band input features: refl=(DN[band]-offset[band])/scale[band] (use 0 if scale min and max in each band to -1.0 and 1.0)", 0.0);
  Optionpk<unsigned short> cv_opt("cv", "cv", "n-fold cross validation mode",0);
  Optionpk<unsigned int> nneuron_opt("\0", "nneuron", "number of neurons in hidden layers in neural network (multiple hidden layers are set by defining multiple number of neurons: -n 15 -n 1, default is one hidden layer with 5 neurons)", 5); 
  Optionpk<float> connection_opt("\0", "connection", "connection reate (default: 1.0 for a fully connected network)", 1.0); 
  Optionpk<float> weights_opt("w", "weights", "weights for neural network. Apply to fully connected network only, starting from first input neuron to last output neuron, including the bias neurons (last neuron in each but last layer)", 0.0); 
  Optionpk<float> learning_opt("l", "learning", "learning rate (default: 0.7)", 0.7); 
  Optionpk<unsigned int> maxit_opt("\0", "maxit", "number of maximum iterations (epoch) (default: 500)", 500); 
  Optionpk<short> verbose_opt("v", "verbose", "set to: 0 (results only), 1 (confusion matrix), 2 (debug)",0);

  bool doProcess;//stop process when program was invoked with help option (-h --help)
  try{
    doProcess=input_opt.retrieveOption(argc,argv);
    output_opt.retrieveOption(argc,argv);
    inputCols_opt.retrieveOption(argc,argv);
    outputCols_opt.retrieveOption(argc,argv);
    training_opt.retrieveOption(argc,argv);
    from_opt.retrieveOption(argc,argv);
    to_opt.retrieveOption(argc,argv);
    band_opt.retrieveOption(argc,argv);
    offset_opt.retrieveOption(argc,argv);
    scale_opt.retrieveOption(argc,argv);
    cv_opt.retrieveOption(argc,argv);
    nneuron_opt.retrieveOption(argc,argv);
    connection_opt.retrieveOption(argc,argv);
    weights_opt.retrieveOption(argc,argv);
    learning_opt.retrieveOption(argc,argv);
    maxit_opt.retrieveOption(argc,argv);
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

  unsigned int ninput=inputCols_opt.size();
  unsigned int noutput=outputCols_opt.size();
  assert(ninput);
  assert(noutput);
  vector< vector<float> > inputUnits;
  vector< vector<float> > trainingUnits;
  vector< vector<float> > trainingOutput;
  FileReaderAscii inputFile;
  unsigned int inputSize=0;
  if(input_opt.size()){
    inputFile.open(input_opt[0]);
    inputFile.setMinRow(from_opt[0]);
    inputFile.setMaxRow(to_opt[0]);
    inputFile.setComment('#');
    inputFile.readData(inputUnits,inputCols_opt,1,0,true,verbose_opt[0]);
    inputFile.close();
    inputSize=inputUnits.size();
  }
  FileReaderAscii trainingFile(training_opt[0]);
  unsigned int sampleSize=0;
  trainingFile.setMinRow(from_opt[0]);
  trainingFile.setMaxRow(to_opt[0]);
  trainingFile.setComment('#');
  trainingFile.readData(trainingUnits,inputCols_opt,1,0,true,verbose_opt[0]);
  trainingFile.readData(trainingOutput,outputCols_opt,1,0,true,verbose_opt[0]);
  trainingFile.close();
  sampleSize=trainingUnits.size();

  if(verbose_opt[0]>1){
    std::cout << "sampleSize: " << sampleSize << std::endl;
    std::cout << "ninput: " << ninput << std::endl;
    std::cout << "noutput: " << noutput << std::endl;
    std::cout << "trainingUnits[0].size(): " << trainingUnits[0].size() << std::endl;
    std::cout << "trainingOutput[0].size(): " << trainingOutput[0].size() << std::endl;
    std::cout << "trainingUnits.size(): " << trainingUnits.size() << std::endl;
    std::cout << "trainingOutput.size(): " << trainingOutput.size() << std::endl;
  }

  assert(ninput==trainingUnits[0].size());
  assert(noutput==trainingOutput[0].size());
  assert(trainingUnits.size()==trainingOutput.size());

  //set scale and offset
  if(offset_opt.size()>1)
    assert(offset_opt.size()==ninput);
  if(scale_opt.size()>1)
    assert(scale_opt.size()==ninput);

  std::vector<float> offset_input(ninput);
  std::vector<float> scale_input(ninput);

  std::vector<float> offset_output(noutput);
  std::vector<float> scale_output(noutput);

  for(int iinput=0;iinput<ninput;++iinput){
    if(verbose_opt[0]>=1)
      cout << "scaling for input feature" << iinput << endl;
    offset_input[iinput]=(offset_opt.size()==1)?offset_opt[0]:offset_opt[iinput];
    scale_input[iinput]=(scale_opt.size()==1)?scale_opt[0]:scale_opt[iinput];
    //search for min and maximum
    if(scale_input[iinput]<=0){
      float theMin=trainingUnits[0][iinput];
      float theMax=trainingUnits[0][iinput];
      for(int isample=0;isample<trainingUnits.size();++isample){
        if(theMin>trainingUnits[isample][iinput])
          theMin=trainingUnits[isample][iinput];
        if(theMax<trainingUnits[isample][iinput])
          theMax=trainingUnits[isample][iinput];
      }
      offset_input[iinput]=theMin+(theMax-theMin)/2.0;
      scale_input[iinput]=(theMax-theMin)/2.0;
      if(verbose_opt[0]>=1){
        std::cout << "Extreme image values for input feature " << iinput << ": [" << theMin << "," << theMax << "]" << std::endl;
        std::cout << "Using offset, scale: " << offset_input[iinput] << ", " << scale_input[iinput] << std::endl;
        std::cout << "scaled values for input feature " << iinput << ": [" << (theMin-offset_input[iinput])/scale_input[iinput] << "," << (theMax-offset_input[iinput])/scale_input[iinput] << "]" << std::endl;
      }
    }
  }

  for(int ioutput=0;ioutput<noutput;++ioutput){
    if(verbose_opt[0]>=1)
      cout << "scaling for output feature" << ioutput << endl;
    //search for min and maximum
    float theMin=trainingOutput[0][ioutput];
    float theMax=trainingOutput[0][ioutput];
    for(int isample=0;isample<trainingOutput.size();++isample){
      if(theMin>trainingOutput[isample][ioutput])
        theMin=trainingOutput[isample][ioutput];
      if(theMax<trainingOutput[isample][ioutput])
        theMax=trainingOutput[isample][ioutput];
    }
    offset_output[ioutput]=theMin+(theMax-theMin)/2.0;
    scale_output[ioutput]=(theMax-theMin)/2.0;
    if(verbose_opt[0]>=1){
      std::cout << "Extreme image values for output feature " << ioutput << ": [" << theMin << "," << theMax << "]" << std::endl;
      std::cout << "Using offset, scale: " << offset_output[ioutput] << ", " << scale_output[ioutput] << std::endl;
      std::cout << "scaled values for output feature " << ioutput << ": [" << (theMin-offset_output[ioutput])/scale_output[ioutput] << "," << (theMax-offset_output[ioutput])/scale_output[ioutput] << "]" << std::endl;
    }
  }



  FANN::neural_net net;//the neural network

  const unsigned int num_layers = nneuron_opt.size()+2;
  const float desired_error = 0.0003;
  const unsigned int iterations_between_reports = (verbose_opt[0])? maxit_opt[0]+1:0;
  if(verbose_opt[0]>=1){
    cout << "creating artificial neural network with " << nneuron_opt.size() << " hidden layer, having " << endl;
    for(int ilayer=0;ilayer<nneuron_opt.size();++ilayer)
      cout << nneuron_opt[ilayer] << " ";
    cout << "neurons" << endl;
  }

  switch(num_layers){
  case(3):{
      unsigned int layers[3];
      layers[0]=ninput;
      layers[1]=nneuron_opt[0];
      layers[2]=noutput;
      net.create_sparse_array(connection_opt[0],num_layers,layers);
    // net.create_sparse(connection_opt[0],num_layers, ninput, nneuron_opt[0], noutput);
    break;
  }
  case(4):{
      unsigned int layers[3];
      layers[0]=ninput;
      layers[1]=nneuron_opt[0];
      layers[2]=nneuron_opt[1];
      layers[3]=noutput;
      net.create_sparse_array(connection_opt[0],num_layers,layers);
    // net.create_sparse(connection_opt[0],num_layers, ninput, nneuron_opt[0], nneuron_opt[1], noutput);
    break;
  }
  default:
    cerr << "Only 1 or 2 hidden layers are supported!" << endl;
    exit(1);
    break;
  }
  if(verbose_opt[0]>=1)
    cout << "network created" << endl;
  
  net.set_learning_rate(learning_opt[0]);

  //   net.set_activation_steepness_hidden(1.0);
  //   net.set_activation_steepness_output(1.0);
    
  net.set_activation_function_hidden(FANN::SIGMOID_SYMMETRIC_STEPWISE);
  net.set_activation_function_output(FANN::SIGMOID_SYMMETRIC_STEPWISE);

  // Set additional properties such as the training algorithm
  //   net.set_training_algorithm(FANN::TRAIN_QUICKPROP);

  // Output network type and parameters
  if(verbose_opt[0]>=1){
    cout << endl << "Network Type                         :  ";
    switch (net.get_network_type())
      {
      case FANN::LAYER:
        cout << "LAYER" << endl;
        break;
      case FANN::SHORTCUT:
        cout << "SHORTCUT" << endl;
        break;
      default:
        cout << "UNKNOWN" << endl;
        break;
      }
    net.print_parameters();
  }
      
  if(verbose_opt[0]>=1){
    cout << "Max Epochs " << setw(8) << maxit_opt[0] << ". "
         << "Desired Error: " << left << desired_error << right << endl;
  }
  bool initWeights=true;

  Vector2d<float> trainingFeatures(sampleSize,ninput);
  for(unsigned int isample=0;isample<sampleSize;++isample){
    for(unsigned int iinput=0;iinput<ninput;++iinput)
      trainingFeatures[isample][iinput]=(trainingUnits[isample][iinput]-offset_input[iinput])/scale_input[iinput];
  }

  Vector2d<float> scaledOutput(sampleSize,noutput);
  for(unsigned int isample=0;isample<sampleSize;++isample){
    for(unsigned int ioutput=0;ioutput<noutput;++ioutput)
      scaledOutput[isample][ioutput]=(trainingOutput[isample][ioutput]-offset_output[ioutput])/scale_output[ioutput];
  }

  if(cv_opt[0]){
    if(verbose_opt[0])
      std::cout << "cross validation" << std::endl;
    std::vector< std::vector<float> > referenceVector;
    std::vector< std::vector<float> > outputVector;
    net.cross_validation(trainingFeatures,
                         scaledOutput,
                         cv_opt[0],
                         maxit_opt[0],
                         desired_error,
                         referenceVector,
                         outputVector);
    assert(referenceVector.size()==outputVector.size());
    vector<double> rmse(noutput);
    for(int isample=0;isample<referenceVector.size();++isample){
      std::cout << isample << " ";
      for(int ioutput=0;ioutput<noutput;++ioutput){
        if(!isample)
          rmse[ioutput]=0;
        double ref=scale_output[ioutput]*referenceVector[isample][ioutput]+offset_output[ioutput];
        double val=scale_output[ioutput]*outputVector[isample][ioutput]+offset_output[ioutput];
        rmse[ioutput]+=(ref-val)*(ref-val);
        std::cout << ref << " " << val;
        if(ioutput<noutput-1)
          std::cout << " ";
        else
          std::cout << std::endl;
      }
    }
    for(int ioutput=0;ioutput<noutput;++ioutput)
      std::cout << "rmse output variable " << ioutput << ": " << sqrt(rmse[ioutput]/referenceVector.size()) << std::endl;
  }


  net.train_on_data(trainingFeatures,
                    scaledOutput,
                    initWeights,
                    maxit_opt[0],
                    iterations_between_reports,
                    desired_error);


  if(verbose_opt[0]>=2){
    net.print_connections();
    vector<fann_connection> convector;
    net.get_connection_array(convector);
    for(int i_connection=0;i_connection<net.get_total_connections();++i_connection)
      cout << "connection " << i_connection << ": " << convector[i_connection].weight << endl;
  }
  //end of training

  ofstream outputStream;
  if(!output_opt.empty())
    outputStream.open(output_opt[0].c_str(),ios::out);
  for(unsigned int isample=0;isample<inputUnits.size();++isample){
    std::vector<float> inputFeatures(ninput);
    for(unsigned int iinput=0;iinput<ninput;++iinput)
      inputFeatures[iinput]=(inputUnits[isample][iinput]-offset_input[iinput])/scale_input[iinput];
    vector<float> result(noutput);
    result=net.run(inputFeatures);

    if(!output_opt.empty())
      outputStream << isample << " ";
    else
      std::cout << isample << " ";
    if(verbose_opt[0]){
      for(unsigned int iinput=0;iinput<ninput;++iinput){
        if(output_opt.size())
          outputStream << inputUnits[isample][iinput] << " ";
        else
          std::cout << inputUnits[isample][iinput] << " ";
      }
    }
    for(unsigned int ioutput=0;ioutput<noutput;++ioutput){
      result[ioutput]=scale_output[ioutput]*result[ioutput]+offset_output[ioutput];
      if(output_opt.size()){
        outputStream << result[ioutput];
        if(ioutput<noutput-1)
          outputStream << " ";
        else
          outputStream << std::endl;
      }
      else{
        std::cout << result[ioutput];
        if(ioutput<noutput-1)
          std::cout << " ";
        else
          std::cout << std::endl;
      }
    }
  }
  if(!output_opt.empty())
    outputStream.close();
}

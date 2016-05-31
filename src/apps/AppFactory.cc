#include <algorithm>
#include <iostream>
#include <string>
#include "imageclasses/ImgReaderOgr.h"
#include "base/Vector2d.h"
#include "base/Optionpk.h"
#include "algorithms/StatFactory.h"
#include "algorithms/Egcs.h"
#include "AppFactory.h"

using namespace std;
using namespace app;

void AppFactory::setOptions(int argc, char* argv[]){
  m_argc=argc;
  m_argv.clear();
  for(int iarg=0;iarg<argc;++iarg)
    m_argv.push_back(argv[iarg]);
}

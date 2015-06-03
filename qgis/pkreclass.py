# -*- coding: utf-8 -*-

"""
***************************************************************************
    pkreclass.py
    ---------------------
    Date                 : April 2015
    Copyright            : (C) 2015 by Pieter Kempeneers
    Email                : kempenep at gmail dot com
***************************************************************************
*                                                                         *
*   This program is free software; you can redistribute it and/or modify  *
*   it under the terms of the GNU General Public License as published by  *
*   the Free Software Foundation; either version 2 of the License, or     *
*   (at your option) any later version.                                   *
*                                                                         *
***************************************************************************
"""

__author__ = 'Pieter Kempeneers'
__date__ = 'April 2015'
__copyright__ = '(C) 2015, Pieter Kempeneers'
# This will get replaced with a git SHA1 when you do a git archive
__revision__ = '$Format:%H$'

import os
from pktoolsUtils import pktoolsUtils
from pktoolsAlgorithm import pktoolsAlgorithm
from processing.core.parameters import ParameterMultipleInput
from processing.core.parameters import ParameterRaster
from processing.core.outputs import OutputRaster
from processing.core.parameters import ParameterSelection
from processing.core.parameters import ParameterNumber
from processing.core.parameters import ParameterString
from processing.core.parameters import ParameterBoolean
from processing.core.parameters import ParameterExtent

class pkreclass(pktoolsAlgorithm):

    INPUT = "INPUT"
    OUTPUT = "OUTPUT"
    CLASS = "CLASS"
    BAND = "BAND"
    RECLASS = "RECLASS"
    MASK = "MASK"
    MSKNODATA = "MSKNODATA"
    NODATA = "NODATA"
    RTYPE = 'RTYPE'
    TYPE = ['none', 'Byte','Int16','UInt16','UInt32','Int32','Float32','Float64','CInt16','CInt32','CFloat32','CFloat64']
    EXTRA = 'EXTRA'

    def cliName(self):
        return "pkreclass"

    def defineCharacteristics(self):
        self.name = "reclass raster datasets"
        self.group = "[pktools] raster"
        self.addParameter(ParameterMultipleInput(self.INPUT, 'Input layer raster data set',ParameterMultipleInput.TYPE_RASTER))
        self.addOutput(OutputRaster(self.OUTPUT, "Output raster data set"))
        self.addParameter(ParameterString(self.BAND, "Band index(es) to replace, e.g., 0;1;2 (other bands are copied to output)", '0'))
        self.addParameter(ParameterRaster(self.MASK, "Mask raster dataset",optional=True))
        self.addParameter(ParameterString(self.MSKNODATA, "Mask value(s) not to consider for classification (e.g., 0;255)","0"))
        self.addParameter(ParameterString(self.CLASS, "list of classes to reclass, in combination with reclass option, e.g., 0;1;2;3",""))
        self.addParameter(ParameterString(self.RECLASS, "list of recoded classes, in combination with class option e.g., 10;11;12;13",""))
        self.addParameter(ParameterNumber(self.NODATA, "nodata value to put in image if not valid",0,None,0))

        self.addParameter(ParameterSelection(self.RTYPE, 'Output raster type (leave as none to keep original type)', self.TYPE, 0))
        self.addParameter(ParameterString(self.EXTRA,
                          'Additional parameters', '-of GTiff', optional=True))

    def processAlgorithm(self, progress):
        cliPath = "\"" + os.path.join(pktoolsUtils.pktoolsPath(), self.cliName()) + "\""
        commands = [cliPath]

        commands.append('-i')
        commands.append(self.getParameterValue(self.INPUT))

        if self.TYPE[self.getParameterValue(self.RTYPE)] != "none":
            commands.append('-ot')
            commands.append(self.TYPE[self.getParameterValue(self.RTYPE)])

        commands.append("-o")
        commands.append(self.getOutputValue(self.OUTPUT))

        commands.append('-nodata')
        commands.append(str(self.getParameterValue(self.NODATA)))
        
        band=str(self.getParameterValue(self.BAND))
        if band != '':
            bandValues = band.split(';')
            for bandValue in bandValues:
                commands.append('-b')
                commands.append(bandValue)

        theclass=str(self.getParameterValue(self.CLASS))
        if theclass != '':
            classValues = theclass.split(';')
            for classValue in classValues:
                commands.append('-c')
                commands.append(classValue)
        reclass=str(self.getParameterValue(self.RECLASS))
        if reclass != '':
            reclassValues = reclass.split(';')
            for reclassValue in reclassValues:
                commands.append('-r')
                commands.append(reclassValue)

        mask = str(self.getParameterValue(self.MASK))
        if mask != "None":
            commands.append('-m')
            commands.append(mask)
            msknodata=str(self.getParameterValue(self.MSKNODATA))
            msknodataValues = msknodata.split(';')
            for msknodataValue in msknodataValues:
                commands.append('-msknodata')
                commands.append(msknodataValue)

        extra = str(self.getParameterValue(self.EXTRA))
        if len(extra) > 0:
            commands.append(extra)

        pktoolsUtils.runpktools(commands, progress)

# -*- coding: utf-8 -*-

"""
***************************************************************************
    pkcrop.py
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

class pkcrop(pktoolsAlgorithm):

    INPUT = "INPUT"
    OUTPUT = "OUTPUT"
    DX = "DX"
    DY = "DY"
    PROJWIN = 'PROJWIN'
    BAND = "BAND"
    NODATA = "NODATA"
    RESAMPLE_OPTIONS = ['near', 'bilinear']
    RESAMPLE = "RESAMPLE"
    RTYPE = 'RTYPE'
    TYPE = ['none', 'Byte','Int16','UInt16','UInt32','Int32','Float32','Float64','CInt16','CInt32','CFloat32','CFloat64']
    EXTRA = 'EXTRA'

    def cliName(self):
        return "pkcrop"

    def defineCharacteristics(self):
        self.name = "crop raster datasets"
        self.group = "[pktools] raster"
        self.addParameter(ParameterMultipleInput(self.INPUT, 'Input layer raster data set',ParameterMultipleInput.TYPE_RASTER))
        self.addOutput(OutputRaster(self.OUTPUT, "Output raster data set"))
        self.addParameter(ParameterSelection(self.RTYPE, 'Output raster type (leave as none to keep original type)', self.TYPE, 0))
        self.addParameter(ParameterNumber(self.DX, "Output resolution in x (leave 0 for no change)",0.0,None,0.0))
        self.addParameter(ParameterNumber(self.DY, "Output resolution in y (leave 0 for no change)",0.0,None,0.0))
        self.addParameter(ParameterExtent(self.PROJWIN,
                          'Georeferenced boundingbox'))
        self.addParameter(ParameterString(self.NODATA, "invalid value(s) for input raster dataset (e.g., 0;255)","none"))
        self.addParameter(ParameterString(self.BAND, "Band(s) in input image to crop, e.g., 0;1;2 (leave empty to retain all bands)",'', optional=True))
        self.addParameter(ParameterSelection(self.RESAMPLE,"resampling method",self.RESAMPLE_OPTIONS, 0))
        self.addParameter(ParameterString(self.EXTRA,
                          'Additional parameters', '-of GTiff', optional=True))

    def processAlgorithm(self, progress):
        cliPath = '"' + os.path.join(pktoolsUtils.pktoolsPath(), self.cliName()) + '"'
        commands = [cliPath]

        input=self.getParameterValue(self.INPUT)
        inputFiles = input.split(';')
        for inputFile in inputFiles:
            commands.append('-i')
            commands.append('"' + inputFile + '"')

        if self.TYPE[self.getParameterValue(self.RTYPE)] != "none":
            commands.append('-ot')
            commands.append(self.TYPE[self.getParameterValue(self.RTYPE)])
        output=self.getParameterValue(self.OUTPUT)
        if output != "":
            commands.append("-o")
            commands.append('"' + output + '"')
        if self.getParameterValue(self.DX) != 0:
            commands.append("-dx")
            commands.append(str(self.getParameterValue(self.DX)))
        if self.getParameterValue(self.DY) != 0:
            commands.append("-dy")
            commands.append(str(self.getParameterValue(self.DY)))

        projwin = str(self.getParameterValue(self.PROJWIN))
        if(str(projwin).find(',')>0):
           regionCoords = projwin.split(',')
           commands.append('-ulx')
           commands.append(regionCoords[0])
           commands.append('-uly')
           commands.append(regionCoords[3])
           commands.append('-lrx')
           commands.append(regionCoords[1])
           commands.append('-lry')
           commands.append(regionCoords[2])

        nodata=self.getParameterValue(self.NODATA)
        if nodata != "none":
            nodataValues = nodata.split(';')
            for nodataValue in nodataValues:
                commands.append('-nodata')
                commands.append(nodataValue)
        
        band=self.getParameterValue(self.BAND)
        if band != '':
            bandValues = band.split(';')
            for bandValue in bandValues:
                commands.append('-b')
                commands.append(bandValue)
        commands.append("-r")
        commands.append(self.RESAMPLE_OPTIONS[self.getParameterValue(self.RESAMPLE)])

        extra = str(self.getParameterValue(self.EXTRA))
        if len(extra) > 0:
            commands.append(extra)

        pktoolsUtils.runpktools(commands, progress)

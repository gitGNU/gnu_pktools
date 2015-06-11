# -*- coding: utf-8 -*-

"""
***************************************************************************
    pklas2img.py
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
from processing.core.parameters import ParameterFile
from processing.core.outputs import OutputRaster
from processing.core.parameters import ParameterSelection
from processing.core.parameters import ParameterNumber
from processing.core.parameters import ParameterString
from processing.core.parameters import ParameterBoolean
from processing.core.parameters import ParameterExtent

class pklas2img(pktoolsAlgorithm):

    INPUT = "INPUT"
    OUTPUT = "OUTPUT"
    ATTRIBUTE_OPTIONS = ["z","intensity", "return", "nreturn"]
    COMPOSITE_OPTIONS = ["last", "min", "max", "median", "mean", "sum", "first", "profile" "percentile", "height", "values", "percentile", "number"]
    FILTER_OPTIONS = ["all","first","last","single","multiple"]

    ATTRIBUTE = "ATTRIBUTE"
    COMPOSITE = "COMPOSITE"
    FILTER = "FILTER"

    PERCENTILE = "PERCENTILE"
    DX = "DX"
    DY = "DY"
#    PROJWIN = 'PROJWIN'
    NODATA = "NODATA"
    RTYPE = 'RTYPE'
    TYPE = ['Float32','Byte','Int16','UInt16','UInt32','Int32','Float64','CInt16','CInt32','CFloat32','CFloat64']
    EXTRA = 'EXTRA'

    def cliName(self):
        return "pklas2img"

    def defineCharacteristics(self):
        self.name = "(Linux only) Create raster dataset from LAS(Z) data point cloud(s)"
        self.group = "[pktools] LiDAR"

        self.addParameter(ParameterFile(self.INPUT, "Input LAS(Z) data set(s)", False, False))

        self.addParameter(ParameterSelection(self.ATTRIBUTE,"name of the point attribute to select",self.ATTRIBUTE_OPTIONS, 0))
        self.addParameter(ParameterSelection(self.COMPOSITE,"composite for multiple points in cell",self.COMPOSITE_OPTIONS, 0))
        self.addParameter(ParameterSelection(self.FILTER,"filter las points",self.FILTER_OPTIONS, 0))
        self.addOutput(OutputRaster(self.OUTPUT, "Output raster data set"))
        self.addParameter(ParameterSelection(self.RTYPE, 'Output raster type', self.TYPE, 0))
        self.addParameter(ParameterNumber(self.PERCENTILE, "Percentile value used for rule percentile",0.0,100.0,95))
        self.addParameter(ParameterNumber(self.DX, "Output resolution in x",0.0,None,1.0))
        self.addParameter(ParameterNumber(self.DY, "Output resolution in y",0.0,None,1.0))
        # self.addParameter(ParameterExtent(self.PROJWIN,
        #                   'Georeferenced boundingbox'))
        self.addParameter(ParameterNumber(self.NODATA, "nodata value to put in image",0,None,0))
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

        output=self.getOutputValue(self.OUTPUT)
        if output != "":
            commands.append("-o")
            commands.append('"' + output + '"')

        commands.append("-n")
        commands.append(self.ATTRIBUTE_OPTIONS[self.getParameterValue(self.ATTRIBUTE)])
        commands.append("-comp")
        commands.append(self.COMPOSITE_OPTIONS[self.getParameterValue(self.COMPOSITE)])
        commands.append("-fir")
        commands.append(self.FILTER_OPTIONS[self.getParameterValue(self.FILTER)])
        if self.getParameterValue(self.DX) != 0:
            commands.append("-dx")
            commands.append(str(self.getParameterValue(self.DX)))
        if self.getParameterValue(self.DY) != 0:
            commands.append("-dy")
            commands.append(str(self.getParameterValue(self.DY)))
        #projwin = str(self.getParameterValue(self.PROJWIN))
        # regionCoords = projwin.split(',')
        # commands.append('-ulx')
        # commands.append(regionCoords[0])
        # commands.append('-uly')
        # commands.append(regionCoords[3])
        # commands.append('-lrx')
        # commands.append(regionCoords[1])
        # commands.append('-lry')
        # commands.append(regionCoords[2])
        commands.append('-nodata')
        commands.append(str(self.getParameterValue(self.NODATA)))

        extra = str(self.getParameterValue(self.EXTRA))
        if len(extra) > 0:
            commands.append(extra)

        pktoolsUtils.runpktools(commands, progress)

# -*- coding: utf-8 -*-

"""
***************************************************************************
    pkfilter_spectral.py
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

class pkfilter_spectral(pktoolsAlgorithm):

    INPUT = "INPUT"
    OUTPUT = "OUTPUT"
    METHOD_OPTIONS = ["none", "median", "var", "min", "max", "sum", "mean", "dilate", "erode", "close", "open", "smooth", "density", "smoothnodata  values", "threshold local filtering", "stdev", "dwt", "dwti", "dwt_cut", "dwt_cut_from", "savgolay", "percentile"]
    METHOD = "METHOD"
    DZ = "DZ"
    NODATA = "NODATA"
    PADDING_OPTIONS = ["symmetric", "replicate", "circular", "zero"]
    PADDING = "PADDING"
    RTYPE = 'RTYPE'
    TYPE = ['none', 'Byte','Int16','UInt16','UInt32','Int32','Float32','Float64','CInt16','CInt32','CFloat32','CFloat64']
    EXTRA = 'EXTRA'

    def cliName(self):
        return "pkfilter"

    def defineCharacteristics(self):
        self.name = "spectral/temporal filter"
        self.group = "[pktools] filter"
        self.addParameter(ParameterRaster(self.INPUT, 'Input layer raster data set',ParameterRaster))
        self.addParameter(ParameterSelection(self.METHOD,"filter rule",self.METHOD_OPTIONS, 0))
        self.addOutput(OutputRaster(self.OUTPUT, "Output raster data set"))
        self.addParameter(ParameterSelection(self.RTYPE, 'Output raster type (leave as none to keep original type)', self.TYPE, 0))
        self.addParameter(ParameterNumber(self.DZ, "Filter kernel size",0.0,None,1.0))
        #for smooth nodata:
        self.addParameter(ParameterString(self.NODATA, "nodata value to smooth(e.g., 0;255)","none"))
        self.addParameter(ParameterSelection(self.PADDING,"Padding (edge effects)",self.PADDING_OPTIONS, 0))
        self.addParameter(ParameterString(self.EXTRA,
                          'Additional parameters', '-of GTiff', optional=True))

    def processAlgorithm(self, progress):
        commands = [os.path.join(pktoolsUtils.pktoolsPath(), self.cliName())]
        input=self.getParameterValue(self.INPUT)
        if input != "":
            commands.append('-i')
            commands.append(input)
        method=self.METHOD_OPTIONS[self.getParameterValue(self.METHOD)]
        if method != "none":
            commands.append("-f")
            commands.append(method)
        commands.append("-pad")
        commands.append(self.PADDING_OPTIONS[self.getParameterValue(self.PADDING)])

        if self.TYPE[self.getParameterValue(self.RTYPE)] != "none":
            commands.append('-ot')
            commands.append(self.TYPE[self.getParameterValue(self.RTYPE)])
        output=self.getParameterValue(self.OUTPUT)
        if output != "":
            commands.append("-o")
            commands.append(self.getOutputValue(self.OUTPUT))
        if self.getParameterValue(self.DZ) != 0:
            commands.append("-dz")
            commands.append(str(self.getParameterValue(self.DZ)))
        nodata=self.getParameterValue(self.NODATA)
        if nodata != "none":
            nodataValues = nodata.split(';')
            for nodataValue in nodataValues:
                commands.append('-nodata')
                commands.append(nodataValue)
        extra = str(self.getParameterValue(self.EXTRA))
        if len(extra) > 0:
            commands.append(extra)

        pktoolsUtils.runpktools(commands, progress)

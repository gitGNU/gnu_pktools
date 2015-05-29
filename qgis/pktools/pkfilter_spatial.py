# -*- coding: utf-8 -*-

"""
***************************************************************************
    pkfilter_spatial.py
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

class pkfilter_spatial(pktoolsAlgorithm):

    INPUT = "INPUT"
    OUTPUT = "OUTPUT"
    METHOD_OPTIONS = ["none", "median", "var", "min", "max", "sum", "mean", "dilate", "erode", "close", "open", "homog ", "heterog ", "sobelx ", "sobely ", "sobelxy ", "sobelyx" , "smooth", "countid", "smoothnodata  values", "threshold local filtering", "ismin", "ismax", "order ", "stdev", "mrf", "dwt", "dwti", "dwt_cut", "dwt_cut_from", "scramble", "shift", "savgolay", "percentile"]
    METHOD = "METHOD"
#    RESAMPLE_OPTIONS = ['near', 'bilinear']
#    RESAMPLE = "RESAMPLE"
    DIM = "DIM"
    NODATA = "NODATA"
    PADDING_OPTIONS = ["symmetric", "replicate", "circular", "zero"]
    PADDING = "PADDING"
    RTYPE = 'RTYPE'
    TYPE = ['none', 'Byte','Int16','UInt16','UInt32','Int32','Float32','Float64','CInt16','CInt32','CFloat32','CFloat64']
    EXTRA = 'EXTRA'

    def cliName(self):
        return "pkfilter"

    def defineCharacteristics(self):
        self.name = "spatial filter"
        self.group = "[pktools] filter"

        self.addParameter(ParameterRaster(self.INPUT, 'Input layer raster data set',ParameterRaster))
        self.addParameter(ParameterSelection(self.METHOD,"filter rule",self.METHOD_OPTIONS, 0))
        self.addOutput(OutputRaster(self.OUTPUT, "Output raster data set"))
        self.addParameter(ParameterSelection(self.RTYPE, 'Output raster type (leave as none to keep original type)', self.TYPE, 0))
        self.addParameter(ParameterNumber(self.DIM, "Filter kernel size (odd value)",0.0,None,3.0))
        #for smooth nodata:
        self.addParameter(ParameterString(self.NODATA, "invalid value(s) for input raster dataset (e.g., 0;255)","none"))
        self.addParameter(ParameterSelection(self.PADDING,"Padding (edge effects)",self.PADDING_OPTIONS, 0))
#        self.addParameter(ParameterSelection(self.RESAMPLE,"resampling method",self.RESAMPLE_OPTIONS, 0))
        self.addParameter(ParameterString(self.EXTRA,
                          'Additional parameters', '', optional=True))

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
        if self.getParameterValue(self.DIM) != 0:
            commands.append("-dx")
            commands.append(str(self.getParameterValue(self.DIM)))
            commands.append("-dy")
            commands.append(str(self.getParameterValue(self.DIM)))
        nodata=self.getParameterValue(self.NODATA)
        if nodata != "none":
            nodataValues = nodata.split(';')
            for nodataValue in nodataValues:
                commands.append('-nodata')
                commands.append(nodataValue)
#        commands.append("-r")
#        commands.append(self.RESAMPLE_OPTIONS[self.getParameterValue(self.RESAMPLE)])
        extra = str(self.getParameterValue(self.EXTRA))
        if len(extra) > 0:
            commands.append(extra)

        pktoolsUtils.runpktools(commands, progress)

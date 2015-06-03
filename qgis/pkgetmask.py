# -*- coding: utf-8 -*-

"""
***************************************************************************
    pkgetmask.py
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

class pkgetmask(pktoolsAlgorithm):

    INPUT = "INPUT"
    BAND = "BAND"
    MIN = "MIN"
    MAX = "MAX"
    OPERATOR_OPTIONS = ["OR", "AND"]
    OPERATOR = "OPERATOR"
    DATA = "DATA"
    NODATA = "NODATA"
    OUTPUT = "OUTPUT"
    RTYPE = 'RTYPE'
    TYPE = ['none', 'Byte','Int16','UInt16','UInt32','Int32','Float32','Float64','CInt16','CInt32','CFloat32','CFloat64']
    EXTRA = 'EXTRA'

    def cliName(self):
        return "pkgetmask"

    def defineCharacteristics(self):
        self.name = "create mask from raster dataset"
        self.group = "[pktools] raster"
        self.addParameter(ParameterRaster(self.INPUT, 'Input layer raster data set',ParameterRaster))
        self.addParameter(ParameterString(self.BAND, "Band(s) used for mask (e.g., 0;1)","0"))
        self.addParameter(ParameterString(self.MIN, "Minimum valid value (one value per band)","none"))
        self.addParameter(ParameterString(self.MAX, "Maximum valid value (one value per band)","none"))
        self.addParameter(ParameterSelection(self.OPERATOR,"getmask rule",self.OPERATOR_OPTIONS, 0))
        self.addParameter(ParameterString(self.DATA, "write value(s) for valid pixels (e.g., 0;255)","1"))
        self.addParameter(ParameterString(self.NODATA, "write value(s) for invalid pixels","0"))
        self.addOutput(OutputRaster(self.OUTPUT, "Output raster data set"))
        self.addParameter(ParameterSelection(self.RTYPE, 'Output raster type (leave as none to keep original type)', self.TYPE, 0))
        self.addParameter(ParameterString(self.EXTRA,
                          'Additional parameters', '-of GTiff', optional=True))

    def processAlgorithm(self, progress):
        cliPath = "\"" + os.path.join(pktoolsUtils.pktoolsPath(), self.cliName()) + "\""
        commands = [cliPath]

        input=self.getParameterValue(self.INPUT)
        commands.append('-i')
        commands.append(input)

        band=self.getParameterValue(self.BAND)
        bandValues = band.split(';')
        for bandValue in bandValues:
            commands.append('-band')
            commands.append(bandValue)
        min=self.getParameterValue(self.MIN)
        if min != "none":
            minValues = min.split(';')
            for minValue in minValues:
                commands.append('-min')
                commands.append(minValue)
        max=self.getParameterValue(self.MAX)
        if max != "none":
            maxValues = max.split(';')
            for maxValue in maxValues:
                commands.append('-max')
                commands.append(maxValue)
        commands.append("-p")
        commands.append(self.OPERATOR_OPTIONS[self.getParameterValue(self.OPERATOR)])
        if data != "none":
            dataValues = data.split(';')
            for dataValue in dataValues:
                commands.append('-data')
                commands.append(dataValue)
        nodata=self.getParameterValue(self.NODATA)
        if nodata != "none":
            nodataValues = nodata.split(';')
            for nodataValue in nodataValues:
                commands.append('-nodata')
                commands.append(nodataValue)
        if self.TYPE[self.getParameterValue(self.RTYPE)] != "none":
            commands.append('-ot')
            commands.append(self.TYPE[self.getParameterValue(self.RTYPE)])
        output=self.getParameterValue(self.OUTPUT)
        if output != "":
            commands.append("-o")
            commands.append(self.getOutputValue(self.OUTPUT))
        data=self.getParameterValue(self.DATA)

        extra = str(self.getParameterValue(self.EXTRA))
        if len(extra) > 0:
            commands.append(extra)

        pktoolsUtils.runpktools(commands, progress)

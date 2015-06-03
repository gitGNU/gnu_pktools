# -*- coding: utf-8 -*-

"""
***************************************************************************
    pkfilterdem.py
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

class pkfilterdem(pktoolsAlgorithm):

    INPUT = "INPUT"
    OUTPUT = "OUTPUT"
    DIM = "DIM"
    RTYPE = 'RTYPE'
    TYPE = ['Float32','Byte','Int16','UInt16','UInt32','Int32','Float64','CInt16','CInt32','CFloat32','CFloat64']
    FILTER_OPTIONS = ["promorph"]
    FILTER = "FILTER"
    EXTRA = 'EXTRA'

    def cliName(self):
        return "pkfilterdem"

    def defineCharacteristics(self):
        self.name = "Create DTM from DEM raster dataset)"
        self.group = "[pktools] LiDAR"

        self.addParameter(ParameterRaster(self.INPUT, 'Input layer raster data set',ParameterRaster))
        self.addParameter(ParameterSelection(self.FILTER,"filter",self.FILTER_OPTIONS, 0))
        self.addParameter(ParameterNumber(self.DIM, "maximum filter kernel size",3,None,17))

        self.addOutput(OutputRaster(self.OUTPUT, "Output raster data set"))
        self.addParameter(ParameterSelection(self.RTYPE, 'Output raster type', self.TYPE, 0))
        self.addParameter(ParameterString(self.EXTRA,
                          'Additional parameters', '-of GTiff', optional=True))

    def processAlgorithm(self, progress):
        cliPath = '"' + os.path.join(pktoolsUtils.pktoolsPath(), self.cliName()) + '"'
        commands = [cliPath]

        input=self.getParameterValue(self.INPUT)
        if input != "":
            commands.append('-i')
            commands.append('"' + input + '"')

        filter=self.FILTER_OPTIONS[self.getParameterValue(self.FILTER)]
        if filter != "none":
            commands.append("-f")
            commands.append(filter)
        if self.getParameterValue(self.DIM) != 0:
            commands.append("-dim")
            commands.append(str(self.getParameterValue(self.DIM)))

        if self.TYPE[self.getParameterValue(self.RTYPE)] != "none":
            commands.append('-ot')
            commands.append(self.TYPE[self.getParameterValue(self.RTYPE)])
        output=self.getOutputValue(self.OUTPUT)
        if output != "":
            commands.append("-o")
            commands.append('"' + output + '"')

        extra = str(self.getParameterValue(self.EXTRA))
        if len(extra) > 0:
            commands.append(extra)

        pktoolsUtils.runpktools(commands, progress)

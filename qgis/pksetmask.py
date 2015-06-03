# -*- coding: utf-8 -*-

"""
***************************************************************************
    pksetmask.py
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

class pksetmask(pktoolsAlgorithm):

    INPUT = "INPUT"
    MASK = "MASK"
    MSKNODATA = "MSKNODATA"
    MSKBAND = "MSKBAND"
    OPERATOR_OPTIONS = ["=","<",">","!"]
    OPERATOR = "OPERATOR"
    NODATA = "NODATA"
    OUTPUT = "OUTPUT"
    RTYPE = 'RTYPE'
    TYPE = ['none', 'Byte','Int16','UInt16','UInt32','Int32','Float32','Float64','CInt16','CInt32','CFloat32','CFloat64']
    EXTRA = 'EXTRA'

    def cliName(self):
        return "pksetmask"

    def defineCharacteristics(self):
        self.name = "apply mask to raster dataset"
        self.group = "[pktools] raster"
        self.addParameter(ParameterRaster(self.INPUT, 'Input layer raster data set',ParameterRaster))
        self.addParameter(ParameterMultipleInput(self.MASK, 'Mask(s) to apply',ParameterMultipleInput.TYPE_RASTER))
        self.addParameter(ParameterString(self.MSKNODATA, "Mask value(s), provide value for each mask (e.g., 250;255)","1"))
        self.addParameter(ParameterString(self.MSKBAND, "Mask band(s) to read, provide band for each mask (e.g., 0;1)","0"))
        self.addParameter(ParameterSelection(self.OPERATOR,"setmask rule",self.OPERATOR_OPTIONS, 0))
        self.addParameter(ParameterString(self.NODATA, "nodata value to put in image if not valid","0"))
        self.addOutput(OutputRaster(self.OUTPUT, "Output raster data set"))
        self.addParameter(ParameterSelection(self.RTYPE, 'Output raster type (leave as none to keep original type)', self.TYPE, 0))
        self.addParameter(ParameterString(self.EXTRA,
                          'Additional parameters', '-of GTiff', optional=True))

    def processAlgorithm(self, progress):
        cliPath = '"' + os.path.join(pktoolsUtils.pktoolsPath(), self.cliName()) + '"'
        commands = [cliPath]

        input=self.getParameterValue(self.INPUT)
        commands.append('-i')
        commands.append('"' + input + '"')

        mask=self.getParameterValue(self.MASK)
        maskFiles = mask.split(';')
        for maskFile in maskFiles:
            commands.append('-m')
            commands.append(maskFile)

        commands.append(str(self.getParameterValue(self.MSKBAND)))
        mskband=self.getParameterValue(self.MSKBAND)
        mskbandValues = mskband.split(';')
        for mskbandValue in mskbandValues:
                commands.append('-mskband')
                commands.append(mskbandValue)
        commands.append(str(self.getParameterValue(self.MSKNODATA)))
        msknodata=self.getParameterValue(self.MSKNODATA)
        msknodataValues = msknodata.split(';')
        for msknodataValue in msknodataValues:
                commands.append('-msknodata')
                commands.append(msknodataValue)

        band=self.getParameterValue(self.BAND)
        commands.append("-p")
        commands.append(self.OPERATOR_OPTIONS[self.getParameterValue(self.OPERATOR)])
        nodata=self.getParameterValue(self.NODATA)
        commands.append('-nodata')
        commands.append(nodataValue)
        if self.TYPE[self.getParameterValue(self.RTYPE)] != "none":
            commands.append('-ot')
            commands.append(self.TYPE[self.getParameterValue(self.RTYPE)])
        output=self.getOutputValue(self.OUTPUT)
        if output != "":
            commands.append("-o")
            commands.append('"' + output + '"')

        data=self.getParameterValue(self.DATA)

        extra = str(self.getParameterValue(self.EXTRA))
        if len(extra) > 0:
            commands.append(extra)

        pktoolsUtils.runpktools(commands, progress)

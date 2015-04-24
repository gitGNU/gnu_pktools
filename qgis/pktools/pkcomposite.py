# -*- coding: utf-8 -*-

"""
***************************************************************************
    pkcomposite.py
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

class pkcomposite(pktoolsAlgorithm):

    INPUT = "INPUT"
    OUTPUT = "OUTPUT"
    CRULE_OPTIONS = ["overwrite", "maxndvi", "maxband", "minband", "validband", "mean", "mode", "median", "sum", "minallbands", "maxallbands","stdev"]
    CRULE = "CRULE"
    DX = "DX"
    DY = "DY"
    BB = "BB"
    ULX = "ULX"
    ULY = "ULY"
    LRX = "LRX"
    LRY = "LRY"
    CB = "CB"
    SRCNODATA = "SRCNODATA"
    BNDNODATA = "BNDNODATA"
    DSTNODATA = "DSTNODATA"
    MINGUI = "MINGUI"
    MAXGUI = "MAXGUI"
    RESAMPLE_OPTIONS = ['near', 'bilinear']
    RESAMPLE = "RESAMPLE"
    RTYPE = 'RTYPE'
    TYPE = ['none', 'Byte','Int16','UInt16','UInt32','Int32','Float32','Float64','CInt16','CInt32','CFloat32','CFloat64']
    EXTRA = 'EXTRA'

    def defineCharacteristics(self):
        self.name = "pkcomposite"
        self.group = "Tools"
        self.addParameter(ParameterMultipleInput(self.INPUT, 'Input layer raster data set',ParameterMultipleInput.TYPE_RASTER))
        self.addParameter(ParameterSelection(self.CRULE,"composite rule",self.CRULE_OPTIONS, 0))
        self.addOutput(OutputRaster(self.OUTPUT, "Output raster data set"))
        self.addParameter(ParameterSelection(self.RTYPE, 'Output raster type (leave as none to keep original type)', self.TYPE, 0))
        self.addParameter(ParameterString(self.DX, "Output resolution in x (leave as none to keep original resolution)","none"))
        self.addParameter(ParameterString(self.DY, "Output resolution in y (leave as none to keep original resolution)","none"))
        self.addParameter(ParameterBoolean(self.BB, "user bounding box (do not select to keep original bounding box)", False))
        self.addParameter(ParameterString(self.ULX, "Upper left x value bounding box (select user bounding box)"))
        self.addParameter(ParameterString(self.ULY, "Upper left y value bounding box (select user bounding box)"))
        self.addParameter(ParameterString(self.LRX, "Lower right x value bounding box (select user bounding box)"))
        self.addParameter(ParameterString(self.LRY, "Lower right y value bounding box (select user bounding box)"))
        self.addParameter(ParameterString(self.CB, "band index(es) used for the composite rule","0"))
        self.addParameter(ParameterString(self.SRCNODATA, "invalid value for input raster dataset","none"))
        self.addParameter(ParameterString(self.BNDNODATA, "Band(s) in input image to check if pixel is valid","0"))
        self.addParameter(ParameterString(self.DSTNODATA, "nodata value(s) to put in output raster dataset if not valid or out of bounds","0"))
        self.addParameter(ParameterString(self.MINGUI, "flag values smaller or equal to this value as invalid","none"))
        self.addParameter(ParameterString(self.MAXGUI, "flag values smaller or equal to this value as invalid","none"))
        self.addParameter(ParameterSelection(self.RESAMPLE,"resampling method",self.RESAMPLE_OPTIONS, 0))
        self.addParameter(ParameterString(self.EXTRA,
                          'Additional parameters', '', optional=True))

    def processAlgorithm(self, progress):
        commands = [os.path.join(pktoolsUtils.pktoolsPath(), "bin", "pkcomposite")]
#        commands = [" echo pkcomposite "]
        input=self.getParameterValue(self.INPUT)
        inputFiles = input.split(';')
        for inputFile in inputFiles:
            commands.append('-i')
            commands.append(inputFile)

        if self.TYPE[self.getParameterValue(self.RTYPE)] != "none":
            commands.append('-ot')
            commands.append(self.TYPE[self.getParameterValue(self.RTYPE)])
        output=self.getParameterValue(self.OUTPUT)
        if output != "":
            commands.append("-o")
            commands.append(self.getOutputValue(self.OUTPUT))
        commands.append("-cr")
        commands.append(self.CRULE_OPTIONS[self.getParameterValue(self.CRULE)])
        if self.getParameterValue(self.DX) != "none":
            commands.append("-dx")
            commands.append(str(self.getParameterValue(self.DX)))
        if self.getParameterValue(self.DY) != "none":
            commands.append("-dy")
            commands.append(str(self.getParameterValue(self.DY)))
        if self.getParameterValue(self.BB):
            commands.append("-ulx")
            commands.append(str(self.getParameterValue(self.ULX)))
            commands.append("-uly")
            commands.append(str(self.getParameterValue(self.ULY)))
            commands.append("-lrx")
            commands.append(str(self.getParameterValue(self.LRX)))
            commands.append("-lry")
            commands.append(str(self.getParameterValue(self.LRY)))
        commands.append("-cb")
        commands.append(str(self.getParameterValue(self.CB)))
        srcnodata=self.getParameterValue(self.SRCNODATA)
        if srcnodata != "none":
            srcnodataValues = srcnodata.split(';')
            for srcnodataValue in srcnodataValues:
                commands.append('-srcnodata')
                commands.append(srcnodataValue)
        bndnodata=self.getParameterValue(self.BNDNODATA)
        bndnodataValues = bndnodata.split(';')
        for bndnodataValue in bndnodataValues:
            commands.append('-bndnodata')
            commands.append(bndnodataValue)
        dstnodata=self.getParameterValue(self.DSTNODATA)
        dstnodataValues = dstnodata.split(';')
        for dstnodataValue in dstnodataValues:
            commands.append('-dstnodata')
            commands.append(dstnodataValue)

        minGUI=self.getParameterValue(self.MINGUI)
        if minGUI != "none":
            minValues = minGUI.split(';')
            for minValue in minValues:
                commands.append('-min')
                commands.append(minValue)
        maxGUI=self.getParameterValue(self.MAXGUI)
        if maxGUI != "none":
            maxValues = maxGUI.split(';')
            for maxValue in maxValues:
                commands.append('-max')
                commands.append(maxValue)
        commands.append("-r")
        commands.append(self.RESAMPLE_OPTIONS[self.getParameterValue(self.RESAMPLE)])
        extra = str(self.getParameterValue(self.EXTRA))
        if len(extra) > 0:
            commands.append(extra)

#        commands.append(" |tee /tmp/a")
        pktoolsUtils.runpktools(commands, progress)

# -*- coding: utf-8 -*-

"""
***************************************************************************
    pkextractogr.py
    ---------------------
    Date                 : October 2016
    Copyright            : (C) 2016 by Pieter Kempeneers
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
from processing.core.parameters import ParameterRaster
from processing.core.parameters import ParameterVector
from processing.core.outputs import OutputVector
from processing.core.parameters import ParameterSelection
from processing.core.parameters import ParameterNumber
from processing.core.parameters import ParameterString
from processing.core.parameters import ParameterBoolean
from processing.core.parameters import ParameterExtent

FORMATS = [
    'ESRI Shapefile',
    'GeoJSON',
    'GeoRSS',
    'SQLite',
    'GMT',
    'MapInfo File',
    'INTERLIS 1',
    'INTERLIS 2',
    'GML',
    'Geoconcept',
    'DXF',
    'DGN',
    'CSV',
    'BNA',
    'S57',
    'KML',
    'GPX',
    'PGDump',
    'GPSTrackMaker',
    'ODS',
    'XLSX',
    'PDF',
]
EXTS = [
    '.shp',
    '.geojson',
    '.xml',
    '.sqlite',
    '.gmt',
    '.tab',
    '.ili',
    '.ili',
    '.gml',
    '.txt',
    '.dxf',
    '.dgn',
    '.csv',
    '.bna',
    '.000',
    '.kml',
    '.gpx',
    '.pgdump',
    '.gtm',
    '.ods',
    '.xlsx',
    '.pdf',
]

class pkextractogr(pktoolsAlgorithm):

    INPUT = "INPUT"
    SAMPLE = "SAMPLE"
    OUTPUT = "OUTPUT"
    RULE_OPTIONS = ['centroid', 'mean', 'stdev', 'median', 'proportion', 'count', 'min', 'max', 'mode', 'sum', 'percentile']
    RULE = "RULE"
    BUFFER = "BUFFER"
    SRCNODATA = "SRCNODATA"
    BNDNODATA = "BNDNODATA"
    EXTRA = 'EXTRA'
    FORMAT = "FORMAT"

    def cliName(self):
        return "pkextractogr"

    def defineCharacteristics(self):
        self.name = "extract vector sample from raster"
        self.group = "[pktools] raster/vector"
        self.addParameter(ParameterRaster(self.INPUT, 'Input raster data set'))
        self.addParameter(ParameterFile(self.SAMPLE, "Sample vector data set", False, optional=False))
        self.addParameter(ParameterString(self.LAYERS, "Layer name(s) in sample (leave empty to select all",'', optional=True))
        self.addParameter(ParameterSelection(self.RULE,"extraction rule",self.RULE_OPTIONS, 0))

        self.addOutput(OutputVector(self.OUTPUT, 'Output vector data set'))
        self.addParameter(ParameterSelection(self.FORMAT,
                          'Destination Format', FORMATS))

        self.addParameter(ParameterNumber(self.BUFFER, "Buffer for calculating statistics for point features",0,19,0))
        self.addParameter(ParameterString(self.SRCNODATA, "invalid value(s) for input raster dataset (e.g., 0;255)","none"))
        self.addParameter(ParameterString(self.BNDNODATA, "Band(s) in input image to check if pixel is valid (e.g., 0;1)","0"))

        self.addParameter(ParameterString(self.EXTRA,
                          'Additional parameters', '', optional=True))

    def processAlgorithm(self, progress):
        cliPath = '"' + os.path.join(pktoolsUtils.pktoolsPath(), self.cliName()) + '"'
        commands = [cliPath]

        input=self.getParameterValue(self.INPUT)
        commands.append('-i')
        commands.append('"' + input + '"')

        sample=self.getParameterValue(self.SAMPLE)
        commands.append('-s')
        commands.append(sample)

        layer=self.getParameterValue(self.LAYERS)
        if layer != '':
            layerValues = layer.split(';')
            for layerValue in layerValues:
                commands.append('-ln')
                commands.append(layerValue)

        commands.append("-r")
        commands.append(self.RULE_OPTIONS[self.getParameterValue(self.RULE)])

        output = self.getOutputFromName(self.OUTPUT)
        outFile = output.value
        formatIdx = self.getParameterValue(self.FORMAT)
        outFormat = '"' + FORMATS[formatIdx] + '"'
        commands.append('-f')
        commands.append(outFormat)
        ext = EXTS[formatIdx]
        if not outFile.endswith(ext):
            outFile += ext
            output.value = outFile
        commands.append('-o')
        commands.append('"' + outFile + '"')

        buffer=self.getParameterValue(self.BUFFER)
        if buffer > 1:
            commands.append("-buf")
            commands.append(str(buffer))

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

        extra = str(self.getParameterValue(self.EXTRA))
        if len(extra) > 0:
            commands.append(extra)

        pktoolsUtils.runpktools(commands, progress)

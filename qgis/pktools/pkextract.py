# -*- coding: utf-8 -*-

"""
***************************************************************************
    pkextract.py
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

class pkextract(pktoolsAlgorithm):

    INPUT = "INPUT"
    SAMPLE = "SAMPLE"
    ITERATE = "ITERATE"
    OUTPUT = "OUTPUT"
    
    RULE_OPTIONS = ['centroid', 'point', 'mean', 'proportion', 'min', 'max', 'mode', 'sum', 'median', 'stdev', 'percentile']

    RULE = "RULE"
    POLYGON = "POLYGON"
    BUFFER = "BUFFER"
    SRCNODATA = "SRCNODATA"
    BNDNODATA = "BNDNODATA"

    EXTRA = 'EXTRA'
    
    FORMAT = "FORMAT"

    def cliName(self):
        return "pkextract"

    def defineCharacteristics(self):
        self.name = "extract vector sample from raster"
        self.group = "[pktools] raster/vector"
        self.addParameter(ParameterRaster(self.INPUT, 'Input raster data set'))
        self.addParameter(ParameterVector(self.SAMPLE, 'Sample vector data set'))
        self.addParameter(ParameterBoolean(self.ITERATE, "Iterate over all layers",True))
        self.addParameter(ParameterSelection(self.RULE,"composite rule",self.RULE_OPTIONS, 0))

        self.addOutput(OutputVector(self.OUTPUT, 'Output vector data set'))
        self.addParameter(ParameterSelection(self.FORMAT,
                          'Destination Format', FORMATS))

        self.addParameter(ParameterBoolean(self.POLYGON, "Create OGRPolygon as geometry instead of OGRPoint",False))
        self.addParameter(ParameterNumber(self.BUFFER, "Buffer for calculating statistics for point features",1,25,1))
        self.addParameter(ParameterString(self.SRCNODATA, "invalid value(s) for input raster dataset (e.g., 0;255)","none"))
        self.addParameter(ParameterString(self.BNDNODATA, "Band(s) in input image to check if pixel is valid (e.g., 0;1)","0"))

        self.addParameter(ParameterString(self.EXTRA,
                          'Additional parameters', '', optional=True))

    def processAlgorithm(self, progress):

        commands = [os.path.join(pktoolsUtils.pktoolsPath(), self.cliName())]

        input=self.getParameterValue(self.INPUT)
        commands.append('-i')
        commands.append(input)

        sample=self.getParameterValue(self.SAMPLE)
        if self.getParameterValue(self.ITERATE):
            samplename=str(sample)[:str(sample).find('|')]
        else:
            samplename=str(sample).replace("|layername"," -ln")
        commands.append('-s')
        commands.append(samplename)

        if self.getParameterValue(self.POLYGON):
            commands.append("-polygon")

        commands.append("-r")
        commands.append(self.RULE_OPTIONS[self.getParameterValue(self.RULE)])

        output = self.getOutputFromName(self.OUTPUT)
        outFile = output.value
        formatIdx = self.getParameterValue(self.FORMAT)
        outFormat = FORMATS[formatIdx].replace(' ','\\ ')
        commands.append('-f')
        commands.append(outFormat)
        ext = EXTS[formatIdx]
        if not outFile.endswith(ext):
            outFile += ext
            output.value = outFile
        commands.append('-o')
        commands.append(outFile)

        if self.getParameterValue(self.POLYGON):
            commands.append("-polygon")
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

        f=open('/tmp/a','w')
        for item in commands:
            print >> f, item
        f.close()
        pktoolsUtils.runpktools(commands, progress)

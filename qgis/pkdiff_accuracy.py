# -*- coding: utf-8 -*-

"""
***************************************************************************
    pkdiff_accuracy.py
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

from PyQt4.QtCore import QVariant
from qgis.core import QgsField

from pktoolsUtils import pktoolsUtils
from pktoolsAlgorithm import pktoolsAlgorithm
from processing.core.parameters import ParameterRaster
from processing.core.parameters import ParameterVector
from processing.core.outputs import OutputVector
from processing.core.outputs import OutputFile
#from processing.core.outputs import OutputTable
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

class pkdiff_accuracy(pktoolsAlgorithm):

    INPUT = "INPUT"
    REFERENCE = "REFERENCE"
    ITERATE = "ITERATE"
    LABELREF = "LABELREF"
    NODATA = "NODATA"
#    TABLE = 'TABLE'
    OUTPUT = "OUTPUT"
    CMOUTPUT = "CMOUTPUT"
    CMFORMAT_OPTIONS = ["ascii", "latex"]
    CMFORMAT = "CMFORMAT"

    FORMAT = "FORMAT"
    LABELCLASS = "LABELCLASS"
    EXTRA = 'EXTRA'

    def cliName(self):
        return "pkdiff"

    def defineCharacteristics(self):
        self.name = "Accuracy assessment with ground reference"
        self.group = "[pktools] supervised classification"
        self.addParameter(ParameterRaster(self.INPUT, 'Classification result (raster map)'))
        self.addParameter(ParameterVector(self.REFERENCE, 'Labeled reference vector data set'))
        self.addParameter(ParameterBoolean(self.ITERATE, "Iterate over all layers",True))
        self.addParameter(ParameterString(self.LABELREF, "Attribute name of the reference label","label"))
        self.addParameter(ParameterString(self.NODATA, "No data value(s) in input or reference dataset to ignore (e.g., 0;255)","0"))
        self.addOutput(OutputFile(self.CMOUTPUT, self.tr("Confusion matrix output file ")))
        self.addParameter(ParameterSelection(self.CMFORMAT,"Format for confusion matrix output",self.CMFORMAT_OPTIONS, 0))

#        self.addOutput(OutputTable(self.TABLE, self.tr('Confusion matrix table')))
        self.addOutput(OutputVector(self.OUTPUT, 'Assessment output vector data set'))
        self.addParameter(ParameterSelection(self.FORMAT,
                          'Assessment output vector Format', FORMATS))
        self.addParameter(ParameterString(self.LABELCLASS, "Attribute name of classified (map) label","class"))
        self.addParameter(ParameterString(self.EXTRA,
                          'Additional parameters', '', optional=True))

    def processAlgorithm(self, progress):
        cliPath = '"' + os.path.join(pktoolsUtils.pktoolsPath(), self.cliName()) + '"'
        commands = [cliPath]
        #outputtable = self.getOutputFromName(self.TABLE)

        input=self.getParameterValue(self.INPUT)
        commands.append('-i')
        commands.append('"' + input + '"')

        reference=self.getParameterValue(self.REFERENCE)
        if self.getParameterValue(self.ITERATE):
            if str(reference).find('|')>0:
                referencename=str(reference)[:str(reference).find('|')]
            else:
                referencename=str(reference)
        else:
            referencename=str(reference).replace("|layername"," -ln")
        commands.append('-ref')
        commands.append(referencename)

        commands.append('-lr');
        commands.append(self.getParameterValue(self.LABELREF))

        nodata=self.getParameterValue(self.NODATA)
        if nodata != "none":
            nodataValues = nodata.split(';')
            for nodataValue in nodataValues:
                commands.append('-nodata')
                commands.append(nodataValue)

        commands.append("-cm")
        commands.append("-cmf")
        commands.append(self.CMFORMAT_OPTIONS[self.getParameterValue(self.CMFORMAT)])
        commands.append("-cmo")
        commands.append(self.getOutputValue(self.CMOUTPUT))

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
        commands.append('-lc');
        commands.append(self.getParameterValue(self.LABELCLASS))

        extra = str(self.getParameterValue(self.EXTRA))
        if len(extra) > 0:
            commands.append(extra)

        pktoolsUtils.runpktools(commands, progress)

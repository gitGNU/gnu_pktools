# -*- coding: utf-8 -*-

"""
***************************************************************************
    pktoolsAlgorithm.py
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

from processing.core.GeoAlgorithm import GeoAlgorithm

import os
from PyQt4 import QtGui
from pktools.pktoolsUtils import pktoolsUtils

from processing.core.parameters import ParameterFile
from processing.core.parameters import ParameterBoolean
from processing.core.parameters import ParameterNumber
from processing.core.parameters import ParameterString
from processing.core.parameters import ParameterSelection
from processing.core.outputs import OutputFile
from processing.core.outputs import OutputRaster
from processing.core.outputs import OutputVector

class pktoolsAlgorithm(GeoAlgorithm):

    def getIcon(self):
        filepath = os.path.dirname(__file__) + "/logo.png"
        return QtGui.QIcon(filepath)

    def checkBeforeOpeningParametersDialog(self):
            path = pktoolsUtils.pktoolsPath()
            if path == "":
                return "pktools folder is not configured.\nPlease configure it before running pktools algorithms."
    
    def help(self):
        helpUrl = 'http://pktools.nongnu.org/'
        return False, helpUrl        

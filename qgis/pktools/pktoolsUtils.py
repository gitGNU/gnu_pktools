# -*- coding: utf-8 -*-

"""
***************************************************************************
    pktoolsUtils.py
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

from PyQt4.QtCore import *
from PyQt4.QtGui import *
import subprocess
from processing.core.ProcessingLog import ProcessingLog
from processing.core.ProcessingConfig import ProcessingConfig

class pktoolsUtils():

    PKTOOLS_FOLDER = "PKTOOLS_FOLDER"

    @staticmethod
    def pktoolsPath():
        folder = ProcessingConfig.getSetting(pktoolsUtils.PKTOOLS_FOLDER)
        if folder == None:
            folder =""

        return folder

    @staticmethod
    def runpktools(commands, progress):
        loglines = []
        loglines.append("pktools execution console output")
        commandline = " ".join(commands)
        proc = subprocess.Popen(commandline, shell=True, stdout=subprocess.PIPE, stdin=subprocess.PIPE,stderr=subprocess.STDOUT, universal_newlines=False).stdout
        for line in iter(proc.readline, ""):
            loglines.append(line)
        ProcessingLog.addToLog(ProcessingLog.LOG_INFO, loglines)

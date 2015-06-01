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
import os
import subprocess
from processing.core.ProcessingLog import ProcessingLog
from processing.core.ProcessingConfig import ProcessingConfig

class pktoolsUtils():

    PKTOOLS_FOLDER = "PKTOOLS_FOLDER"

    @staticmethod
    def pktoolsPath():
        folder = ProcessingConfig.getSetting(pktoolsUtils.PKTOOLS_FOLDER)
        if folder == None:
            folder = unicode('/usr/local/bin')
        return os.path.abspath(folder)

    @staticmethod
    def runpktools(commands, progress):
        settings = QSettings()#from gdal
        loglines = []
        loglines.append("pktools execution console output")
        loglines.append(commands)
        progress.setInfo('pktools command:')
        commandline = " ".join(commands)
        progress.setCommand(commandline)
        proc = subprocess.Popen(
            commandline,
            shell=True,
            stdout=subprocess.PIPE,
            stdin=open(os.devnull),
            stderr=subprocess.STDOUT,
            universal_newlines=True,
        ).stdout
        progress.setInfo('pktools command output:')

        #original
        #proc = subprocess.Popen(commandline, shell=True, stdout=subprocess.PIPE, stdin=subprocess.PIPE,stderr=subprocess.STDOUT, universal_newlines=True).stdout
        #1
        #proc = subprocess.Popen(commandline, shell=True, stdout=subprocess.PIPE, stdin=subprocess.PIPE,stderr=subprocess.PIPE, universal_newlines=True).stdout
        #from gdal
        
        for line in iter(proc.readline, ""):
            progress.setConsoleInfo(line)
            loglines.append(line)
        ProcessingLog.addToLog(ProcessingLog.LOG_INFO, loglines)
#        except Exception, e:
#            ProcessingLog.addToLog(ProcessingLog.LOG_ERROR,
#                self.tr('Error in pktools algorithm: %s\n%s' % (descriptionFile, str(e))))


        ProcessingLog.addToLog(ProcessingLog.LOG_INFO, commandline)
        pktoolsUtils.consoleOutput = loglines

#    @staticmethod
#    def getConsoleOutput():
#        return pktoolsUtils.consoleOutput

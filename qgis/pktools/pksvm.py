# -*- coding: utf-8 -*-

"""
***************************************************************************
    pksvm.py
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
from processing.core.parameters import ParameterVector
from processing.core.parameters import ParameterRaster
from processing.core.outputs import OutputRaster
from processing.core.parameters import ParameterSelection
from processing.core.parameters import ParameterNumber
from processing.core.parameters import ParameterString
from processing.core.parameters import ParameterBoolean
from processing.core.parameters import ParameterExtent

class pksvm(pktoolsAlgorithm):

    INPUT = "INPUT"
    TRAINING = "TRAINING"
    LABEL = "LABEL"
    CV = "CV"
    GAMMA = "GAMMA"
    COST = "COST"
    OUTPUT = "OUTPUT"
    MASK = "MASK"
    MSKNODATA = "MSKNODATA"
    NODATA = "NODATA"

#    SVM_TYPE_OPTIONS = ["C_SVC", "nu_SVC,one_class", "epsilon_SVR", "nu_SVR"]
#    KERNEL_TYPE_OPTIONS = ["linear", "polynomial", "radial", "sigmoid"]
    EXTRA = 'EXTRA'

    def defineCharacteristics(self):
        self.name = "Support vector machine"
        self.group = "[pktools] supervised classification"
        self.addParameter(ParameterRaster(self.INPUT, 'Input layer raster data set',ParameterRaster,""))
        self.addParameter(ParameterVector(self.TRAINING, 'Training vector file. A single vector file contains all training features (must be set as: b0, b1, b2,...) for all classes (class numbers identified by label option). Use multiple training files for bootstrap aggregation (alternative to the bag and bsize options, where a random subset is taken from a single training file'))
        self.addParameter(ParameterString(self.LABEL, "Attribute name for class label in training vector file","label"))
        self.addParameter(ParameterBoolean(self.CV, "Two-fold cross validation mode",False))
        self.addParameter(ParameterNumber(self.GAMMA, "Gamma in kernel function",1.0))
        self.addParameter(ParameterNumber(self.COST, "The parameter C of C_SVC",1000.0))
        self.addParameter(ParameterRaster(self.MASK, "Use the first band of the specified file as a validity mask. Nodata values can be set with the option msknodata.",""))
        self.addParameter(ParameterString(self.MSKNODATA, "Mask value(s) not to consider for classification (use negative values if only these values should be taken into account). Values will be taken over in classification image (e.g., 0;255)","0"))
        self.addOutput(OutputRaster(self.OUTPUT, "Output raster data set"))
        self.addParameter(ParameterString(self.NODATA, "Nodata value to put where image is masked as nodata","0"))
        self.addParameter(ParameterString(self.EXTRA,
                          'Additional parameters', '', optional=True))

#        self.addParameter(ParameterSelection(self.KERNEL_TYPE,"Type of kernel function (linear,polynomial,radial,sigmoid)",self.KERNEL_TYPE_OPTIONS, 2))
#        self.addParameter(ParameterSelection(self.SVM_TYPE,"Type of SVM (C_SVC, nu_SVC,one_class, epsilon_SVR, nu_SVR)",self.SVM_TYPE_OPTIONS, 0))

    def processAlgorithm(self, progress):
        commands = [os.path.join(pktoolsUtils.pktoolsPath(), "pksvm")]
#        commands = [" echo pksvm "]
        input=self.getParameterValue(self.INPUT)
        if input != "":
            commands.append('-i')
            commands.append(input)
        commands.append('-t')
        commands.append(self.getParameterValue(self.TRAINING))
        commands.append('-label')
        commands.append(self.getParameterValue(self.LABEL))
        if self.getParameterValue(self.CV):
            commands.append("-cv 2")
        commands.append('-g')
        commands.append(self.getParameterValue(self.GAMMA))
        commands.append('-cc')
        commands.append(self.getParameterValue(self.COST))

        mask = str(self.getParameterValue(self.MASK))
        if len(mask) > 0:
            commands.append('-m')
            commands.append(mask)
            msknodata=self.getParameterValue(self.MSKNODATA)
            msknodataValues = msknodata.split(';')
            for msknodataValue in msknodataValues:
                commands.append('-msknodata')
                commands.append(msknodataValue)
        commands.append('-nodata')
        commands.append(self.getParameterValue(self.NODATA))

        output=self.getParameterValue(self.OUTPUT)
        if output != "":
            commands.append("-o")
            commands.append(self.getOutputValue(self.OUTPUT))
        extra = str(self.getParameterValue(self.EXTRA))
        if len(extra) > 0:
            commands.append(extra)

#        commands.append(" |tee /tmp/a")
        pktoolsUtils.runpktools(commands, progress)

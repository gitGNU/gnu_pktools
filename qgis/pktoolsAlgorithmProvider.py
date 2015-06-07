# -*- coding: utf-8 -*-

"""
***************************************************************************
    pktoolsAlgorithmProvider.py
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


#from pktools.ExampleAlgorithm import ExampleAlgorithm
#raster utilities
from pktools.pkcomposite import pkcomposite
from pktools.pkcrop import pkcrop
from pktools.pkreclass import pkreclass
from pktools.pkgetmask import pkgetmask
from pktools.pksetmask import pksetmask
#raster/vector utilities
from pktools.pkextract import pkextract
from pktools.pkextract_grid import pkextract_grid
from pktools.pkextract_random import pkextract_random
#Supervised classification utilities
from pktools.pksvm import pksvm
from pktools.pkdiff_accuracy import pkdiff_accuracy
#LiDAR utilities
from pktools.pklas2img import pklas2img
from pktools.pkfilterdem import pkfilterdem
#filter utilities
from pktools.pkfilter_spectral import pkfilter_spectral
from pktools.pkfilter_spatial import pkfilter_spatial

from processing.core.AlgorithmProvider import AlgorithmProvider
from processing.core.ProcessingConfig import Setting, ProcessingConfig
import os
from PyQt4 import QtGui
from pktools.pktoolsUtils import pktoolsUtils


class pktoolsAlgorithmProvider(AlgorithmProvider):

    MY_DUMMY_SETTING = "MY_DUMMY_SETTING"

    def __init__(self):
        AlgorithmProvider.__init__(self)
        # deactivate provider by default
        self.activate = True
        # load algorithms
#        self.alglist = [pkinfo()]
        self.alglist = [pkreclass(),pkcrop(),pkcomposite(),pkgetmask(),pksetmask(),pkextract(),pkextract_grid(),pkextract_random(),pksvm(),pkdiff_accuracy(),pklas2img(),pkfilterdem(),pkfilter_spectral(),pkfilter_spatial()]
        # pktools = [pkinfo()]
        # for alg in pktools:
        #     alg.group = "pktools"
        #     self.alglist.extend(pktools)
        for alg in self.alglist:
            alg.provider = self

    def initializeSettings(self):
        '''In this method we add settings needed to configure our provider.
        Do not forget to call the parent method, since it takes care or
        automatically adding a setting for activating or deactivating the
        algorithms in the provider
        '''
        AlgorithmProvider.initializeSettings(self)
        ProcessingConfig.addSetting(Setting(self.getDescription(), pktoolsUtils.PKTOOLS_FOLDER, "pktools folder", pktoolsUtils.pktoolsPath()))

#        ProcessingConfig.addSetting(Setting("Example algorithms", pktoolsAlgorithmProvider.MY_DUMMY_SETTING, "Example setting", "Default value"))
 #       '''To get the parameter of a setting parameter, use
#        ProcessingConfig.getSetting(name_of_parameter)
#        '''

    def unload(self):
        '''Setting should be removed here, so they do not appear anymore
        when the plugin is unloaded'''
        AlgorithmProvider.unload(self)
        ProcessingConfig.removeSetting(pktoolsAlgorithmProvider.MY_DUMMY_SETTING)

    def getName(self):
        '''This is the name that will appear on the toolbox group.
        It is also used to create the command line name of all the algorithms
        from this provider
        '''
        return "pktools"

    def getDescription(self):
        '''This is the provired full name.
        '''
        return "Utilities for remote sensing image processing"

    def getIcon(self):
        filepath = os.path.dirname(__file__) + "/logo.png"
        return QtGui.QIcon(filepath)

    def _loadAlgorithms(self):
        '''Here we fill the list of algorithms in self.algs.
        This method is called whenever the list of algorithms should be updated.
        If the list of algorithms can change
        (for instance, if it contains algorithms from user-defined scripts and
        a new script might have been added), you should create the list again
        here.
        In this case, since the list is always the same, we assign from the pre-made list.
        This assignment has to be done in this method even if the list does not change,
        since the self.algs list is cleared before calling this method
        '''
        self.algs = self.alglist

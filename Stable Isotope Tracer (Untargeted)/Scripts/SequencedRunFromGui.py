# -*- coding: utf-8 -*-
"""
Created on Wed Dec 13 13:52:48 2023

@author: Evan Larson, Hegeman & Cohen Labs - University of Minnesota
"""

import os, time, datetime
from Scripts.OverallRun import seqGUIRUN as Run
def SequenceRun():
    start = time.time()
    numfolders = 0
    for folders in os.listdir('.//Sequenced'):
        print('Analyzing Datafiles in folder: {}'.format(str(folders)))
        numfolders = numfolders + 1
        Run(folders, './/Sequenced/{}'.format(str(folders)))
    end = time.time()
    print('Completed! Total Run Time = {} for {} datasets'.format(str(datetime.timedelta(seconds = (end-start))).split(".")[0], str(numfolders)))

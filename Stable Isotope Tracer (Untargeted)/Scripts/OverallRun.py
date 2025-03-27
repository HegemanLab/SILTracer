# -*- coding: utf-8 -*-
"""
Created on Tue Feb 27 15:52:28 2024

@author: Evan Larson, Hegeman & Cohen Labs - University of Minnesota
"""

import datetime, time, pandas as pd, os,sys

from Scripts.FeatureFinderv5multicsv import findfeaturesbyscanparallel as generate
from Scripts.FeatureFinderv5multicsv import MakePkls as mkpkls
from Scripts.Regressionv8multicsvdb import paralleldatabaseforms as merge
from Scripts.Regressionv8multicsvdb import c13c12regressparallel as regress
from Scripts.Regressionv8multicsvdb import cleanupoutputchecklabeling as clean
from Scripts.Regressionv8multicsvdb import targetlist as targeted
from Scripts.UnlabeledLabeledCompare import UnlabeledLabeledCompare as unlabcompare
from Scripts.UnlabeledLabeledCompare import preplabeledfeaturelistforunlabeledcompare as unlabprep
from Scripts.UnlabeledLabeledCompare import removepkls
from Scripts.Regressionv8multicsvdb import combinerts 
def GUIRUN(featuretype, premadelist, premadexlsx, mincarbons, maxcarbons, masswindow, rtbin, massbin, rtwindow, removeplots, C12C13Window, formppm, fully, minR2, minint, label, targetforms, target, minscan, preferCHNO, UnlabeledControl):
    
    wd = './'
    start = time.time()
    
    if featuretype == 'Generate New Feature List':
        
        generate(wd = wd, masswindow = C12C13Window, mincarbons = mincarbons, maxcarbons = maxcarbons, minint = minint, label = label, rtbin = rtbin, massbin = massbin, minscan = minscan, unlabeledcontrol = 'No')
        regress(wd = wd, featurelist = 'GeneratedFeatureList', masswindow = masswindow, rtwindow = rtwindow, minR2 = minR2, removeplots = removeplots, label = label,featuretype = featuretype, unlabeledcontrol = 'No', minint = minint)
      
        merge('GeneratedFeatureList.csv', wd = wd,formppm = formppm, full = fully, label = label, massbin = massbin, preferCHNO = preferCHNO, unlabeledcontrol = 'No')
    
        clean('{}/Output/PeakInfo_RegressResults.xlsx'.format(wd), targetforms, label, wd = wd , unlabeledcontrol = 'No')
    
        if targetforms == 'Yes':
            targeted(target,wd = wd, unlabeledcontrol = 'No', ppm = formppm)
    
        
        if UnlabeledControl == 'No':
            combinerts(wd = wd, rtbin = rtbin)
            removepkls(wd, unlabeledcontrol= UnlabeledControl)
            try:
                os.mkdir('{}/Output/Peaks/OtherOutputFiles'.format(wd))
            except FileExistsError:
                pass
            for filename in ['FullRegressionDataOutput.xlsx','GeneratedFeatureList.csv', 'MeanStandardDevOutput.xlsx', 'PeakInfo_RegressResults.xlsx', 'TrimmedFeatureList.xlsx']:
                try:
                    os.rename('{}/Output/{}'.format(wd, filename), '{}/Output/Peaks/OtherOutputFiles/{}'.format(wd, filename))
                except FileExistsError:
                    os.replace('{}/Output/{}'.format(wd, filename), '{}/Output/Peaks/OtherOutputFiles/{}'.format(wd, filename))
            end = time.time()
            print('Completed! Total Run Time = {}'.format(str(datetime.timedelta(seconds = (end-start))).split(".")[0]))
    if featuretype == 'Work From Premade Feature List':   
        mkpkls(wd = wd, minint = minint, unlabeledcontrol = UnlabeledControl)    
        regress(wd = wd,featurelist = premadelist, masswindow = masswindow, rtwindow = rtwindow, minR2 = minR2, removeplots = removeplots, label = label,featuretype = featuretype, unlabeledcontrol = 'No', minint = minint)
    
        merge(premadexlsx,wd = wd, formppm = formppm, full = fully, label = label, massbin = massbin, preferCHNO = preferCHNO, unlabeledcontrol = 'No')
    
        clean('{}/Output/PeakInfo_RegressResults.xlsx'.format(wd), targetforms, label, wd = wd, unlabeledcontrol = 'No')
    
        if targetforms == 'Yes':
            targeted(target, wd = wd, unlabeledcontrol = 'No', ppm = formppm)
    
        
        if UnlabeledControl == 'No':
            
            combinerts(wd = wd, rtbin = rtbin)
            removepkls(wd, unlabeledcontrol= UnlabeledControl)
            try:
                os.mkdir('{}/Output/Peaks/OtherOutputFiles'.format(wd))
            except FileExistsError:
                pass
            for filename in ['FullRegressionDataOutput.xlsx','GeneratedFeatureList.csv', 'MeanStandardDevOutput.xlsx', 'PeakInfo_RegressResults.xlsx', 'TrimmedFeatureList.xlsx']:
                try:
                    os.rename('{}/Output/{}'.format(wd, filename), '{}/Output/Peaks/OtherOutputFiles/{}'.format(wd, filename))
                except FileExistsError:
                    os.replace('{}/Output/{}'.format(wd, filename), '{}/Output/Peaks/OtherOutputFiles/{}'.format(wd, filename))
            end = time.time()
            print('Completed! Total Run Time = {}'.format(str(datetime.timedelta(seconds = (end-start))).split(".")[0]))
    if UnlabeledControl == 'Yes':
        wd = './'
        start = time.time()
        
        if featuretype == 'Generate New Feature List':
            
            
            unlabprep(wd)
            mkpkls(wd = wd, minint = minint, unlabeledcontrol = UnlabeledControl)
            regress(wd = wd, featurelist = 'GeneratedFeatureList', masswindow = masswindow, rtwindow = rtwindow, minR2 = minR2, removeplots = removeplots, label = label,featuretype = featuretype, unlabeledcontrol = UnlabeledControl, minint = minint)
          
            merge('GeneratedFeatureList.csv', wd = wd,formppm = formppm, full = fully, label = label, massbin = massbin, preferCHNO = preferCHNO, unlabeledcontrol = UnlabeledControl)
        
            clean('{}/Output/Unlabeled Output/PeakInfo_RegressResults.xlsx'.format(wd), target, label, wd = wd , unlabeledcontrol = UnlabeledControl)
        
            if targetforms == 'Yes':
                targeted(target, wd = wd, unlabeledcontrol = UnlabeledControl, ppm = formppm)
            unlabcompare(wd = wd, target = targetforms, rtbin = rtbin, massbin = massbin)
            combinerts(wd = wd, rtbin = rtbin)
            removepkls(wd, unlabeledcontrol= UnlabeledControl)
            try:
                os.mkdir('{}/Output/Peaks/OtherOutputFiles'.format(wd))
            except FileExistsError:
                pass
            for filename in ['FullRegressionDataOutput.xlsx','GeneratedFeatureList.csv', 'MeanStandardDevOutput.xlsx', 'PeakInfo_RegressResults.xlsx', 'TrimmedFeatureList.xlsx']:
                try:
                    os.rename('{}/Output/{}'.format(wd, filename), '{}/Output/Peaks/OtherOutputFiles/{}'.format(wd, filename))
                except FileExistsError:
                    os.replace('{}/Output/{}'.format(wd, filename), '{}/Output/Peaks/OtherOutputFiles/{}'.format(wd, filename))
            try:
                os.mkdir('{}/Output/Unlabeled Output/Peaks/OtherOutputFiles'.format(wd))
            except FileExistsError:
                pass
            for filename in ['FullRegressionDataOutput.xlsx','GeneratedFeatureList.csv', 'MeanStandardDevOutput.xlsx', 'PeakInfo_RegressResults.xlsx', 'TrimmedFeatureList.xlsx']:
                try:
                    os.rename('{}/Output/Unlabeled Output/{}'.format(wd, filename), '{}/Output/Unlabeled Output/Peaks/OtherOutputFiles/{}'.format(wd, filename))
                except FileExistsError:
                    os.replace('{}/Output/Unlabeled Output/{}'.format(wd, filename), '{}/Output/Unlabeled Output/Peaks/OtherOutputFiles/{}'.format(wd, filename))
            end = time.time()
            print('Completed! Total Run Time = {}'.format(str(datetime.timedelta(seconds = (end-start))).split(".")[0]))
        if featuretype == 'Work From Premade Feature List':   
            unlabprep(wd)    
            mkpkls(wd = wd, minint = minint, unlabeledcontrol = UnlabeledControl)
            regress(wd = wd,featurelist = premadelist, masswindow = masswindow, rtwindow = rtwindow, minR2 = minR2, removeplots = removeplots, label = label,featuretype = featuretype, unlabeledcontrol = UnlabeledControl, minint = minint)
            
            merge(premadexlsx,wd = wd, formppm = formppm, full = fully, label = label, massbin = massbin, preferCHNO = preferCHNO, unlabeledcontrol = UnlabeledControl)
        
            clean('{}/Output/Unlabeled Output/PeakInfo_RegressResults.xlsx'.format(wd), target, label, wd = wd, unlabeledcontrol = UnlabeledControl)
        
            if targetforms == 'Yes':
                targeted(target, wd = wd, unlabeledcontrol = UnlabeledControl, ppm = formppm)
            unlabcompare(wd = wd, target = targetforms, rtbin = rtbin, massbin = massbin)
            combinerts(wd = wd, rtbin = rtbin)
            removepkls(wd, unlabeledcontrol= UnlabeledControl)
            try:
                os.mkdir('{}/Output/Peaks/OtherOutputFiles'.format(wd))
            except FileExistsError:
                pass
            for filename in ['FullRegressionDataOutput.xlsx','GeneratedFeatureList.csv', 'MeanStandardDevOutput.xlsx', 'PeakInfo_RegressResults.xlsx', 'TrimmedFeatureList.xlsx']:
                try:
                    os.rename('{}/Output/{}'.format(wd, filename), '{}/Output/Peaks/OtherOutputFiles/{}'.format(wd, filename))
                except FileExistsError:
                    os.replace('{}/Output/{}'.format(wd, filename), '{}/Output/Peaks/OtherOutputFiles/{}'.format(wd, filename))
            try:
                os.mkdir('{}/Output/Unlabeled Output/Peaks/OtherOutputFiles'.format(wd))
            except FileExistsError:
                pass
            for filename in ['FullRegressionDataOutput.xlsx','GeneratedFeatureList.csv', 'MeanStandardDevOutput.xlsx', 'PeakInfo_RegressResults.xlsx', 'TrimmedFeatureList.xlsx']:
                try:
                    os.rename('{}/Output/Unlabeled Output/{}'.format(wd, filename), '{}/Output/Unlabeled Output/Peaks/OtherOutputFiles/{}'.format(wd, filename))
                except FileExistsError:
                    os.replace('{}/Output/Unlabeled Output/{}'.format(wd, filename), '{}/Output/Unlabeled Output/Peaks/OtherOutputFiles/{}'.format(wd, filename))
            end = time.time()
            print('Completed! Total Run Time = {}'.format(str(datetime.timedelta(seconds = (end-start))).split(".")[0]))

def seqGUIRUN(folder, wd):
    start = time.time()
    try:
        config = pd.read_csv('{}/Config_{}.csv'.format(wd, str(folder)), index_col = 0)
    except FileNotFoundError:
        print('Config File in folder: {} not properly named, please change to Config_{}.csv'.format(folder,folder))
        sys.exit()
    featuretype = str(config.loc['Program Function','Parameters'])
    label = str(config.loc['Select Isotopic Label','Parameters'])
    mincarbons = int(config.loc['Minimum Number of Labels','Parameters'])
    maxcarbons = int(config.loc['Maximum Number of Labels','Parameters'])
    C12C13Window = float(config.loc['Mass Window for Label','Parameters'])
    rtbin = float(config.loc['Retention Time Bin Size','Parameters'])
    massbin = float(config.loc['Mass Bin Size','Parameters'])
    minint = float(config.loc['Minimum Intensity for Features','Parameters'])
    minscan = int(config.loc['Minimum Binned Scans for Features','Parameters'])
    rtwindow = float(config.loc['Retention Time Window For Peak Height','Parameters'])
    masswindow = float(config.loc['Regression Mass Window','Parameters'])
    minR2 = float(config.loc['Minimum R2','Parameters'])
    formppm = float(config.loc['PPM Error','Parameters'])
    target = str(config.loc['TargetFormulas','Parameters'])
    removeplots = str(config.loc['RemovePlots','Parameters'])
    fully = str(config.loc['FullyorPartially','Parameters'])
    targetforms = str(config.loc['FormulaListFile','Parameters'])
    UnlabeledControl = str(config.loc['Unlabeled Control', 'Parameters'])
    premadelist = str(config.loc['PremadeList','Parameters']).replace('.csv','')
    premadexlsx = str(config.loc['PremadeList','Parameters'])
    preferCHNO = str(config.loc['PreferCHNO','Parameters'])
    start = time.time()
    
    if featuretype == 'Generate New Feature List':
        
        generate(wd = wd, masswindow = C12C13Window, mincarbons = mincarbons, maxcarbons = maxcarbons, minint = minint, label = label, rtbin = rtbin, massbin = massbin, minscan = minscan, unlabeledcontrol = 'No')
        regress(wd = wd, featurelist = 'GeneratedFeatureList', masswindow = masswindow, rtwindow = rtwindow, minR2 = minR2, removeplots = removeplots, label = label,featuretype = featuretype, unlabeledcontrol = 'No', minint = minint)
      
        merge('GeneratedFeatureList.csv', wd = wd,formppm = formppm, full = fully, label = label, massbin = massbin, preferCHNO = preferCHNO, unlabeledcontrol = 'No')
    
        clean('{}/Output/PeakInfo_RegressResults.xlsx'.format(wd), targetforms, label, wd = wd , unlabeledcontrol = 'No')
    
        if targetforms == 'Yes':
            targeted(target,wd = wd, unlabeledcontrol = 'No', ppm = formppm)
    
        
        if UnlabeledControl == 'No':
            combinerts(wd = wd, rtbin = rtbin)
            removepkls(wd, unlabeledcontrol= UnlabeledControl)
            try:
                os.mkdir('{}/Output/Peaks/OtherOutputFiles'.format(wd))
            except FileExistsError:
                pass
            for filename in ['FullRegressionDataOutput.xlsx','GeneratedFeatureList.csv', 'MeanStandardDevOutput.xlsx', 'PeakInfo_RegressResults.xlsx', 'TrimmedFeatureList.xlsx']:
               try:
                   os.rename('{}/Output/{}'.format(wd, filename), '{}/Output/Peaks/OtherOutputFiles/{}'.format(wd, filename))
               except FileExistsError:
                   os.replace('{}/Output/{}'.format(wd, filename), '{}/Output/Peaks/OtherOutputFiles/{}'.format(wd, filename))
            end = time.time()
            print('Completed! Total Run Time = {}'.format(str(datetime.timedelta(seconds = (end-start))).split(".")[0]))
    if featuretype == 'Work From Premade Feature List':   
        mkpkls(wd = wd, minint = minint, unlabeledcontrol = UnlabeledControl)    
        regress(wd = wd,featurelist = premadelist, masswindow = masswindow, rtwindow = rtwindow, minR2 = minR2, removeplots = removeplots, label = label,featuretype = featuretype, unlabeledcontrol = 'No', minint = minint)
    
        merge(premadexlsx,wd = wd, formppm = formppm, full = fully, label = label, massbin = massbin, preferCHNO = preferCHNO, unlabeledcontrol = 'No')
    
        clean('{}/Output/PeakInfo_RegressResults.xlsx'.format(wd), targetforms, label, wd = wd, unlabeledcontrol = 'No')
    
        if targetforms == 'Yes':
            targeted(target, wd = wd, unlabeledcontrol = 'No', ppm = formppm)
    
        
        if UnlabeledControl == 'No':
            combinerts(wd = wd, rtbin = rtbin)
            removepkls(wd, unlabeledcontrol= UnlabeledControl)
            try:
                os.mkdir('{}/Output/Peaks/OtherOutputFiles'.format(wd))
            except FileExistsError:
                pass
            for filename in ['FullRegressionDataOutput.xlsx','GeneratedFeatureList.csv', 'MeanStandardDevOutput.xlsx', 'PeakInfo_RegressResults.xlsx', 'TrimmedFeatureList.xlsx']:
                try:
                    os.rename('{}/Output/{}'.format(wd, filename), '{}/Output/Peaks/OtherOutputFiles/{}'.format(wd, filename))
                except FileExistsError:
                    os.replace('{}/Output/{}'.format(wd, filename), '{}/Output/Peaks/OtherOutputFiles/{}'.format(wd, filename))
            end = time.time()
            print('Completed! Total Run Time = {}'.format(str(datetime.timedelta(seconds = (end-start))).split(".")[0]))
    if UnlabeledControl == 'Yes':
        start = time.time()
        
        if featuretype == 'Generate New Feature List':
            
            
            unlabprep(wd)
            mkpkls(wd = wd, minint = minint, unlabeledcontrol = UnlabeledControl)
            regress(wd = wd, featurelist = 'GeneratedFeatureList', masswindow = masswindow, rtwindow = rtwindow, minR2 = minR2, removeplots = removeplots, label = label,featuretype = featuretype, unlabeledcontrol = UnlabeledControl, minint = minint)
          
            merge('GeneratedFeatureList.csv', wd = wd,formppm = formppm, full = fully, label = label, massbin = massbin, preferCHNO = preferCHNO, unlabeledcontrol = UnlabeledControl)
        
            clean('{}/Output/Unlabeled Output/PeakInfo_RegressResults.xlsx'.format(wd), target, label, wd = wd , unlabeledcontrol = UnlabeledControl)
        
            if targetforms == 'Yes':
                targeted(target, wd = wd, unlabeledcontrol = UnlabeledControl, ppm = formppm)
            unlabcompare(wd = wd, target = targetforms, rtbin = rtbin, massbin = massbin)
            combinerts(wd = wd, rtbin = rtbin)
            removepkls(wd, unlabeledcontrol= UnlabeledControl)
            try:
                os.mkdir('{}/Output/Peaks/OtherOutputFiles'.format(wd))
            except FileExistsError:
                pass
            for filename in ['FullRegressionDataOutput.xlsx','GeneratedFeatureList.csv', 'MeanStandardDevOutput.xlsx', 'PeakInfo_RegressResults.xlsx', 'TrimmedFeatureList.xlsx']:
                try:
                    os.rename('{}/Output/{}'.format(wd, filename), '{}/Output/Peaks/OtherOutputFiles/{}'.format(wd, filename))
                except FileExistsError:
                    os.replace('{}/Output/{}'.format(wd, filename), '{}/Output/Peaks/OtherOutputFiles/{}'.format(wd, filename))
            try:
                os.mkdir('{}/Output/Unlabeled Output/Peaks/OtherOutputFiles'.format(wd))
            except FileExistsError:
                pass
            for filename in ['FullRegressionDataOutput.xlsx','GeneratedFeatureList.csv', 'MeanStandardDevOutput.xlsx', 'PeakInfo_RegressResults.xlsx', 'TrimmedFeatureList.xlsx']:
                try:
                    os.rename('{}/Output/Unlabeled Output/{}'.format(wd, filename), '{}/Output/Unlabeled Output/Peaks/OtherOutputFiles/{}'.format(wd, filename))
                except FileExistsError:
                    os.replace('{}/Output/Unlabeled Output/{}'.format(wd, filename), '{}/Output/Unlabeled Output/Peaks/OtherOutputFiles/{}'.format(wd, filename))
            end = time.time()
            print('Completed! Total Run Time = {}'.format(str(datetime.timedelta(seconds = (end-start))).split(".")[0]))
        if featuretype == 'Work From Premade Feature List':   
            unlabprep(wd)
            mkpkls(wd = wd, minint = minint, unlabeledcontrol = UnlabeledControl)
            regress(wd = wd,featurelist = premadelist, masswindow = masswindow, rtwindow = rtwindow, minR2 = minR2, removeplots = removeplots, label = label,featuretype = featuretype, unlabeledcontrol = UnlabeledControl, minint = minint)
        
            merge(premadexlsx,wd = wd, formppm = formppm, full = fully, label = label, massbin = massbin, preferCHNO = preferCHNO, unlabeledcontrol = UnlabeledControl)
        
            clean('{}/Output/Unlabeled Output/PeakInfo_RegressResults.xlsx'.format(wd), target, label, wd = wd, unlabeledcontrol = UnlabeledControl)
        
            if targetforms == 'Yes':
                targeted(target, wd = wd, unlabeledcontrol = UnlabeledControl, ppm = formppm)
            unlabcompare(wd = wd, target = targetforms, rtbin = rtbin, massbin = massbin)
            combinerts(wd = wd, rtbin = rtbin)
            removepkls(wd, unlabeledcontrol= UnlabeledControl)
            try:
                os.mkdir('{}/Output/Peaks/OtherOutputFiles'.format(wd))
            except FileExistsError:
                pass
            for filename in ['FullRegressionDataOutput.xlsx','GeneratedFeatureList.csv', 'MeanStandardDevOutput.xlsx', 'PeakInfo_RegressResults.xlsx', 'TrimmedFeatureList.xlsx']:
                try:
                    os.rename('{}/Output/{}'.format(wd, filename), '{}/Output/Peaks/OtherOutputFiles/{}'.format(wd, filename))
                except FileExistsError:
                    os.replace('{}/Output/{}'.format(wd, filename), '{}/Output/Peaks/OtherOutputFiles/{}'.format(wd, filename))
            try:
                os.mkdir('{}/Output/Unlabeled Output/Peaks/OtherOutputFiles'.format(wd))
            except FileExistsError:
                pass
            for filename in ['FullRegressionDataOutput.xlsx','GeneratedFeatureList.csv', 'MeanStandardDevOutput.xlsx', 'PeakInfo_RegressResults.xlsx', 'TrimmedFeatureList.xlsx']:
                try:
                    os.rename('{}/Output/Unlabeled Output/{}'.format(wd, filename), '{}/Output/Unlabeled Output/Peaks/OtherOutputFiles/{}'.format(wd, filename))
                except FileExistsError:
                    os.replace('{}/Output/Unlabeled Output/{}'.format(wd, filename), '{}/Output/Unlabeled Output/Peaks/OtherOutputFiles/{}'.format(wd, filename))
            end = time.time()
            print('Completed! Total Run Time = {}'.format(str(datetime.timedelta(seconds = (end-start))).split(".")[0]))
         
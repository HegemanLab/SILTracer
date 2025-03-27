# -*- coding: utf-8 -*-
"""
Created on Tue May 21 14:23:40 2024

@author: Evan Larson, Hegeman & Cohen Labs - University of Minnesota
"""


import pandas as pd, os, re
from functools import partial


    
def UnlabeledLabeledCompare(wd, target,rtbin, massbin):
    if target == 'Yes':
        LabeledOutput = pd.read_excel('{}/Output/FinalOutput.xlsx'.format(wd))
        Labeled = LabeledOutput.copy()
        Unlabeled = pd.read_excel('{}/Output/Unlabeled Output/FinalOutput.xlsx'.format(wd))
    if target == 'No':
        LabeledOutput = pd.read_excel('{}/Output/FinalOutput.xlsx'.format(wd))
        Labeled = LabeledOutput.copy()
        Unlabeled = pd.read_excel('{}/Output/Unlabeled Output/FinalOutput.xlsx'.format(wd))
    
    Merged = pd.merge(Labeled,Unlabeled, how = 'left', on = 'Unlabeled_Labeled Features', suffixes=['_Labeled','_Unlabeled'])
    
    Merged['SlopeDiff'] = Merged['Slope_Labeled'] - Merged['Slope_Unlabeled']
    
    Merged = Merged[['Unlabeled_Labeled Features', 'SlopeDiff']]
    
    NewOutput = pd.merge(LabeledOutput, Merged, how = 'left', on = 'Unlabeled_Labeled Features')
    NewOutputcopy = NewOutput.copy()
    NewOutputcopy['Plot Hyperlink'] = '=HYPERLINK("Plots\\' + NewOutputcopy['Feature'] + '","Link")'
    NewOutputcopy['Formula Hyperlink'] = '=HYPERLINK("Formulas\\' + NewOutputcopy['MassFeature'] + '.xlsx","Link")'
    NewOutputCleaned = NewOutputcopy.loc[:, ~NewOutputcopy.columns.str.match('Unnamed')].copy()
    NewOutputCleaned.drop_duplicates(inplace = True)
    NewOutputCleaned.to_excel('{}/Output/FinalOutput.xlsx'.format(wd))
    




def roundmass(mass, massbin):
    return(round(round(mass / massbin) * massbin,3))

def preplabeledfeaturelistforunlabeledcompare(wd):
    masslist = pd.read_csv('{}/Output/{}.csv'.format(wd,'GeneratedFeatureList'))
    masslistout = masslist.loc[masslist['NumLabels'] == 1]
    try:
        os.mkdir('{}/Output/Unlabeled Output'.format(wd))
    except FileExistsError:
        pass
    masslistout.to_csv('{}/Output/Unlabeled Output/GeneratedFeatureList.csv'.format(wd))
    
def singletomultfeatlist(file, label, labellow, labelhigh, massbin):
    if label == '13C':
        md = 1.003355
    elif label == '15N':
        md = 0.997035
    elif label == '18O':
        md = 2.004245
    elif label == '34S':
        md = 1.99579
    else:
        md = label
    featurelist = pd.read_csv('.//{}'.format(str(file)))
    masslistout = featurelist[['mzUnlabeled','rtUnlabeled', 'm/zUnlabeledround', 'rtUnlabeledround','pol', 'NumLabels']].loc[featurelist['NumLabels'] == 1].drop('NumLabels', axis = 1)
    masslistoutall = pd.DataFrame()
    massbinpartial = partial(roundmass, massbin = massbin)
    for num in range(labellow, labelhigh + 1):
        temp = masslistout.copy()
        temp['NumLabels'] = num
        temp['rtLabeled'] = temp['rtUnlabeled']
        temp['mzLabeled'] = temp['mzUnlabeled'] + (num * md)
        temp['Feature'] = temp['m/zUnlabeledround'].astype(str) + '_' + temp['rtUnlabeledround'].astype(str) + '_' + temp['mzLabeled'].apply(massbinpartial).astype(str) + '_' + temp['rtUnlabeledround'].astype(str)
        temp['MassFeature'] = temp['m/zUnlabeledround'].astype(str) + '_' + temp['mzLabeled'].apply(massbinpartial).astype(str)
        masslistoutall = pd.concat([masslistoutall, temp], axis = 0)
        del temp
    masslistoutall.to_csv('.//{}_addLabels.csv'.format(str(file.strip('.csv'))))
    
def removepkls(wd, unlabeledcontrol):
    pklregex = re.compile('(.*)(.pkl)')
    for files in os.listdir('{}/Data'.format(str(wd))):
        match = pklregex.search(files)
        if match == None:
            continue
        else:
            os.remove('{}/Data/{}'.format(str(wd),str(files)))
    if unlabeledcontrol == 'Yes':
        for files in os.listdir('{}/Data/Unlabeled Control'.format(str(wd))):
            match = pklregex.search(files)
            if match == None:
                continue
            else:
                os.remove('{}/Data/Unlabeled Control/{}'.format(str(wd),str(files)))
    


# -*- coding: utf-8 -*-
"""
Created on Mon May 15 09:40:49 2023

@author: Evan Larson, Hegeman & Cohen Labs - University of Minnesota
"""

import multiprocessing as mp
from functools import partial
import pandas as pd, os, re, time,shutil, datetime, math
import warnings
from sklearn.exceptions import UndefinedMetricWarning
warnings.filterwarnings("ignore", category=UndefinedMetricWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)
from sklearn.linear_model import LinearRegression
import numpy as np
from matplotlib import pyplot as plt
from Scripts.ChemicalCompositionCalculatorV2 import searchdatabase
# from ChemicalCompositionCalculatorV2 import searchdatabase2
from p_tqdm import p_umap
from statistics import fmean, mean
from scipy import stats
# from numpy import random

    # Clusters features which are highly coorelated within the same RT bin
def combinerts(wd, rtbin):
    # Reads in data from previous output
    data = pd.read_excel('{}/Output/FinalOutput.xlsx'.format(wd))
    rtroundpartial = partial(roundrt, rtbin = rtbin)
    # Rounds RT to the same RT bin to group clusters
    data['rtforclustering'] = data['rtUnlabeled'].apply(rtroundpartial)
    unrts = data['rtforclustering']
    unrtsunique = np.unique(unrts)
    peaks = []
    combinertsparallelpartial = partial(combinertsparallel, wd = wd, threshold = 0.6)
    # Separates each RT cluster into individual groups for submitting to multiprocessing
    for rts in unrtsunique:
        trimdata = data.loc[data['rtforclustering'] == rts]
        trimdict = {str(rts):trimdata}
        peaks.append(trimdict)
    print('Performing Correlation on {} Peaks'.format(str(len(peaks))))
    datalist = p_umap(combinertsparallelpartial, peaks)
    data = pd.concat(datalist)
    
           
    
    data['avg'] = data[['UnlabeledHeight', 'LabeledHeight']].mean(axis = 1)
    data['20%R2_R2 AvgRound'] = data['20%R2_R2 Avg'].round(2)
    data.sort_values(by = ['20%R2_R2 AvgRound','avg'], ascending = False, inplace = True)
    data['rtlink'] = data['rtforclustering'].astype(str).str.replace('.','_')
    data['Plot Hyperlink'] = '=HYPERLINK("Peaks\\Plots\\' + data['Feature'] + '","Link")'
    data['Formula Hyperlink'] = '=HYPERLINK("Peaks\\Formulas\\' + data['MassFeature'] + '.xlsx","Link")'
    data['Peak Hyperlink'] = '=HYPERLINK("Peaks\\' + data['rtlink'].astype(str) + '.xlsx","Link")'
    data = data.drop(labels = ['avg','20%R2_R2 AvgRound', 'rtlink', 'rtforclustering'], axis = 1)
    
    datareduce = data.drop_duplicates(subset = 'PeakFeatureGroup')
    datasort = data.sort_values(by = 'Unlabeled_Labeled Features')
    datareducesort = datareduce.sort_values(by = 'Unlabeled_Labeled Features')
    data = data.drop(labels = 'FeaturesInPeak', axis = 1)
    datasort = datasort.loc[:, ~datasort.columns.str.match('Unnamed')]
    datareducesort = datareducesort.loc[:, ~datareducesort.columns.str.match('Unnamed')]
     
    
    datasort.to_excel('{}/Output/FinalOutput.xlsx'.format(wd))
    datareducesort.to_excel('{}/Output/FinalOutput_PeakReduced.xlsx'.format(wd))    


def combinertsparallel(datadict,wd, threshold):    
    #take each peak, find joined csvs, cross corr, export any that are well correlated, make new df with link to correct file 
    peak = list(datadict.keys())[0]
    data = datadict[peak]

    peakout = peak.replace('.','_')
    data['avg'] = data[['UnlabeledHeight', 'LabeledHeight']].mean(axis = 1)
    data['20%R2_R2 AvgRound'] = data['20%R2_R2 Avg'].round(2)
    datasort = data.sort_values(by = ['20%R2_R2 AvgRound','avg'], ascending = False)
    
    features = data['Unlabeled_Labeled Features']
    featuressort = datasort['Unlabeled_Labeled Features'].to_list()
    data = data.drop(labels = ['avg','20%R2_R2 AvgRound'], axis = 1)
    if len(features) == 1:
        data['HighestIntensityFeature'] = 'N/A'
        data['PeakFeatureGroup'] = 'N/A'
        data['FeaturesInPeak'] = 'N/A'
        data['PeakCorrelationScore'] = 'N/A'
        
        return data
    else: 
        rtdf = pd.DataFrame()
        for feat in features:
            featdf = pd.DataFrame()
            for file in os.listdir('{}/Output/Peaks/JoinedCSV/'.format(wd)):
                if file.startswith(str(feat)) and file.endswith('joined.csv'):
                    read = pd.read_csv('{}/Output/Peaks/JoinedCSV/{}'.format(wd,str(file)), usecols = ['scan', 'int_Unlabeled', 'int_Labeled'])
                    read.set_index('scan', inplace = True)
                    read['int_Sum'] = read['int_Unlabeled'] + read['int_Labeled']
                    read.drop(['int_Unlabeled', 'int_Labeled'], inplace = True, axis = 1)
                    filename = str(file).replace(str(feat)+'_','').replace('joined.csv', '')
                    multindex = pd.MultiIndex.from_tuples([(str(feat),str(filename))])
                    read.columns = multindex
                    featdf = pd.concat([featdf, read], axis = 1, join = 'outer')
            featdfsum = featdf.groupby(axis = 1, level = 0).sum()
            
            rtdf = pd.concat([rtdf, featdfsum], axis = 1, join = 'outer')
            rtdf.sort_index(axis = 0, inplace = True)
            rtdf.fillna(0, inplace = True)

        fullpeakout = pd.DataFrame()
        featcorrdf = rtdf.corr('pearson')
        
        i = 0
        for feat in featuressort:
            
            try:
                featcorr = featcorrdf.loc[featcorrdf[str(feat)]>threshold].index.to_list()
                featcorrdfout = pd.DataFrame(featcorrdf[str(feat)].loc[featcorrdf[str(feat)]>threshold])
                
                
                
            except KeyError:
                continue
            featcorrdfout['Unlabeled_Labeled Features'] = featcorrdfout.index
            featcorrdfout['HighestIntensityFeature'] = str(feat)
            featcorrdfout['PeakFeatureGroup'] = peak + '_' + str(i)
            featcorrdfout['FeaturesInPeak'] = len(featcorrdfout)
            featcorrdfout.rename({str(feat):'PeakCorrelationScore'}, inplace = True, axis = 1)
            featcorrdfout.reset_index(inplace = True, drop = True)
            
            featcorrdf.drop(columns = featcorr, inplace = True)
            featcorrdf.drop(index = featcorr, inplace = True)
            dataout = data.merge(featcorrdfout, how = 'right', on = 'Unlabeled_Labeled Features')
            fullpeakout = pd.concat([fullpeakout, dataout], axis = 0)
            i = i + 1
            
          
        fullpeakout = fullpeakout.copy()
        fullpeakout['Plot Hyperlink'] = '=HYPERLINK("Plots\\' + fullpeakout['Feature'] + '","Link")'
        fullpeakout['Formula Hyperlink'] = '=HYPERLINK("Formulas\\' + fullpeakout['MassFeature'] + '.xlsx","Link")'
        
        fullpeakout.to_excel('{}/Output/Peaks/{}.xlsx'.format(wd,str(peakout)))
        return fullpeakout

def cleanupoutputchecklabeling(datapath, target, label, wd, unlabeledcontrol):
    # Reads in formula/regression/mass data from previous steps and sorts by both feature label and R2 while dropping duplicate mass features with the exact height ratio to remove double counted features that fell outside binning
    data = pd.read_excel(datapath)
    data['sum'] = data['UnlabeledHeight'] + data['LabeledHeight']
    data.sort_values(by = ['sum'], ascending = False, inplace = True)
    data.drop_duplicates(subset = ['MassFeature', 'HeightRatio'], inplace = True)
    data = data.drop('sum', axis = 1)
    datacleaned = data.loc[:, ~data.columns.str.match('Unnamed')]
    datacleanedcopy = data
    datacleaned = data
    # If there is no target list to search, makes a hyperlink in the output file to link directly to the plot folder for each feature
    if target == 'No' and unlabeledcontrol == 'No':
        
        datacleanedcopy['Plot Hyperlink'] = '=HYPERLINK("Peaks\\Plots\\' + datacleaned['Feature'] + '","Link")'
        datacleanedcopy['Formula Hyperlink'] = '=HYPERLINK("Peaks\\Formulas\\' + datacleaned['MassFeature'] + '.xlsx","Link")'
    if unlabeledcontrol == 'No':
        datacleanedcopy.to_excel('{}/Output/FinalOutput.xlsx'.format(wd))
    if unlabeledcontrol == 'Yes':
        datacleanedcopy.to_excel('{}/Output/Unlabeled Output/FinalOutput.xlsx'.format(wd))
        

def targetlist(targetlist, wd, unlabeledcontrol, ppm):
    # Reads in the target list and the output data
    if unlabeledcontrol == 'No':
        output = 'Output'
    if unlabeledcontrol == 'Yes':
        output = 'Output/Unlabeled Output'
    targets = pd.read_excel('.//{}'.format(str(targetlist)))
    data = pd.read_excel('{}/{}/FinalOutput.xlsx'.format(wd, output))
    # os.remove('{}/{}/FinalOutput.xlsx'.format(wd, output))
    # Fills unmatched formulas with an empty space
    targets['M+H'] = targets['Neutral Mass'] + 1.007276
    targets['M+Na'] = targets['Neutral Mass'] + 22.989769
    targets['M-H'] =targets['Neutral Mass'] - 1.007276
    data['Formula'].fillna('', inplace = True)
    datatomatch = data[['Unlabeled_Labeled Features','mzUnlabeled','pol']]
    ppmmil = ppm / 10**6
    Allmatches = pd.DataFrame()
    for rows in range(len(targets)):
        MH = targets.loc[rows,'M+H']
        MNa = targets.loc[rows, 'M+Na']
        MHNeg = targets.loc[rows,'M-H']
        Name = targets.loc[rows,'Name']
        Hmatches = datatomatch.loc[((((datatomatch['mzUnlabeled'].values) - (ppmmil*datatomatch['mzUnlabeled'].values)) <= MH) & (((datatomatch['mzUnlabeled'].values) + (ppmmil*datatomatch['mzUnlabeled'].values)) >= MH) & (datatomatch['pol'] == 'POSITIVE'))] 
        Namatches = datatomatch.loc[((((datatomatch['mzUnlabeled'].values) - (ppmmil*datatomatch['mzUnlabeled'].values)) <= MNa) & (((datatomatch['mzUnlabeled'].values) + (ppmmil*datatomatch['mzUnlabeled'].values)) >= MNa) & (datatomatch['pol'] == 'POSITIVE'))]
        Hnegmatches = datatomatch.loc[((((datatomatch['mzUnlabeled'].values) - (ppmmil*datatomatch['mzUnlabeled'].values)) <= MHNeg) & (((datatomatch['mzUnlabeled'].values) + (ppmmil*datatomatch['mzUnlabeled'].values)) >= MHNeg) & (datatomatch['pol'] == 'NEGATIVE'))]
        Hmatches = Hmatches.copy() 
        Hmatches.loc[:,'[M+H] Tentative Match'] = [Name] * len(Hmatches)
        if len(Hmatches) == 1:
            Hmatches['[M+H] PPM Error'] = round(float((((MH - Hmatches['mzUnlabeled'].values)*10**6)/MH)),3)
        elif len(Hmatches) > 1:
            Hmatches['[M+H] PPM Error'] = round((((MH - Hmatches['mzUnlabeled'].astype(float))*10**6)/MH),3)
        else:
            Hmatches['[M+H] PPM Error'] = [''] * len(Hmatches)
        Namatches = Namatches.copy() 
        Namatches.loc[:,'[M+Na] Tentative Match'] = [Name] * len(Namatches)
        if len(Namatches) == 1:
            Namatches['[M+Na] PPM Error'] = round(float((((MNa - Namatches['mzUnlabeled'].values)*10**6)/MNa)),3)
        elif len(Namatches) > 1:
            Namatches['[M+Na] PPM Error'] = round((((MNa - Namatches['mzUnlabeled'].astype(float))*10**6)/MNa),3)
        else:
            Namatches['[M+Na] PPM Error'] = [''] * len(Namatches)
        Hnegmatches = Hnegmatches.copy() 
        Hnegmatches.loc[:,'[M-H] Tentative Match'] = [Name] * len(Hnegmatches)
        if len(Hnegmatches) == 1:
            Hnegmatches['[M-H] PPM Error'] = round(float((((MHNeg - Hnegmatches['mzUnlabeled'].values)*10**6)/MHNeg)),3)
        elif len(Hnegmatches) > 1:
            Namatches['[M-H] PPM Error'] = round((((MNa - Hnegmatches['mzUnlabeled'].astype(float))*10**6)/MNa),3)
        else:
            Hnegmatches['[M-H] PPM Error'] = [''] * len(Hnegmatches)
        Hmatches.drop(labels = ['mzUnlabeled', 'pol'], axis = 1, inplace = True)
        Hmatches.fillna('', inplace = True)
        Namatches.drop(labels = ['mzUnlabeled', 'pol'], axis = 1, inplace = True)
        Namatches.fillna('', inplace = True)
        Hnegmatches.drop(labels = ['mzUnlabeled', 'pol'], axis = 1, inplace = True)
        Hnegmatches.fillna('', inplace = True)
        Namatches['[M+Na] Tentative Match'] = Namatches['[M+Na] Tentative Match'].astype(str) + '*'
        Hnegmatches['[M-H] Tentative Match'] = Hnegmatches['[M-H] Tentative Match'].astype(str) + '**'
        Namatches['[M+Na] Tentative Match'] = Namatches['[M+Na] Tentative Match'].replace('*','')
        Hnegmatches['[M-H] Tentative Match'] = Hnegmatches['[M-H] Tentative Match'].replace('**','')
        # matches = Hmatches.merge(Namatches,on = 'Unlabeled_Labeled Features' ).merge(Hnegmatches, on = 'Unlabeled_Labeled Features')
        
        
        Allmatches = pd.concat([Allmatches,Hmatches,Hnegmatches,Namatches], axis = 0, ignore_index = True )
        
    # Merges the target list with the output data based upon formulas, adds in a hyperlink to the plot folders and exports the final data output
    
    Allmatches['Tentative Matches'] = Allmatches['[M+H] Tentative Match'].fillna('').astype(str) + '' + Allmatches['[M+Na] Tentative Match'].fillna('').astype(str) + ' ' + Allmatches['[M-H] Tentative Match'].fillna('').astype(str)
    Allmatches['Tentative Matches PPM Error'] = Allmatches['[M+H] PPM Error'].fillna('').astype(str) + '' + Allmatches['[M+Na] PPM Error'].fillna('').astype(str) + '' + Allmatches['[M-H] PPM Error'].fillna('').astype(str)
    try:
        Allmatches['Tentative Matches PPM Error'] = pd.to_numeric(Allmatches['Tentative Matches PPM Error'])
    except TypeError:
        pass
    Allmatches.drop(labels = ['[M+H] Tentative Match','[M+Na] Tentative Match','[M-H] Tentative Match','[M+H] PPM Error','[M+Na] PPM Error','[M-H] PPM Error'], axis = 1, inplace = True)
    merged = pd.merge(data, Allmatches, how = 'left', on = 'Unlabeled_Labeled Features')
    mergedcleaned = merged.loc[:, ~merged.columns.str.match('Unnamed')]
    # This line is unecessary but pandas gets mad if you dont have it
    mergedcleaned = mergedcleaned.copy()
    if unlabeledcontrol == 'No':
        mergedcleaned['Plot Hyperlink'] = '=HYPERLINK("Peaks\\Plots\\' + mergedcleaned['Feature'] + '","Link")'
        mergedcleaned['Formula Hyperlink'] = '=HYPERLINK("Peaks\\Formulas\\' + mergedcleaned['MassFeature'] + '.xlsx","Link")'
    mergedcleaned.to_excel('{}/{}/FinalOutput.xlsx'.format(wd, output))



def databaseformsearch(joinedround, formppm, full, label, massbin, preferCHNO, wd):
    # Dataframe to collect all results
    summedformula = pd.DataFrame()
    count = 0
    # Resets index for iteration
    joinedround.reset_index(inplace = True, drop = True)
    # Iterates through all submited regression/mass results to find database matches, if split into individual features per dataframe then this is not necessary. However, if splitting the regression/mass data differently this will cover that possibility
    for i in range(len(joinedround)):
        count = count + 1
        # Submits individual mass information for current feature to find any database matches and assignes the unique mass feature label for future joining, returns the 
        formulas = searchdatabase(mass = joinedround['mzUnlabeled'][i], c13mass = joinedround['mzLabeled'][i], polarity = joinedround['pol'][i], carbons = joinedround['NumLabels'][i], ppm = formppm, massbin = massbin, preferCHNO = preferCHNO, label = label, full = full, wd = wd)
        formulas['MassFeature'] = joinedround['MassFeature'][i]
        summedformula = pd.concat([summedformula, formulas], axis = 0)
    return summedformula


def paralleldatabaseforms( featurelist, formppm, full, label, massbin, preferCHNO, wd, unlabeledcontrol):
    start = time.time()
    # Makes an output folder for all the formula calculations
    if unlabeledcontrol == 'No':
        output = 'Output'
    if unlabeledcontrol == 'Yes':
        output = 'Output/Unlabeled Output'
    try:
        os.mkdir('{}/{}/Peaks/Formulas'.format(wd, output))
    except FileExistsError:
        pass
    # Reads in the feature list and regression resulsts and merges the two on the unique rounded feature labels
    featlist = pd.read_csv('{}/{}/{}'.format(wd,output,str(featurelist)))
    results = pd.read_excel('{}/{}/TrimmedFeatureList.xlsx'.format(wd, output), index_col = 0)
    featlist['Unlabeled_Labeled Features'] = featlist['Feature']
    joined = pd.merge(results, featlist, on = 'Unlabeled_Labeled Features', how = 'left')
    joinedround = joined.copy()
    joinedround.drop_duplicates(inplace = True) 
        
    # Sorts the results by R2 and then drops any duplicate mass features (drops isomers with different retention times as identifying formula matches for each is unncessary)
    joinedround.sort_values(by = 'R2', inplace = True, ignore_index = True)
    joinedround.drop_duplicates(subset = ['MassFeature'], inplace = True)
    joinedround.reset_index(inplace = True)
    joinedroundfiltersort = joinedround.sort_values(by = 'R2', ignore_index = True)
    
    # Splits all results into individual dataframes for submitting to multiprocessing
    split = np.array_split(joinedroundfiltersort, len(joinedroundfiltersort))
    # Using partial function, locks in parameters for all database searching
    parallelformsetvars = partial(databaseformsearch, formppm = formppm, full = full, label = label, massbin = massbin, preferCHNO = preferCHNO, wd = wd)
    print('Generating Formulas:')
    # Submits all split result data for database searching using multiprocess, then concats all results into a single output and exports it to excel
    formresults = p_umap(parallelformsetvars, split)
    formresults = pd.concat(formresults)
    # formresults.to_excel('IntialFormulas.xlsx')
    # Combines formula results with regression results based on mass feature, will produce replicated results for isomers
    joinedform = pd.merge(joined, formresults, on = 'MassFeature', how = 'left')
    joinedform.to_excel('{}/{}/PeakInfo_RegressResults.xlsx'.format(wd,output))
    end = time.time()
    print("\n\nTook {} to determine formulas for {} compounds.".format(str(datetime.timedelta(seconds = (end-start))).split(".")[0], str(len(joinedround))))

def c13c12regressparallel(wd,featurelist, masswindow, rtwindow, minR2, removeplots, label, featuretype, unlabeledcontrol, minint):
    # Creates Output folders for regression outputs, redundant checking for some folders but minimal calculation time
    joinedoutfull = pd.DataFrame()
    if unlabeledcontrol == 'No':
        datapath = '/Data/Raw Files/'
        pklpath = '/Data/'
        output = 'Output'
    if unlabeledcontrol == 'Yes':
        datapath = '/Data/Unlabeled Control/Raw Files/'
        pklpath = '/Data/Unlabeled Control/'
        output = 'Output/Unlabeled Output'
    try:
        os.mkdir('{}/{}'.format(wd, output))
    except FileExistsError:
        pass
    try:
        os.mkdir('{}/{}/Peaks'.format(wd,output))
    except FileExistsError:
        pass    
    try:
        os.mkdir('{}/{}/Peaks/Plots'.format(wd, output))
    except FileExistsError:
        pass
    try:
        os.mkdir('{}/{}/Peaks/JoinedCSV'.format(wd,output))
    except FileExistsError:
        pass
    
    start = time.time()
    # Loads feature list from feature finding steps or premade
    if featuretype == 'Work From Premade Feature List':
        try:
            masslist = pd.read_csv('{}/{}.csv'.format(wd,str(featurelist)))
        except FileNotFoundError:
            pass
        try:
            masslist = pd.read_csv('{}/{}/{}.csv'.format(wd,output,str(featurelist)))
        except FileNotFoundError:
            pass
    else:
        masslist = pd.read_csv('{}/{}/{}.csv'.format(wd,output,str(featurelist)))
    # Dataframe to collect regression output
    allfilestats = pd.DataFrame()
    # Generates indexed list of all features to create multiindex during regression steps as well as provides a count of compounds
    c12c13 = list(masslist['Feature'])
    countcmpds = len(c12c13)
    countfiles = 0
    # Regular expression to identify pickle files generated in the last steps
    pklregex = re.compile('(.*)(.pkl)\Z')
    # Search through data folders and count pickled files 
    for files in os.listdir('{}{}'.format(wd,pklpath)):
        match = pklregex.search(files)
        if match == None:
            continue
        else:
            countfiles = countfiles + 1
    # Splits masslist into chunks equal to the number of cpus for splitting for multiprocessing
    masslistsplit = np.array_split(masslist, mp.cpu_count())   
    # Resets the index for each masslist so that iteration works
    for items in masslistsplit:
        items.reset_index(inplace = True)    
    # Sets all the regression parameters using the partial function  
    regressfilesparallelpartial = partial(regressfilesparallelloc, masswindow = masswindow, rtwindow = rtwindow, minR2 = minR2, label = label, wd = wd, unlabeledcontrol = unlabeledcontrol)     
    print('Regression on {} compounds in {} files'.format(str(countcmpds), str(countfiles)))
    # Runs regression sending each chunk of the masslist to separate CPUs, then concats all results into a single dataframe and fills NA with zeros
    filestats = p_umap(regressfilesparallelpartial, masslistsplit)   
    # joins = []
    # filestats = []
    # for items in out:
    #     joins.append(items['joined'])
    #     filestats.append(items['filestats'])
    # joinedoutfull = pd.concat(joins)  
    allfilestats = pd.concat(filestats)
    dataexportdf = allfilestats.fillna(0)
    # Outputs all regression data
    dataexportdf.to_excel('{}/{}/FullRegressionDataOutput.xlsx'.format(wd,output))
    # Groupby to find mean and standard deviation between files for each feature
    mean = allfilestats.groupby('Unlabeled_Labeled Features').mean(numeric_only = True)
    mean.sort_values(by = 'Unlabeled_Labeled Features', inplace = True)
    std = allfilestats.groupby('Unlabeled_Labeled Features').std(numeric_only = True)
    std.sort_values(by = 'Unlabeled_Labeled Features', inplace = True)
    stdmeanexport = pd.concat([mean,std], keys = ['Mean','Std'])
    stdmeanexport.to_excel('{}/{}/MeanStandardDevOutput.xlsx'.format(wd,output))
    # Trims the output datafiles based upon the minimum R2 provided by the user and then outputs all the data 
    meantrim = mean[mean['R2'] >= minR2]
    meanresindex = mean.reset_index()
    dffulltrim = meanresindex[meanresindex['Unlabeled_Labeled Features'].isin(meantrim.index.tolist())]
    dffulltrim.to_excel('{}/{}/TrimmedFeatureList.xlsx'.format(wd,output))
    # joinedoutfull.sort_index(inplace = True)
    # joinedoutfinal = joinedoutfull.groupby(joinedoutfull.index).sum()
    # joinedoutfinal.to_csv('{}/{}/JoinedEICS.csv'.format(wd,output))
    # Removes any plots which fall below the mininimum R2, determined by the user 
    
    if removeplots == 'Remove':
        means = mean.index.to_series()
        trims = meantrim.index.to_series()
        sub = means[~means.isin(trims)]
        for features in sub:
            shutil.rmtree('{}/{}/Peaks/Plots/{}'.format(wd,output,str(features)))
    end = time.time()
    print("\n\nTook {} for regression on {} compounds in {} files.".format(str(datetime.timedelta(seconds = (end-start))).split(".")[0], str(int(countcmpds)), str(countfiles)))

def regressfilesparallelloc(masslist, masswindow, rtwindow, minR2, label, wd, unlabeledcontrol):
    # Dataframe to collect all regression data from this masslist chunk 
    joinedout = pd.DataFrame()
    summedfilestatsdf = pd.DataFrame()
    # Regular expression to identify pickled datafiles
    pklregex = re.compile('(.*)(.pkl)\Z')
    
    
    # Iterates through data folder and identifies each pickled datafile
    if unlabeledcontrol == 'No':
        pklpath = '/Data/'
        output = 'Output'
    if unlabeledcontrol == 'Yes':
        pklpath = '/Data/Unlabeled Control/'
        output = 'Output/Unlabeled Output'
    for files in os.listdir('{}{}'.format(wd,pklpath)):
        match = pklregex.search(files)
        if match == None:
            continue
        else:

            # Loads in pickled dataset
            dataset = pd.read_pickle('{}{}{}.pkl'.format(wd,pklpath,str(match.group(1))))
            # Iterates through each line of the masslist to get mass and retention time for unlabeled and labeled mass as well as the unique feature pair name
            for tup in masslist.itertuples():        
                
                rt12 = tup.rtUnlabeled
                mz12 = tup.mzUnlabeled
                rt13 = tup.rtLabeled
                mz13 = tup.mzLabeled

                featpair = tup.Feature
                # Trys to make the folder to collect plots unless it already exists
                try:
                    os.mkdir('{}/{}/Peaks/Plots/{}'.format(wd,output,featpair))
                except FileExistsError:
                    pass
                # Using the user provided retention time window, calculates the max and minimum retention time for the labeled and unlabeled pairs
                C12RTMin = rt12 - rtwindow
                C12RTMax = rt12 + rtwindow
                C13RTMin = rt13 - rtwindow
                C13RTMax = rt13 + rtwindow
                # Filters full picked dataset to include only the masses of interest within the user provided mass window for the labeled and unlabeled datasets, additionally sums intensity of any mass within the given mass window
                C12scaninfoinitial = dataset.loc[(dataset['m/z'].values >= (float(mz12) - float(masswindow))) & (dataset['m/z'].values <= (float(mz12) + float(masswindow)))].sort_values(by = ['scan', 'int'], ascending = [True,False])
                C12scaninfo = pd.merge(C12scaninfoinitial[['scan','rt','m/z','polarity']], C12scaninfoinitial[['scan','int']].groupby(by = 'scan').sum(), on = 'scan', how = 'right').drop_duplicates(subset = 'scan', keep = 'first')
                C13scaninfoinitial = dataset.loc[(dataset['m/z'].values >= (float(mz13) - float(masswindow))) & (dataset['m/z'].values <= (float(mz13) + float(masswindow)))].sort_values(by = ['scan', 'int'], ascending = [True,False])
                C13scaninfo = pd.merge(C13scaninfoinitial[['scan','rt','m/z','polarity']], C13scaninfoinitial[['scan','int']].groupby(by = 'scan').sum(), on = 'scan', how = 'right').drop_duplicates(subset = 'scan', keep = 'first')
                # pd.concat([C12scaninfo, C13scaninfo], axis=1).to_excel('.//CalibrationTest/EICS{}_{}.xlsx'.format(str(featpair),files.strip('.pkl')))
                # Trims the extracted data within the given retention time range for the labeled and unlabeled compounds, if no data was extracted 0s are substituted for output data
                try:
                   C12scaninfotrim1 = C12scaninfo.loc[(C12scaninfo['rt'].values >= float(C12RTMin)) & (C12scaninfo['rt'].values <= float(C12RTMax))]
                except KeyError:
                    C12scaninfotrim1 = pd.DataFrame({'rt':[0], 'scan':[0], 'int':[0]})
                    C12scaninfo = pd.DataFrame({'rt':[0], 'scan':[0], 'int':[0]})
                try:
                    C13scaninfotrim1 = C13scaninfo.loc[(C13scaninfo['rt'].values >= float(C13RTMin)) & (C13scaninfo['rt'].values <= float(C13RTMax))]
                except KeyError:
                    C13scaninfotrim1 = pd.DataFrame({'rt':[0], 'scan':[0], 'int':[0]})
                    C13scaninfo = pd.DataFrame({'rt':[0], 'scan':[0], 'int':[0]})
                
                # Sorts trimmed data by intensity to find the scan with maximum intentsity within range
                C12scaninfo_sorted = C12scaninfotrim1.sort_values('int', ascending = False)
                # print(C12scaninfo_sorted, C12scaninfo_sorted['scan'].values[:int(0)])
                C12top = C12scaninfo_sorted['scan'].values[:int(1)]
                C13scaninfo_sorted = C13scaninfotrim1.sort_values('int', ascending = False)
                C13top = C13scaninfo_sorted['scan'].values[:int(1)]
                # Relic of older methods, used to use an average of multiple top scans to find the top, for now just an "average" of a single value
                try:
                    C12avg = round(np.average(C12top))
                except ValueError:
                    C12avg = 0
                try:
                    C13avg = round(np.average(C13top))
                except ValueError:
                    C13avg = 0
                # If the label is deuterium, finds the difference in scans between the 1H and 2H peaks which is then used to adjust the scan numbers in the extracted data to align the peaks for regression (Includes the chance that the deuterated peak elutes after the unlabled peak, should not happen but includes it just to cover exceptions) 
                if label == '2H':
                    if C13avg > C12avg:
                        avgadj = C13avg - C12avg
                    elif C12avg > C13avg:
                        avgadj = C12avg - C13avg
                    if C13avg > C12avg: 
                        C13scaninfo ['scan'] = C13scaninfo ['scan'] - float(avgadj)
                        adjmax = C13avg - float(avgadj)
                    elif C12avg > C13avg:
                        C13scaninfo ['scan'] = C13scaninfo ['scan'] + float(avgadj)
                        adjmax = C13avg + float(avgadj)
                else: 
                    adjmax = 0
                
                # Uses findidealscanwidth function (see below) to optimize the scan window to extract the data for regression, then uses that optimized scan width to exatract the exact data used for the regression
                try:
                    bestscanwindow, results, joined = findidealscanwidth(C12scaninfo, C13scaninfo, C12avg, C13avg, label, adjmax, featpair)
                    
                    C12scaninfo = C12scaninfo.loc[(C12scaninfo['scan'].values >= (float(C12avg) - float(bestscanwindow))) & (C12scaninfo['scan'].values <= (float(C12avg) + float(bestscanwindow)))]
                # Catches an error where the dataframe is empty, sets the ideal scan window to avoid later error
                except KeyError:
                    bestscanwindow = 0
                # Extractes labeled data depending on the label type
                if label == '2H':
                    C13scaninfo = C13scaninfo.loc[(C13scaninfo['scan'].values >= (float(adjmax) - float(bestscanwindow))) & (C13scaninfo['scan'].values <= (float(adjmax) + float(bestscanwindow)))]
                else:
                    C13scaninfo = C13scaninfo.loc[(C13scaninfo['scan'].values >= (float(C13avg) - float(bestscanwindow))) & (C13scaninfo['scan'].values <= (float(C13avg) + float(bestscanwindow)))]
                
                fullEIC = pd.concat([C12scaninfo.add_suffix('_Unlabeled'), C13scaninfo.add_suffix('_Labeled')], axis=1)
                fullEIC.to_csv('{}/{}/Peaks/JoinedCSV/{}_{}FullEIC.csv'.format(wd,output,featpair, str(match.group(1))))
                # Joins the labeled and unlabeled data based on scan number, joins based on intersetion of the two datasets 
                joined = pd.merge(C12scaninfo, C13scaninfo, suffixes=("_Unlabeled","_Labeled"), on = 'scan', how = 'inner')
                # Outputs joined data to csv for later inspection
                joined.to_csv('{}/{}/Peaks/JoinedCSV/{}_{}joined.csv'.format(wd,output,featpair, str(match.group(1))))
                
                                                      
                
                
                
                
                
                # Resphapes the unlabeled (x) data to align with the labaled data (y) for regression
                x = np.array(joined['int_Unlabeled']).reshape(-1,1)
                y = np.array(joined['int_Labeled'])
                x_int = np.array(joined['int_Unlabeled'])
                sumint = np.add(np.array(joined['int_Unlabeled']), np.array(joined['int_Labeled'])).reshape(-1,1)
                y_12peakint = np.array(C12scaninfo['int'])
                x_12peakrt = np.array(C12scaninfo['rt']).reshape(-1,1)
                x_13peakrt = np.array(C13scaninfo['rt']).reshape(-1,1)
                y_13peakint = np.array(C13scaninfo['int'])
                # Pulls some stats for further validation of regression
                # std_mean = np.mean(x)
                # std_median = np.median(x)
                # sil_mean = np.mean(y)
                # sil_median = np.median(y)
                # Extacts the peak height for the labeled and unlabeled as well as the number of points used for regression (based on the length of the joined dataframe)
                C12Height = np.amax(C12scaninfo['int'])
                C13Height = np.amax(C13scaninfo['int'])
                points = len(joined)
                # Calculates the peak height ratio or returns zero if unlabled height is zero
                try:
                    HeightRatio = C13Height/C12Height
                except ZeroDivisionError:
                    HeightRatio = 0
                # Regression Steps, any errors in these steps sets all output to zero and deletes the 
                try:
                    # Creates regression model, fits the data to it and extracts the slope, R2, and intercept
                    LR = LinearRegression()
                    LR.fit(x,y)
                    slope = LR.coef_[0]
                    r2 = LR.score(x,y)
                    intercept = LR.intercept_
                    # Removes the top 20,40,60% of points from the joined un/labeled data and does regression on these sets, if the peak overlap/regression is good then removing a percentage of points will have little effect on the R2
                    # xlog = np.log(x_int).reshape(-1,1)
                    # ylog = np.log(y)
                    # LRlog = LinearRegression()
                    # LRlog.fit(xlog,ylog)
                    # slopelog = LRlog.coef_[0]
                    # r2log = LRlog.score(xlog,ylog)
                    
                    
                    
                    # sr2 = stats.spearmanr(x_int,y).statistic
                    # spvalue = stats.spearmanr(x_int,y).pvalue
                    joinedsort = joined.sort_values(by = ['int_Unlabeled'], ascending = True)
                    length = len(joinedsort)
                    L1 = math.floor(length*0.2)
                    # L2 = math.floor(length*0.4)
                    # L3 = math.floor(length*0.6)
                    joineddrop1 = joinedsort.iloc[:-L1]
                    # joineddrop2 = joinedsort.iloc[:-L2]
                    # joineddrop3 = joinedsort.iloc[:-L3]
                    x1 = np.array(joineddrop1['int_Unlabeled']).reshape(-1,1)
                    # x1_int = np.array(joineddrop1['int_Unlabeled'])
                    y1 = np.array(joineddrop1['int_Labeled'])
                    # x2 = np.array(joineddrop2['int_Unlabeled']).reshape(-1,1)
                    # y2 = np.array(joineddrop2['int_Labeled'])
                    # x3 = np.array(joineddrop3['int_Unlabeled']).reshape(-1,1)
                    # y3 = np.array(joineddrop3['int_Labeled'])
                    LRLabelIncorp = LinearRegression()
                    LRLabelIncorp.fit(sumint, y)
                    slopelabelincorp = (LRLabelIncorp.coef_[0])*100
                    try:
                        LR1 = LinearRegression()
                        LR1.fit(x1,y1)
                        r21 = LR1.score(x1,y1)
                        
                    except ValueError:
                        r21 = 0
                        
                    
                    # Sets R2 equal to zero if slope is negative 
                    if slope < 0:
                        r2 = 0
                    # Sets R2 and slope = 0 if regression is on only 1 or two points
                    if points < 3: 
                        slope = 0
                        r2 = 0
                    # If any of the conditons are met, the exported csv is removed as the regression is "not good"
                    avgr2 = (r2 + r21)/2
                    
                    
                    
                    if avgr2 < minR2 or r2 < minR2 or r21 < minR2 or slope <= 0 or points < 3:
                        os.remove('{}/{}/Peaks/JoinedCSV/{}_{}joined.csv'.format(wd,output,featpair, str(match.group(1))))
                        os.remove('{}/{}/Peaks/JoinedCSV/{}_{}FullEIC.csv'.format(wd,output,featpair, str(match.group(1))))
                        continue 
                    
                    # Using matplot plots the regression and exports to the plot folder
                    plt.scatter(x,y)
                    plt.plot(x, LR.predict(x))
                    axregress = plt.gca()
                    plt.title('File: {} Feat: {}'.format(files.strip('.mzML'), featpair))
                    plt.xlabel('Unlabeled Abundance')
                    plt.ylabel('Labeled Abundance')
                    plt.text(0.75, 0, 'Slope = {} \n R^2 = {}'.format(str(round(LR.coef_[0], 3)),str(round(LR.score(x,y),3))), transform = axregress.transAxes)
                    plt.savefig('{}/{}/Peaks/Plots/{}/Regression{} {}.png'.format(wd,output,featpair, files.strip('.pkl'), featpair), bbox_inches = 'tight')
                    plt.close()
                    
                    fig, ax1 = plt.subplots()                   
                    axtwin = ax1.twinx()
                    ax1.plot(x_12peakrt,y_12peakint, color = 'b', label = 'unlabeled')
                    axtwin.plot(x_13peakrt,y_13peakint, color = 'r', label = 'labeled')
                    lines, labels = ax1.get_legend_handles_labels()
                    lines2, labels2 = axtwin.get_legend_handles_labels()
                    axtwin.legend(lines + lines2, labels + labels2, loc=0)
                    fig.suptitle('File: {} Feat: {}'.format(files.strip('.mzML'), featpair))
                    ax1.set_xlabel('Retention Time')
                    ax1.set_ylabel('Unlabeled Intensity')
                    axtwin.set_ylabel('Labeled Intensity')
                    
                    fig.savefig('{}/Output/Peaks/Plots/{}/Peak{} {}.png'.format(wd,featpair, files.strip('.pkl'), featpair), bbox_inches='tight')
                    plt.close(fig)
                except (ValueError):
                    slope = 0
                    r2 = 0
                    intercept = 0
                    pr2 = 0
                    pr2over = 0 
                   
                    continue
                try:
                    
                    C12intsortnoshape = C12scaninfo['int']
                    C13intsortnoshape = C13scaninfo['int']
                    C12calibrationpeak = pd.DataFrame(peaksimulation(np.max(C12intsortnoshape), np.std(C12intsortnoshape), len(C12intsortnoshape)), columns = ['int'])
                    C13calibrationpeak = pd.DataFrame(peaksimulation(np.max(C13intsortnoshape), np.std(C13intsortnoshape), len(C13intsortnoshape)), columns = ['int'])
                    C12intsortnoshape = pd.DataFrame(C12scaninfo['int'], columns = ['int'])
                    C13intsortnoshape = pd.DataFrame(C13scaninfo['int'], columns = ['int'])
                    
                    C12Merged = alignsimpeak(C12intsortnoshape, C12calibrationpeak)
                    C13Merged = alignsimpeak(C13intsortnoshape, C13calibrationpeak)
                    
                    C12intsortnoshape = np.array(C12Merged['int_ex'])
                    C13intsortnoshape = np.array(C13Merged['int_ex'])
                    C12calibrationpeak = np.array(C12Merged['int_sim'])
                    C13calibrationpeak = np.array(C13Merged['int_sim'])
                    
                    C12CalPearR2 = stats.pearsonr(C12intsortnoshape,C12calibrationpeak).statistic
                    C13CalPearR2 = stats.pearsonr(C13intsortnoshape,C13calibrationpeak).statistic
                    pd.concat([C12Merged, C13Merged], axis=1).to_csv('{}/{}/Peaks/JoinedCSV/{}_{}GaussianFit.csv'.format(wd,output,featpair, str(match.group(1))))
                    
                    SimulatedPeakAvgScore = fmean([C12CalPearR2, C13CalPearR2])
                    SimulatedPeakVarScore = np.var([C12CalPearR2, C13CalPearR2])
                    CompleteScore = fmean([r2, r21, C12CalPearR2, C13CalPearR2]) 
                    CompleteScorewVar = fmean([r2, r21, C12CalPearR2, C13CalPearR2]) * (1 - abs(SimulatedPeakVarScore))
                except ValueError:
                    
                    C12CalPearR2 = 0
                    C13CalPearR2 = 0
                    SimulatedPeakAvgScore = 0
                    SimulatedPeakVarScore = 0
                    CompleteScore = 0 
                    CompleteScorewVar = 0
                if len(y_12peakint) == len(y_13peakint):
                    pr2 = stats.pearsonr(y_12peakint,y_13peakint).statistic
                else:
                    pr2 = 0
                pr2over = stats.pearsonr(x_int, y).statistic
                
                filestats = { 'Unlabeled_Labeled Features':featpair, 'DataFile':[files.strip('.pkl')],'Y-Intercept':[intercept], 'Slope':[slope], 'R2':[r2], 'UnlabeledHeight':[C12Height], 'LabeledHeight':[C13Height],'HeightRatio':[HeightRatio],'20% R2':[r21], 'RegressPoints':[points], 'OptimizedScanWindow':[bestscanwindow], 'SlopeDerivedLabelIncorp':[slopelabelincorp], '20%R2_R2 Avg':[avgr2],'C12CalPearR2':[C12CalPearR2],'C13CalPearR2':[C13CalPearR2], 'SimPeakAvg':[SimulatedPeakAvgScore], 'CompleteScore':[CompleteScore]}
                filestatsdf = pd.DataFrame.from_dict(filestats, orient = 'columns')
                summedfilestatsdf = pd.concat([summedfilestatsdf, filestatsdf], axis = 0)
                
    return summedfilestatsdf

def findidealscanwidth(C12scaninfo, C13scaninfo, C12avg, C13avg, label, adjmax, featpair):
    # Set of scan ranges to test to find ideal scan range
    scanranges = [10,12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50]
    r2results = []
    scanrangeout = []
    joinlen = []
    # Iterates through scan ranges to find the best R2
    for scanwindow in scanranges:
        # Extracts all data based on the previously identified peak height and the current scan window
        C12scaninfo2 = C12scaninfo.loc[(C12scaninfo['scan'].values >= (float(C12avg) - float(scanwindow))) & (C12scaninfo['scan'].values <= (float(C12avg) + float(scanwindow)))]
        if label == '2H':
            C13scaninfo2 = C13scaninfo.loc[(C13scaninfo['scan'].values >= (float(adjmax) - float(scanwindow))) & (C13scaninfo['scan'].values <= (float(adjmax) + float(scanwindow)))]
        else:
            C13scaninfo2 = C13scaninfo.loc[(C13scaninfo['scan'].values >= (float(C13avg) - float(scanwindow))) & (C13scaninfo['scan'].values <= (float(C13avg) + float(scanwindow)))]    
        # Merges data based on scan number, on interesecting scans, then arranges the data for regression
        joined = pd.merge(C12scaninfo2, C13scaninfo2, suffixes=("_Unlabeled","_Labeled"), on = 'scan', how = 'inner')
        joinedlength = len(joined)
        x = np.array(joined['int_Unlabeled']).reshape(-1,1)
        y = np.array(joined['int_Labeled'])
        # Performs regression on the current scan range data and finds the R2 and appends to output, except for an empty dataframe which gives an r2 of 0
        try:
            LRtest = LinearRegression()
            LRtest.fit(x,y)
            r2test = LRtest.score(x,y)
            joinedsort = joined.sort_values(by = ['int_Unlabeled'], ascending = True)
            length = len(joinedsort)
            L1 = math.floor(length*0.2)
            joineddrop1 = joinedsort.iloc[:-L1]

            x1 = np.array(joineddrop1['int_Unlabeled']).reshape(-1,1)
            y1 = np.array(joineddrop1['int_Labeled'])
            try:
                LR1 = LinearRegression()
                LR1.fit(x1,y1)
                r21 = LR1.score(x1,y1)
            except ValueError:
                r21 = 0
            r2testout = mean([r2test,r21])
        except ValueError:
            r2testout = 0
            
        r2results.append(r2testout)
        scanrangeout.append(scanwindow)
        joinlen.append(joinedlength)
    # Sorts the R2 results and returns the best scan range
    results = pd.DataFrame({'Scans':scanrangeout, 'R2':r2results, 'len':joinlen})
    results['R2Round'] = results['R2'].round(1)
    # results.to_excel('.//Test Lengths/Test_{}.xlsx'.format(featpair))
    results.sort_values(by = ['R2Round','Scans'], ascending = [False, False], ignore_index = True, inplace = True)
    results = results.drop('R2Round', axis = 1)
    # results.sort_values(by = ['R2'], ascending = False, ignore_index = True, inplace = True)
    bestwindow = results['Scans'][0]
    
    return bestwindow, results, joined


def peaksimulation (peakmax, peakstd, length):
    var = math.sqrt(peakstd)
    # var = peakstd/2
    dist = np.linspace(peakmax - 3*var, peakmax + 3*var, length)
    
    output = stats.norm.pdf(dist,peakmax,var)
    outputadj = (output/np.max(output))*peakmax
    return outputadj
    
def alignsimpeak(exdata, simdata):
    exdata.reset_index(inplace = True)
    simdata.reset_index(inplace = True)
    exdata['Index'] = exdata.index
    simdata['Index'] = simdata.index
    exdatamax = exdata['int'].idxmax()
    simdatamax = simdata['int'].idxmax()
    if exdatamax > simdatamax:
        adj = exdatamax - simdatamax
        simdata['Index'] = simdata['Index'] + adj
    elif simdatamax > exdatamax:
        adj = simdatamax - exdatamax 
        simdata['Index'] = simdata['Index'] - adj
    merged = pd.merge(exdata,simdata, how = 'right', on = 'Index', suffixes = ['_ex','_sim']).fillna(0)
    return merged

def roundrt(rt, rtbin):
    if rtbin == 0:
        return rt
    else:
        return(round(round(rt / rtbin) * rtbin,2))



# -*- coding: utf-8 -*-
"""
Created on Thu May 11 17:14:54 2023

@author: Evan Larson, Hegeman & Cohen Labs - University of Minnesota
"""
from functools import partial
import pandas as pd, os, re, time, datetime, sys
from pyteomics import mzml
from p_tqdm import p_umap
import numpy as np
import multiprocessing as mp


# Function takes each provided datafile in the data folders, and extracts all labeled feature pairs using the provided user parameters 
def findfeaturesbyscanparallel(wd,mincarbons, maxcarbons, masswindow, minint, label, rtbin, massbin, minscan, unlabeledcontrol):
    start = time.time()
    
    # Makes Output Directory or passes if it already exists
    try:
        os.mkdir('{}/Output'.format(wd))
    except FileExistsError:
        pass
    if unlabeledcontrol == 'Yes':
        try:
            os.mkdir('{}/Output/Unlabeled Output'.format(wd))
        except FileExistsError:
            pass
    # Regular Expression to identify mzml files within the data folder and regular expression to extract scan number within loaded mzml files
    mzmlregex = re.compile('(.*)(.mzML)\Z')
    scanregex = re.compile('scan=(\d*)')
    # Dataframe to collect final results
    results = pd.DataFrame()
    # Counts number for files in data folder 
    filenum = 0

    try:
        labelmd = float(label)
        userlabel = True
    except ValueError:
        userlabel = False
        if not label in {'13C', '15N', '2H', '18O', '34S'}:
            print('Incorrect Label Name, Make sure to use 13C, 15N, 18O, 34S or provide your own mass difference')
            sys.exit()
    if unlabeledcontrol == 'No':
        datapath = '/Data/Raw Files'
        pklpath = '/Data/'
        output = 'Output'
    if unlabeledcontrol == 'Yes':
        datapath = '/Data/Unlabeled Control/Raw Files'
        pklpath = '/Data/Unlabeled Control/'
        output = 'Output/Unlabeled Output'
    for files in os.listdir('{}{}'.format(wd, datapath)):
        filenum = filenum + 1 
        # Regular Expression search to identify MZML files within Datafolder
        match = mzmlregex.search(files)
        
        if match == None:
            continue
        else:
            # Reads mzml file into program using pyteomics module
            if unlabeledcontrol == 'No':
                data = mzml.read('{}/Data/Raw Files/{}'.format(wd,str(files)))
            if unlabeledcontrol == 'Yes':
                data = mzml.read('{}/Data/Unlabeled Control/Raw Files/{}'.format(wd,str(files)))
            spectra = 0
            # List to collect all spectra
            spectralist = []
            # Within loaded mzml files only selects mass spectra at ms1 level
            for spectrum in data:
                if not "m/z array" in spectrum:
                    continue
                elif (spectrum["ms level"] == 1):
                    spectra = spectra + 1
                    # For each spectrum collects mass, intensity, retention time, scan number and polarity into a dataframe. Polarity, retention time and scan number are the same for each scan
                    mzdata = list(spectrum["m/z array"])
                    intensitydata = list(spectrum["intensity array"])
                    rt = spectrum['scanList']['scan'][0]['scan start time']
                    scan = re.search(scanregex, spectrum['id']).group(1).replace('scan=','')
                    rtdata = len(mzdata) * [rt]
                    if 'negative scan' in spectrum:
                        pol = 'NEGATIVE'
                    if 'positive scan' in spectrum:
                        pol = 'POSITIVE'
                    polarity = [pol] * len(mzdata)
                    specdata = {'scan':int(scan), 'm/z':mzdata, 'int':intensitydata, 'rt':rtdata, 'polarity':polarity}
                    specdfinit = pd.DataFrame.from_dict(specdata)
                    # filters out masses that fall below the minimum intensity provided by the user. Index is reset for later indexing steps since low intensity things are removed
                    specdf = specdfinit[specdfinit['int'] > minint]
                    specdf.reset_index(inplace = True, drop = True)
                    # appends this spectrum to the collecting list
                    spectralist.append(specdf)
            # Concats list of spectra dataframes into a single dataframe and then outputs to a pickle format for easy access for later regression steps
            pd.concat(spectralist).to_pickle('{}{}{}.pkl'.format(wd,pklpath, match.group(1)))
            print('Finding Features in {} scans for {}'.format(str(spectra), str(match.group(1))))
            # Chooses between using feature extraction using algorithm for deuterated or undeuterated data
            if label == '2H':
                # Splits list of spectra for deuterated algorithm into overlapping chunks with sufficient length to split evenly across all cpu cores, overlap is 120. 
                length = int((len(spectralist)+(120*mp.cpu_count()))/mp.cpu_count())
                spectralistDeuterated = [spectralist[i:i+length] for i in range(0,len(spectralist), length - 120)]
                # spectralistDeuterated = [spectralist[i : i + 240] for i in range(0, len(spectralist), 120)]
                # Uses partial function to set unchanged parameters for feature finding steps for deuterated data
                findfeaturesparallelpartdeuterated = partial(findfeaturesDeuteratednew, mincarbons = mincarbons, maxcarbons = maxcarbons, minint = minint, masswindow = masswindow, massbin = massbin, rtbin = rtbin)
                # Submits overlapped chunked spectra lists of deuterated data to multiprocessing function using find deuterated features algorithm (see below)
                specresults = p_umap(findfeaturesparallelpartdeuterated, spectralistDeuterated)
                # Concats output from feature finding to a single dataframe containing all identified features
                specresults = pd.concat(specresults)
                
            else:
                # Uses partial function to set unchanged parameters for feature finding steps
                if userlabel == False:
                    findfeaturesparallelpartlocnew = partial(findfeaturesparallellocnew, mincarbons = mincarbons, maxcarbons = maxcarbons, minint = minint, label = label, masswindow = masswindow)
                if userlabel == True:
                    findfeaturesparallelpartlocnew = partial(findfeaturesparallellocnew, mincarbons = mincarbons, maxcarbons = maxcarbons, minint = minint, label = labelmd, masswindow = masswindow)
                # Submits list of individual spectral dataframes to multiprocessing function using find features algorithm (see below)
                specresults = p_umap(findfeaturesparallelpartlocnew, spectralist)
                # Concats output from feature finding to a single dataframe containing all identified features
                specresults = pd.concat(specresults)
            
            # Outputs all identified features for current datafile
            # specresults.to_csv('{}/{}/AllExtractedFeatures_{}.csv'.format(wd, output, str(files).replace('.mzML','')))
            # Partial function locks the bin for retention time and mass rounding 
            rtroundpartial = partial(roundrt, rtbin = rtbin)
            massbinpartial = partial(roundmass, massbin = massbin)
            # Adds new columns based on rounded masses and retention times to group features within those bins
            specresults['m/zUnlabeledround'] = specresults['mzUnlabeled'].apply(massbinpartial)
            specresults['m/zLabeledround'] = specresults['mzLabeled'].apply(massbinpartial)
            specresults['rtUnlabeledround'] = specresults['rtUnlabeled'].apply(rtroundpartial)
            specresults['rtLabeledround'] = specresults['rtLabeled'].apply(rtroundpartial)
            # specresults.loc[(specresults['mzUnlabeled'] > 124) & (specresults['mzUnlabeled'] < 125) & (specresults['rtUnlabeled'] < 2)].to_excel('test_beforeConsolidation_{}.xlsx'.format(str(files).replace('.mzML','')))
            # Groupby on rounded columns to count number of instances of each set which equates to the number of spectra where the feature sets are identified
            count = specresults.groupby(['rtUnlabeledround','rtLabeledround','m/zUnlabeledround', 'm/zLabeledround']).size()
            count = count.reset_index()                       
            count.rename({0:'spectral count'}, axis = 'columns', inplace = True)
            # Merges the spectral count based on retention time and mass for labeled and unlabeled feature sets
            mergedwithcount = pd.merge(specresults, count, on = ['rtUnlabeledround','rtLabeledround','m/zUnlabeledround', 'm/zLabeledround'])
            # Sorts full set by the intensity of the unlabled peak
            mergedwithcount.sort_values(by = 'm/zUnlabeled int', inplace = True, ascending = False)
            # Drops duplicates based on based on retention time and mass for labeled and unlabeled feature sets, leaving only the highest intensity 
            mergedwithcount.drop_duplicates(subset = ['rtUnlabeledround','rtLabeledround','m/zUnlabeledround', 'm/zLabeledround'], inplace = True)
            #  Resorts based on mass
            mergedwithcount.sort_values(by = 'mzUnlabeled', inplace = True)
            # Adds column to output based on file in which the feature was found
            mergedwithcount['file'] = (len(mergedwithcount) * [str(match.group(1))])
            # mergedwithcount.loc[(mergedwithcount['mzUnlabeled'] > 124) & (mergedwithcount['mzUnlabeled'] < 125) & (mergedwithcount['rtUnlabeled'] < 2)].to_excel('test_AfterConsolidation_{}.xlsx'.format(str(files).replace('.mzML','')))
            # Removes any features that present in fewer than the minimum number of scans provided by the user and any features that are present in more than 10% of all spectra, to try to remove background
            mergedfiltered = mergedwithcount.loc[(mergedwithcount['spectral count'] >= minscan) & (mergedwithcount['spectral count'] <= (0.1*spectra))]
            # Adds individual file results to overall results output
            results = pd.concat([results, mergedfiltered],axis = 0, ignore_index = True)
    
    # Creates new feature and mass feature column to label each feature identified
    results['Feature'] = results['m/zUnlabeledround'].astype(str) + '_' + results['rtUnlabeledround'].astype(str) + '_' + results['m/zLabeledround'].astype(str) + '_' + results['rtLabeledround'].astype(str)
    results['MassFeature'] = results['m/zUnlabeledround'].astype(str) + '_' + results['m/zLabeledround'].astype(str) 
    #  Groupby on the feature column to count the number of files each feature is present in which is then added as a file count column in the output folder
    # results.loc[(results['mzUnlabeled'] > 124) & (results['mzUnlabeled'] < 125) & (results['rtUnlabeled'] < 2)].to_excel('test_AfterConsolidationBeforeFileGrouping.xlsx')
    fileresults = results.groupby(['Feature']).size()
    fileresults = fileresults.reset_index()
    fileresults.rename({0:'filecount'}, axis = 'columns', inplace = True)
    fileresultscount = pd.merge(results, fileresults, on = ['Feature'])
    # Removes features based on the number of files those features are present, based on the number of files used to find features, present in greater than or equal to the number of files - 1 
    fileresultscount = fileresultscount.loc[(fileresultscount['filecount'] >= (filenum - 1))]
    # Results are sorted by unlabeled intensity and then duplicate features are dropped
    fileresultscount.sort_values(by = 'm/zUnlabeled int', inplace = True, ascending = False)
    fileresultscount.drop_duplicates(subset = ['Feature'], inplace = True)
    end = time.time()
    # Results are output as a feature list csv which are used in later steps
    fileresultscount.drop(labels = ['m/zUnlabeled int','m/zLabeled int', 'Labeled/Unlabeled Height', 'file'], axis = 1, inplace = True)
    fileresultscount.to_csv('{}/{}/GeneratedFeatureList.csv'.format(wd,output))
    print("\n\nTook {} to analyze {} files.".format(str(datetime.timedelta(seconds = (end-start))).split(".")[0], str(filenum)))
    

def findfeaturesparallellocnew(specdf, mincarbons, maxcarbons, minint, label, masswindow):   
    # Sets the standard mass difference between unlabeled and labeled peaks, all matched features will be a multiple of these
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
    # Dataframe to collect all results
    specresults = pd.DataFrame()
    # Resets index for iterating through the dataframe
    specdf.reset_index(inplace = True, drop = True)
    for rows1 in range(len(specdf)):
        # Sets mass, intensity, retention time and polarity variables for first mass
        mz1 = specdf.loc[rows1,'m/z']
        mz1int = specdf.loc[rows1,'int']
        rt = specdf['rt'][0]
        pol = specdf['polarity'][0]
        # if rt > 2:
        #     break
        # If inteisty is above the minimum intensity provided by the user continue (may be redundnadt based on previous filtering)
        if mz1int >= minint:
            # Filters values and finds matches based on a set of conditions: 
            # 1. Finds if the difference between the current mass and any other masses is a multiple of the label mass difference within a user given range 
            # 2. Checks that the polarity of the masses are the same (should not be a problem as it is scan by scan but check adds little computation time)
            # 3. Makes sure that potential matches are larger in mass than the current mass
            # 4. Finally confirms that the number of label multiple is within the range given by the user 
            spectramatch = specdf.loc[((((specdf['m/z'].values - mz1) % md) <= masswindow) | (((md - (specdf['m/z'].values - mz1)) % md) <= masswindow)) & (specdf['polarity'] == pol) & (mz1 < specdf['m/z'].values) & ((np.round((specdf['m/z'].values - mz1)/md) >= mincarbons) & (np.round((specdf['m/z'].values - mz1)/md) <= maxcarbons))]
            # Continues as long as there are any matches
            if len(spectramatch) > 0:
                # Resets index for iteration
                spectramatch.reset_index(inplace = True, drop = True)
                for matches in range(len(spectramatch)):
                    # Dataframe to collect matched compounds
                    addmatch = pd.DataFrame()
                    # Creates match dataframe which are filled with the mass and intensity values for each match
                    mz2 = spectramatch.loc[matches, 'm/z']
                    mz2int = spectramatch.loc[matches, 'int']
                    addmatch = pd.DataFrame({'mzUnlabeled':mz1, 'm/zUnlabeled int':mz1int, 'mzLabeled':mz2, 'm/zLabeled int':mz2int, 'rtUnlabeled':rt,'rtLabeled':rt, 'Labeled/Unlabeled Height':((mz2int)/mz1int),'pol':pol, 'NumLabels':round((mz2 - mz1)/md), 'exact':(mz2 - mz1)/md}, index = [0])
                    # Adds match dataframe to output dataframe
                    specresults = pd.concat([specresults, addmatch], axis = 0, ignore_index = True)
    return specresults

def findfeaturesparallellocnewfromunlabeledtarget(specdf, mincarbons, maxcarbons, minint, label, masswindow,unlabeledtarget,rtbin):   
    # Sets the standard mass difference between unlabeled and labeled peaks, all matched features will be a multiple of these
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
    # Dataframe to collect all results
    specresults = pd.DataFrame()
    # Resets index for iterating through the dataframe
    specdf.reset_index(inplace = True, drop = True)
    for rows1 in range(len(unlabeledtarget)):
        # Sets mass, intensity, retention time and polarity variables for first mass
        mz1 = unlabeledtarget[rows1,'m/z']
        rt = unlabeledtarget[rows1, 'rt']
        pol = unlabeledtarget[rows1, 'polarity']
        
        # If inteisty is above the minimum intensity provided by the user continue (may be redundnadt based on previous filtering)
        
        # Filters values and finds matches based on a set of conditions: 
        # 1. Finds if the difference between the current mass and any other masses is a multiple of the label mass difference within a user given range 
        # 2. Checks that the polarity of the masses are the same (should not be a problem as it is scan by scan but check adds little computation time)
        # 3. Makes sure that potential matches are larger in mass than the current mass
        # 4. Finally confirms that the number of label multiple is within the range given by the user 
        spectramatch = specdf.loc[((((specdf['m/z'].values - mz1) % md) <= masswindow) | (((md - (specdf['m/z'].values - mz1)) % md) <= masswindow)) & (specdf['polarity'] == pol) & (mz1 < specdf['m/z'].values) & ((np.round((specdf['m/z'].values - mz1)/md) >= mincarbons) & (np.round((specdf['m/z'].values - mz1)/md) <= maxcarbons)) & ((specdf['rt'] >= rt - rtbin) & (specdf['rt'] <= rt + rtbin)) ]
        # Continues as long as there are any matches
        if len(spectramatch) > 0:
            # Resets index for iteration
            spectramatch.reset_index(inplace = True, drop = True)
            for matches in range(len(spectramatch)):
                # Dataframe to collect matched compounds
                addmatch = pd.DataFrame()
                # Creates match dataframe which are filled with the mass and intensity values for each match
                mz2 = spectramatch.loc[matches, 'm/z']
                mz2int = spectramatch.loc[matches, 'int']
                addmatch = pd.DataFrame({'mzUnlabeled':mz1, 'mzLabeled':mz2, 'm/zLabeled int':mz2int, 'rtUnlabeled':rt,'rtLabeled':rt,'pol':pol, 'NumLabels':round((mz2 - mz1)/md), 'exact':(mz2 - mz1)/md}, index = [0])
                # Adds match dataframe to output dataframe
                specresults = pd.concat([specresults, addmatch], axis = 0, ignore_index = True)
    return specresults





def roundrt(rt, rtbin):
    return(round(round(rt / rtbin) * rtbin,1))

def roundmass(mass, massbin):
    return(round(round(mass / massbin) * massbin,3))

def MakePkls(wd, minint, unlabeledcontrol):
    
    # Makes Output Directory or passes if it already exists
    try:
        os.mkdir('{}/Output'.format(wd))
    except FileExistsError:
        pass
    if unlabeledcontrol == 'Yes':
        try:
            os.mkdir('{}/Output/Unlabeled Output'.format(wd))
        except FileExistsError:
            pass
    # Regular Expression to identify mzml files within the data folder and regular expression to extract scan number within loaded mzml files
    mzmlregex = re.compile('(.*)(.mzML)\Z')
    scanregex = re.compile('scan=(\d*)')
    # Dataframe to collect final results
    # Counts number for files in data folder 
    filenum = 0

    if unlabeledcontrol == 'No':
        datapath = '/Data/Raw Files'
        pklpath = '/Data/'
        output = 'Output'
    if unlabeledcontrol == 'Yes':
        datapath = '/Data/Unlabeled Control/Raw Files'
        pklpath = '/Data/Unlabeled Control/'
        output = 'Output/Unlabeled Output'
    for files in os.listdir('{}{}'.format(wd, datapath)):
        filenum = filenum + 1 
        # Regular Expression search to identify MZML files within Datafolder
        match = mzmlregex.search(files)
        
        if match == None:
            continue
        else:
            # Reads mzml file into program using pyteomics module
            if unlabeledcontrol == 'No':
                data = mzml.read('{}/Data/Raw Files/{}'.format(wd,str(files)))
            if unlabeledcontrol == 'Yes':
                data = mzml.read('{}/Data/Unlabeled Control/Raw Files/{}'.format(wd,str(files)))
            spectra = 0
            # List to collect all spectra
            spectralist = []
            # Within loaded mzml files only selects mass spectra at ms1 level
            for spectrum in data:
                if not "m/z array" in spectrum:
                    continue
                elif (spectrum["ms level"] == 1):
                    spectra = spectra + 1
                    # For each spectrum collects mass, intensity, retention time, scan number and polarity into a dataframe. Polarity, retention time and scan number are the same for each scan
                    mzdata = list(spectrum["m/z array"])
                    intensitydata = list(spectrum["intensity array"])
                    rt = spectrum['scanList']['scan'][0]['scan start time']
                    scan = re.search(scanregex, spectrum['id']).group(1).replace('scan=','')
                    rtdata = len(mzdata) * [rt]
                    if 'negative scan' in spectrum:
                        pol = 'NEGATIVE'
                    if 'positive scan' in spectrum:
                        pol = 'POSITIVE'
                    polarity = [pol] * len(mzdata)
                    specdata = {'scan':int(scan), 'm/z':mzdata, 'int':intensitydata, 'rt':rtdata, 'polarity':polarity}
                    specdfinit = pd.DataFrame.from_dict(specdata)
                    # filters out masses that fall below the minimum intensity provided by the user. Index is reset for later indexing steps since low intensity things are removed
                    specdf = specdfinit[specdfinit['int'] > minint]
                    specdf.reset_index(inplace = True, drop = True)
                    # appends this spectrum to the collecting list
                    spectralist.append(specdf)
            # Concats list of spectra dataframes into a single dataframe and then outputs to a pickle format for easy access for later regression steps
            pd.concat(spectralist).to_pickle('{}{}{}.pkl'.format(wd,pklpath, match.group(1)))


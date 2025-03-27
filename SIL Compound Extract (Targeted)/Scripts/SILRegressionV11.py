# -*- coding: utf-8 -*-
"""
Created on Fri Sep 30 13:01:23 2022

@author: Evan Larson, Hegeman & Cohen Labs - University of Minnesota
"""


import pandas as pd, os, re, logging, math
import warnings
from sklearn.exceptions import UndefinedMetricWarning
warnings.filterwarnings("ignore", category=UndefinedMetricWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)
from sklearn.linear_model import LinearRegression
import numpy as np

from matplotlib import pyplot as plt
from pyteomics import mzml

def makepkls():
    mzmlregex = re.compile('(.*)(.mzML)\Z')
    scanregex = re.compile('scan=(\d*)')
    for files in os.listdir('.//SIM'):
        match = mzmlregex.search(files)
        
        if match == None:
            continue
        else:
            data = mzml.read('.//SIM/{}'.format(str(files)))
            spectra = 0
            spectralist = []
            for spectrum in data:
                if not "m/z array" in spectrum:
                    continue
                elif (spectrum["ms level"] == 1):
                    spectra = spectra + 1
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
                    specdata = {'scan':int(scan), 'm/z':mzdata, 'i':intensitydata, 'rt':rtdata, 'polarity':polarity}
                    specdf = pd.DataFrame.from_dict(specdata)
                    specdf.reset_index(inplace = True, drop = True)
                    spectralist.append(specdf)
            pd.concat(spectralist).to_pickle('.//SIM/{}.pkl'.format(match.group(1)))
    for files in os.listdir('.//MS1'):
        match = mzmlregex.search(files)
        
        if match == None:
            continue
        else:
            data = mzml.read('.//MS1/{}'.format(str(files)))
            spectra = 0
            spectralist = []
            for spectrum in data:
                if not "m/z array" in spectrum:
                    continue
                elif (spectrum["ms level"] == 1):
                    spectra = spectra + 1
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
                    specdata = {'scan':int(scan), 'm/z':mzdata, 'i':intensitydata, 'rt':rtdata, 'polarity':polarity}
                    specdf = pd.DataFrame.from_dict(specdata)
                    specdf.reset_index(inplace = True, drop = True)
                    spectralist.append(specdf)
            pd.concat(spectralist).to_pickle('.//MS1/{}.pkl'.format(match.group(1)))
    for files in os.listdir('.//MS2'):
        match = mzmlregex.search(files)
        
        if match == None:
            continue
        else:
            data = mzml.read('.//MS2/{}'.format(str(files)))
            spectra = 0
            spectralist = []
            precursorlist = []
            for spectrum in data:
                if not "m/z array" in spectrum:
                    continue
                elif (spectrum["ms level"] == 2):
                    spectra = spectra + 1
                    mzdata = list(spectrum["m/z array"])
                    intensitydata = list(spectrum["intensity array"])
                    rt = spectrum['scanList']['scan'][0]['scan start time']
                    precursor = spectrum['precursorList']['precursor'][0]['isolationWindow']['isolation window target m/z']
                    scan = re.search(scanregex, spectrum['id']).group(1).replace('scan=','')
                    rtdata = len(mzdata) * [rt]
                    if 'negative scan' in spectrum:
                        pol = 'NEGATIVE'
                    if 'positive scan' in spectrum:
                        pol = 'POSITIVE'
                    polarity = [pol] * len(mzdata)
                    specdata = {'scan':int(scan), 'm/z':mzdata, 'i':intensitydata, 'rt':rtdata, 'polarity':polarity, 'precursor':round(float(precursor),4)}
                    precursorlist.append(str(round(float(precursor),4)))
                    specdf = pd.DataFrame.from_dict(specdata)
                    specdf.reset_index(inplace = True, drop = True)
                    spectralist.append(specdf)
            pd.concat(spectralist).to_pickle('.//MS2/{}.pkl'.format(match.group(1)))
            # precursors = pd.DataFrame({'Precursors':precursorlist})
            # precursors.drop_duplicates(inplace = True)
            # precursors.to_excel('.//Output/Precursors_{}.xlsx'.format(str(files)))
def removepkls():
    pklregex = re.compile('(.*)(.pkl)')
    for types in ['MS1', 'MS2', 'SIM']:
        for files in os.listdir('.//{}'.format(str(types))):
            match = pklregex.search(files)
            if match == None:
                continue
            else:
                os.remove('.//{}/{}'.format(str(types),str(files)))

def SILCmpdExtract (masslistlink = 'Masslist', optimize = False, replicates = True, quant = 'n', perincorp = 'n'):        
    

    
    count = 0 
    if optimize == False:
        masslist = pd.read_excel(masslistlink)
    else:
        masslist = pd.read_excel(masslistlink)
    dataexportdf = pd.DataFrame()
    fulldatasetdf = pd.DataFrame()
    try:
        os.mkdir('.//Output')
    except FileExistsError:
        pass
    logging.basicConfig(filename = './Output/ErrorLog.log', encoding = 'utf-8', level = logging.WARNING)
    try:
        os.mkdir('.//Output/Plots')
    except FileExistsError:
        pass
    try:
        os.mkdir('.//Output/JoinedCSV')
    except FileExistsError:
        pass    
    try:
        os.mkdir('.//Output/Plots/MS1')
    except FileExistsError:
        pass
    try:
        os.mkdir('.//Output/Plots/MS2')
    except FileExistsError:
        pass
    try:
        os.mkdir('.//Output/Plots/SIM')
    except FileExistsError:
        pass
    # Regex for recognizing the mzML files
    if replicates == False:
        mzmlregex = re.compile('(.*)(.pkl)\Z')
    elif replicates == True:
        mzmlregex = re.compile('(.*)(_\d)(.pkl)\Z')
    if quant == 'y':
        QuantFile = pd.read_excel('Quantconfig.xlsx',header = None, index_col = 0)
    # Iterates through possible scan types which will change data format
    for types in ['MS1', 'MS2', 'SIM']: 
        # Scans for files in the path which are .mzML
        for files in os.listdir('.//{}'.format(str(types))):
            fulldataset = {}
            filecmpdlist =[]
            match = mzmlregex.search(files)
            if match == None:
                continue
            else:
                dataset = pd.read_pickle('.//{}/{}'.format(str(types),str(files)))
                count = count + 1
                # Loads parameters for extracting scan info
                for tup in masslist.itertuples():
                    
                    regressstats ={}
                    scantype = tup.ScanType
                    if not str(scantype) == str(types):
                        continue
                    else:
                        pass
                    Compound = tup.Cmpds
                    StdMass = tup.Std
                    StdMass2 = tup.Std2
                    SILMass = tup.SIL
                    SILMass2 = tup.SIL2
                    SILpre = round(tup.SILPre,4)
                    Stdpre = round(tup.StdPre,4)
                    RTMin = tup.StdRtMin
                    RTMax = tup.StdRtMax
                    SILRTMin = tup.SILRtMin
                    SILRTMax = tup.SILRtMax
                    masswindow = tup.MassWindow
                    if '_' in Compound:
                       cmpdshort = Compound.split('_',1)[0]
                    else:
                        cmpdshort = Compound

                    
    
                    # WidthMax = 3
                    RTAdj = tup.RTAdj
                    # Extractes data from Raw Pkl file
                    if types == 'MS1' or types == 'SIM':
                        C12scaninfoinitialmass1 = dataset.loc[(dataset['m/z'].values >= (float(StdMass) - float(masswindow))) & (dataset['m/z'].values <= (float(StdMass) + float(masswindow)))].sort_values(by = ['scan', 'i'], ascending = [True,False])
                        C12scaninfomass1 = pd.merge(C12scaninfoinitialmass1[['scan','rt','m/z','polarity']], C12scaninfoinitialmass1[['scan','i']].groupby(by = 'scan').sum(), on = 'scan', how = 'right').drop_duplicates(subset = 'scan', keep = 'first')
                        
                        if math.isnan(StdMass2):
                            Stdscaninfo = C12scaninfomass1
                        else:
                            C12scaninfoinitialmass2 = dataset.loc[(dataset['m/z'].values >= (float(StdMass2) - float(masswindow))) & (dataset['m/z'].values <= (float(StdMass2) + float(masswindow)))].sort_values(by = ['scan', 'i'], ascending = [True,False])
                            C12scaninfomass2 = pd.merge(C12scaninfoinitialmass2[['scan','rt','m/z','polarity']], C12scaninfoinitialmass2[['scan','i']].groupby(by = 'scan').sum(), on = 'scan', how = 'right').drop_duplicates(subset = 'scan', keep = 'first')
                            Stdsum = pd.concat([C12scaninfomass1[['scan', 'i']], C12scaninfomass2[['scan', 'i']]]).groupby(by = 'scan').sum()
                            Stdscaninfo = pd.merge(pd.concat([C12scaninfomass1[['scan','rt']], C12scaninfomass2[['scan','rt']]]),Stdsum, on = 'scan', how = 'right').sort_values(by = ['scan', 'i'], ascending = [True,False]).drop_duplicates(subset = 'scan', keep = 'first')
                        
                        
                        C13scaninfoinitialmass1 = dataset.loc[(dataset['m/z'].values >= (float(SILMass) - float(masswindow))) & (dataset['m/z'].values <= (float(SILMass) + float(masswindow)))].sort_values(by = ['scan', 'i'], ascending = [True,False])
                        C13scaninfomass1 = pd.merge(C13scaninfoinitialmass1[['scan','rt','m/z','polarity']], C13scaninfoinitialmass1[['scan','i']].groupby(by = 'scan').sum(), on = 'scan', how = 'right').drop_duplicates(subset = 'scan', keep = 'first')
                        if math.isnan(SILMass2):
                            SILscaninfo = C13scaninfomass1
                        else:
                            C13scaninfoinitialmass2 = dataset.loc[(dataset['m/z'].values >= (float(SILMass2) - float(masswindow))) & (dataset['m/z'].values <= (float(SILMass2) + float(masswindow)))].sort_values(by = ['scan', 'i'], ascending = [True,False])
                            C13scaninfomass2 = pd.merge(C13scaninfoinitialmass2[['scan','rt','m/z','polarity']], C13scaninfoinitialmass2[['scan','i']].groupby(by = 'scan').sum(), on = 'scan', how = 'right').drop_duplicates(subset = 'scan', keep = 'first')
                            Silsum = pd.concat([C13scaninfomass1[['scan', 'i']], C13scaninfomass2[['scan', 'i']]]).groupby(by = 'scan').sum()
                            SILscaninfo = pd.merge(pd.concat([C13scaninfomass1[['scan','rt']], C13scaninfomass2[['scan','rt']]]),Silsum, on = 'scan', how = 'right').sort_values(by = ['scan', 'i'], ascending = [True,False]).drop_duplicates(subset = 'scan', keep = 'first')
                    
                    if types == 'MS2':
                        dataset2 = dataset.loc[dataset['precursor'] == Stdpre]
                        C12scaninfoinitialmass1 = dataset2.loc[(dataset2['m/z'].values >= (float(StdMass) - float(masswindow))) & (dataset2['m/z'].values <= (float(StdMass) + float(masswindow)))].sort_values(by = ['scan', 'i'], ascending = [True,False])
                        C12scaninfomass1 = pd.merge(C12scaninfoinitialmass1[['scan','rt','m/z','polarity']], C12scaninfoinitialmass1[['scan','i']].groupby(by = 'scan').sum(), on = 'scan', how = 'right').drop_duplicates(subset = 'scan', keep = 'first')
                        if math.isnan(StdMass2):
                            Stdscaninfo = C12scaninfomass1
                        else:
                            C12scaninfoinitialmass2 = dataset2.loc[(dataset2['m/z'].values >= (float(StdMass2) - float(masswindow))) & (dataset2['m/z'].values <= (float(StdMass2) + float(masswindow)))].sort_values(by = ['scan', 'i'], ascending = [True,False])
                            C12scaninfomass2 = pd.merge(C12scaninfoinitialmass2[['scan','rt','m/z','polarity']], C12scaninfoinitialmass2[['scan','i']].groupby(by = 'scan').sum(), on = 'scan', how = 'right').drop_duplicates(subset = 'scan', keep = 'first')
                            Stdsum = pd.concat([C12scaninfomass1[['scan', 'i']], C12scaninfomass2[['scan', 'i']]]).groupby(by = 'scan').sum()
                            Stdscaninfo = pd.merge(pd.concat([C12scaninfomass1[['scan','rt']], C12scaninfomass2[['scan','rt']]]),Stdsum, on = 'scan', how = 'right').sort_values(by = ['scan', 'i'], ascending = [True,False]).drop_duplicates(subset = 'scan', keep = 'first')
                        
                        
                        dataset3 = dataset.loc[dataset['precursor'] == SILpre]
                        C13scaninfoinitialmass1 = dataset3.loc[(dataset3['m/z'].values >= (float(SILMass) - float(masswindow))) & (dataset3['m/z'].values <= (float(SILMass) + float(masswindow)))].sort_values(by = ['scan', 'i'], ascending = [True,False])
                        C13scaninfomass1 = pd.merge(C13scaninfoinitialmass1[['scan','rt','m/z','polarity']], C13scaninfoinitialmass1[['scan','i']].groupby(by = 'scan').sum(), on = 'scan', how = 'right').drop_duplicates(subset = 'scan', keep = 'first')
                        if math.isnan(SILMass2):
                            SILscaninfo = C13scaninfomass1
                        else:
                            C13scaninfoinitialmass2 = dataset3.loc[(dataset3['m/z'].values >= (float(SILMass2) - float(masswindow))) & (dataset3['m/z'].values <= (float(SILMass2) + float(masswindow)))].sort_values(by = ['scan', 'i'], ascending = [True,False])
                            C13scaninfomass2 = pd.merge(C13scaninfoinitialmass2[['scan','rt','m/z','polarity']], C13scaninfoinitialmass2[['scan','i']].groupby(by = 'scan').sum(), on = 'scan', how = 'right').drop_duplicates(subset = 'scan', keep = 'first')
                            Silsum = pd.concat([C13scaninfomass1[['scan', 'i']], C13scaninfomass2[['scan', 'i']]]).groupby(by = 'scan').sum()
                            SILscaninfo = pd.merge(pd.concat([C13scaninfomass1[['scan','rt']], C13scaninfomass2[['scan','rt']]]),Silsum, on = 'scan', how = 'right').sort_values(by = ['scan', 'i'], ascending = [True,False]).drop_duplicates(subset = 'scan', keep = 'first')

                    # pd.concat([Stdscaninfo,SILscaninfo]).to_excel('{}_{}.xlsx'.format(str(Compound), str(files).replace('.pkl','')))
                    fullEIC = pd.concat([Stdscaninfo.reset_index(drop = True).add_suffix('_Unlabeled'), SILscaninfo.reset_index(drop = True).add_suffix('_Labeled')], axis=1)
                    fullEIC.to_csv('./Output/JoinedCSV/{}_{}FullEIC.csv'.format(str(match.group(1)) + str(match.group(2)),str(Compound) + '_' + str(scantype)))
                    try:
                        if types == 'MS2' and RTAdj == 'y':
                            Stdscaninfo['ScanActualStd'] = Stdscaninfo['scan']
                            SILscaninfo['ScanActualSIL'] = SILscaninfo['scan']
                            Stdscaninfo['scan'] = range(1, len(Stdscaninfo)+1)
                            SILscaninfo['scan'] = range(1, len(SILscaninfo)+1)
                    except KeyError:
                       Stdscaninfo = pd.DataFrame({'rt':[0], 'scan':[0], 'i':[0]})
                    
                    # Stdscaninfo.to_csv('Std_{}.csv'.format(Compound))
                    # SILscaninfo.to_csv('SIL_{}.csv'.format(Compound))
                    try:
                        Stdscaninfotrim1 = Stdscaninfo[(Stdscaninfo['rt'] >= float(RTMin)) & (Stdscaninfo['rt'] <= float(RTMax))]
                    except KeyError:
                        Stdscaninfotrim1 = pd.DataFrame({'rt':[0], 'scan':[0], 'i':[0]})
                        Stdscaninfo = pd.DataFrame({'rt':[0], 'scan':[0], 'i':[0]})
                    try:
                        SILscaninfotrim1 = SILscaninfo[(SILscaninfo['rt'] >= float(SILRTMin)) & (SILscaninfo['rt'] <= float(SILRTMax))]
                    except KeyError:
                        SILscaninfotrim1 = pd.DataFrame({'rt':[0], 'scan':[0], 'i':[0]})
                        SILscaninfo = pd.DataFrame({'rt':[0], 'scan':[0], 'i':[0]})
                    # Identifies the top intensity peaks within the RT range, if more scans are selected than the RT range, the maximum will be used
                    bestscanwindow, WidthMax = findidealscanwidth(Stdscaninfo, SILscaninfo, RTAdj, types, RTMin, RTMax)
                    Stdscaninfo_sorted = Stdscaninfotrim1.sort_values('i', ascending = False)
                    Stdtop = Stdscaninfo_sorted['scan'].values[:int(WidthMax)]
                    
                        
                    # Averages the scan number and rounds to find the peak maximum
                    try:
                        Stdavg = round(np.average(Stdtop))
                    except ValueError:
                        Stdavg = 0
                    SILscaninfo_sorted = SILscaninfotrim1.sort_values('i', ascending = False)
                    SILtop = SILscaninfo_sorted['scan'].values[:int(WidthMax)]
                    # pd.DataFrame(SILtop).to_csv('Test_{}top.csv'.format(str(Compound)))
                    try:
                        SILavg = round(np.average(SILtop))
                    except ValueError:
                        logging.warning("\n\nNo extracted values for {} in {}, \n check mass accuracy and scan numbers\n".format(str(Compound) + '_' + str(scantype),str(files)))
                        SILavg = 0
                    # Adjusts the SIL peak based on the difference between averages
                    # print(str(SILtop) + '_' + str(Stdtop))
                    if RTAdj == 'y':
                        if SILavg > Stdavg:
                            avgadj = SILavg - Stdavg
                        elif Stdavg > SILavg:
                            avgadj = Stdavg - SILavg
                    elif RTAdj == 'n':
                        avgadj = 0
                        adjmax = SILavg
                    
                    if RTAdj == 'y':
                        if SILavg > Stdavg: 
                            # SILscaninfo ['scan'] = SILscaninfo ['scan'] - float(avgadj)
                            adjmax = SILavg - float(avgadj)
                        elif Stdavg > SILavg:
                            # SILscaninfo ['scan'] = SILscaninfo ['scan'] + float(avgadj)
                            adjmax = SILavg + float(avgadj)
                        else:
                            adjmax = SILavg
                
                        
                    # Trims the full dataset to the number of scans left and right of the peak
                    # Stdpeakleft, Stdpeakright = findrtminmax(Stdscaninfotrim1)
                    # SILpeakleft, SILpeakright = findrtminmax(SILscaninfotrim1)
                    # print(Stdpeakleft, Stdpeakright)
                    # Stdscaninfoint = Stdscaninfo[((Stdscaninfo['rt'] >= float(Stdpeakleft)) & (Stdscaninfo['rt'] <= float(Stdpeakright)))]
                    # SILscaninfoint = SILscaninfo[((SILscaninfo['rt'] >= float(SILpeakleft)) & (SILscaninfo['rt'] <= float(SILpeakright)))]
                    
                    # Stdscaninfoint.to_excel('StdTest_{}.xlsx'.format(str(Compound)))
                    # SILscaninfoint.to_excel('SILTest_{}.xlsx'.format(str(Compound)))
                    Stdscaninfo = Stdscaninfo[(Stdscaninfo['scan'] >= (float(Stdavg) - float(bestscanwindow))) & (Stdscaninfo['scan'] <= (float(Stdavg) + float(bestscanwindow)))]
                    SILscaninfo = SILscaninfo[(SILscaninfo['scan'] >= (float(adjmax) - float(bestscanwindow))) & (SILscaninfo['scan'] <= (float(adjmax) + float(bestscanwindow)))]
                    Stdscaninfo.reset_index(inplace = True, drop = True)
                    SILscaninfo.reset_index(inplace = True, drop = True)
                    
                    # try:
                    #     SilNPTrapz = np.trapz(SILscaninfoint['i'], SILscaninfoint['rt'])
                    #     StdNPTrapz = np.trapz(Stdscaninfoint['i'], Stdscaninfoint['rt'])
                    # except (ValueError,IndexError):
                    #     SilNPTrapz = 0
                    #     StdNPTrapz = 0
                    # try:
                    #     SilSciTrapz = integrate.trapezoid(SILscaninfoint['i'], SILscaninfoint['rt'])
                    #     StdSciTrapz = integrate.trapezoid(Stdscaninfoint['i'], Stdscaninfoint['rt'])
                    # except (ValueError,IndexError):
                    #     SilSciTrapz = 0
                    #     StdSciTrapz = 0
                    # try:
                    #     SilSciSimp = integrate.simpson(SILscaninfoint['i'], SILscaninfoint['rt'])
                    #     StdSciSimp = integrate.simpson(Stdscaninfoint['i'], Stdscaninfoint['rt'])
                    # except (ValueError,IndexError):
                    #     SilSciSimp = 0
                    #     StdSciSimp = 0
                    
                    
                    
                    
                    # Joins the SIL and analyte datasets based on scan number based on previously made adjustment 
                    if len(Stdscaninfo.index) == len (SILscaninfo):
                        joined = pd.merge(Stdscaninfo, SILscaninfo, suffixes=("_Std","_SIL"), on = 'scan', how = 'inner')
                        
                        if optimize == False:
                            joined.to_csv('.//Output/JoinedCSV/{}{}joined.csv'.format(str(match.group(1)) + str(match.group(2)),str(Compound) + '_' + str(scantype)))
                    else:
                        pass
                        if replicates == False:
                            joined = pd.merge(Stdscaninfo, SILscaninfo, suffixes=("_Std","_SIL"), on = 'scan', how = 'inner')
                            
                            if optimize == False:
                                joined.to_csv('.//Output/JoinedCSV/{}_{}joined.csv'.format(str(match.group(1)) + str(match.group(2)),str(Compound) + '_' + str(scantype)))
                        elif replicates == True:
                            joined = pd.merge(Stdscaninfo, SILscaninfo, suffixes=("_Std","_SIL"), on = 'scan', how = 'inner')
                            
                            if optimize == False:
                                joined.to_csv('.//Output/JoinedCSV/{}_{}joined.csv'.format(str(match.group(1)) + str(match.group(2)),str(Compound) + '_' + str(scantype)))
                    if replicates == True:
                        filecmpdlist.append(tuple([str(match.group(1)), str(match.group(2).replace('_','')),str(Compound) + '_' + str(scantype)]))
                    elif replicates == False:
                        filecmpdlist.append(tuple([str(match.group(1)), str(Compound) + '_' + str(scantype)]))
                    # Performs linear regression on the intesity values for the SIL and analyte, error exceptions are necessary for empty datasets from poor extraction
                    x = np.array(joined['i_Std']).reshape(-1,1)
                    y = np.array(joined['i_SIL'])
                    x_int = np.array(joined['i_Std'])
                    x_rt = np.array(joined['rt_SIL']).reshape(-1,1)
                    try:
                        SILHeight = np.amax(SILscaninfo['i'])
                    except ValueError:
                        SILHeight = 0
                    try:
                        StdHeight = np.amax(Stdscaninfo['i'])
                    except ValueError:
                        StdHeight = 0
                    try:
                        HeightRatio = SILHeight/StdHeight
                    except ZeroDivisionError:
                        HeightRatio = 0
                    LR = LinearRegression()
                    try:
                        LR.fit(x,y)
                        if optimize == False:
                            plt.scatter(x,y)
                            plt.plot(x, LR.predict(x))
                            ax = plt.gca()
                            plt.title('{}_{}joined'.format(str(match.group(1)) + str(match.group(2)),str(Compound) + '_' + str(scantype)))
                            plt.xlabel('Analyte Abundance')
                            plt.ylabel('SIL Abundance')
                            plt.text(0.75, 0, 'Slope = {} \n R^2 = {}'.format(str(round(LR.coef_[0], 3)),str(round(LR.score(x,y),3))), transform = ax.transAxes)
                            if types == 'MS1':
                                plt.savefig('.//Output/Plots/MS1/{}{}joined.png'.format(str(match.group(1)) + str(match.group(2)),str(Compound) + '_' + str(scantype)))
                            if types == 'MS2':
                                plt.savefig('.//Output/Plots/MS2/{}{}joined.png'.format(str(match.group(1)) + str(match.group(2)),str(Compound) + '_' + str(scantype)))
                            if types == 'SIM':
                                plt.savefig('.//Output/Plots/SIM/{}{}joined.png'.format(str(match.group(1)) + str(match.group(2)),str(Compound) + '_' + str(scantype)))
                            plt.close()
                            
                            fig, ax1 = plt.subplots()                   
                            axtwin = ax1.twinx()
                           
                            ax1.plot(x_rt,x_int, color = 'b', label = 'unlabeled')
                            axtwin.plot(x_rt,y, color = 'r', label = 'labeled')
                            lines, labels = ax1.get_legend_handles_labels()
                            lines2, labels2 = axtwin.get_legend_handles_labels()
                            axtwin.legend(lines + lines2, labels + labels2, loc=0)
                            fig.suptitle('File: {} Feat: {}'.format(files.strip('.mzML'), str(Compound) + '_' + str(scantype)))
                            ax1.set_xlabel('Retention Time')
                            ax1.set_ylabel('Unlabeled Intensity')
                            axtwin.set_ylabel('Labeled Intensity')
                            
                            fig.savefig('.//Output/Plots/{}/Peak{} {}.png'.format(str(scantype), str(match.group(1)) + str(match.group(2)), str(Compound) + '_' + str(scantype)), bbox_inches='tight')
                            plt.close(fig)
                    except ValueError:
                        slope = 0
                        r2 = 0
                        if optimize == True:
                            pass
                        else:
                            logging.warning("\n\nNo extracted values for {} in {}, \n check mass accuracy and scan numbers\n".format(str(Compound) + '_' + str(scantype),str(files)))
                    else:
                        slope = LR.coef_[0]
                        r2 = LR.score(x,y)
                    if slope < 0:
                        r2 = 0
                        if optimize == True:
                            pass
                        else:
                            logging.warning("\n\nSlope is negative for {} in {}, \n check RT extraction range\n".format(str(Compound) + '_' + str(scantype),str(files)))
                    if len(joined) < 3: 
                        slope = 0
                        r2 = 0
                    try:         
                        if quant == 'y':
                           UnConc = QuantFile.at['Un_{}'.format(str(cmpdshort)),1]
                           StdSlope = QuantFile.at['{}_{}'.format(str(cmpdshort),str(types)),1]
                           SampleMass = QuantFile.at['{}'.format(match.group(1)),1]
                           SpikeVolume = QuantFile.at['Spike Volume (uL)',1]
                           AmountSpiked = StdSlope*UnConc*SpikeVolume
                           totalng = AmountSpiked/slope
                           ngmg = totalng/SampleMass
                        else:
                           ngmg = 0
                    except ZeroDivisionError:
                        ngmg = 0
                    try:
                        if perincorp == 'y':
                            percentincorporated = (SILHeight*100/(SILHeight+StdHeight))
                        else:
                            percentincorporated = 0
                    except ZeroDivisionError:
                        percentincorporated = 0
                    # Creates a results dictionary that is then added to the data export dataframe
                    regressstats = { 'Slope':slope, 'R2':r2, 'SILHeight':SILHeight, 'StdHeight':StdHeight, 'HeightRatio':HeightRatio, 'ng/mg Plant Tissue':ngmg, 'PercentLabelIncorp':percentincorporated, 'Unlabeled Mass':StdMass, 'Unlabeled Mass 2':StdMass2, 'Labeled Mass':SILMass, 'Labeled Mass 2':SILMass2, 'ScanWindowUsed':bestscanwindow, 'WidthUsed':WidthMax}
    
                    fulldataset[str(Compound) + '_' + str(scantype)] = regressstats
                    index = pd.MultiIndex.from_tuples(filecmpdlist)
                    fulldatasetdf = pd.DataFrame(fulldataset, index = ['Slope', 'R2', 'SILHeight', 'StdHeight', 'HeightRatio', 'ng/mg Plant Tissue', 'PercentLabelIncorp', 'Unlabeled Mass', 'Unlabeled Mass 2', 'Labeled Mass', 'Labeled Mass 2','ScanWindowUsed', 'WidthUsed'])
                    fulldatasetdf.columns = index
    # Following steps export data dependant on if there is an optimiztion running or not
            dataexportdf = pd.concat([dataexportdf, fulldatasetdf], axis = 1)
            dataexportdf.fillna(value = 0, inplace = True)
    if replicates == True:
        dataexportdf.to_excel('.//Output/RegressionDataOutput.xlsx')
        DEDFT = dataexportdf.T
        DEDFT.to_csv('.//Output/DataforCalc.csv')
        data = pd.read_csv('.//Output/DataForCalc.csv', names = ['File', 'Replicate', 'Compound', 'Slope', 'R2', 'SILHeight', 'StdHeight', 'HeightRatio', 'ng/mg Plant Tissue','PercentLabelIncorp','Unlabeled Mass', 'Unlabeled Mass 2', 'Labeled Mass', 'Labeled Mass 2','ScanWindowUsed', 'WidthUsed'], header = 0 )
        # For data with replicates, takes the mean and std for each condition set
        meanout = data.groupby(['File','Compound']).mean()
        std = data.groupby(['File','Compound']).std()
        dataexport = pd.concat([meanout,std], keys = ['Mean','Std'])
        dataexport.drop('Replicate', axis = 1, inplace = True)
        dataexport.to_excel('.//Output/ReplicateStats.xlsx')
        os.remove('.//Output/DataforCalc.csv')
    elif replicates == False:
        dataexportdf.T.to_excel('.//Output/RegressionDataOutput.xlsx')
        
        
    logging.shutdown()
    return count
    


def extractfilenames():
    filenames = []
    for types in ['MS1', 'MS2', 'SIM']: 
        # Scans for files in the path which are .mzML
        mzmlregex = re.compile('(.*)(_\d)(.mzML)\Z')
        for files in os.listdir('.//{}'.format(str(types))):
            match = mzmlregex.search(files)
            if match == None:
                continue
            else:
                filenames.append(str(match.group(1)))
    filenamedict = {'Files':filenames}
    Output = pd.DataFrame.from_dict(filenamedict, orient = 'columns')
    Output.drop_duplicates(inplace = True)
    Output.to_excel('FileNames.xlsx')
        
# SILCmpdExtract(masslist = 'Test')
def findidealscanwidth(C12scaninfo, C13scaninfo, RTAdj, types, RTMin, RTMax):
    # Set of scan ranges to test to find ideal scan range
    Widths = [1,2,3,4,5]
    if types == 'MS1':
        scanranges = [10,12, 14, 16, 18, 20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50,52,54,56,58,60]
    else:
        scanranges = [10,12, 14, 16, 18, 20,22,24,26,28,30,32,34,36,38,40,42,44,48,50,52,54,56,58,60,62,64,66,68,70,72,74,76,78,80,90,100,110,120,130,140,150,160,170,180,190,200]
    r2results = []
    scanrangeout = []
    joinlen = []
    widthout = []
    # Iterates through scan ranges to find the best R2
    for width in Widths:
        C12scaninfotrim1 = C12scaninfo[(C12scaninfo['rt'] >= float(RTMin)) & (C12scaninfo['rt'] <= float(RTMax))]
        C13scaninfotrim1 = C13scaninfo[(C13scaninfo['rt'] >= float(RTMin)) & (C13scaninfo['rt'] <= float(RTMax))]
        Stdscaninfo_sorted = C12scaninfotrim1.sort_values('i', ascending = False)
        Stdtop = Stdscaninfo_sorted['scan'].values[:int(width)]
        
            
        # Averages the scan number and rounds to find the peak maximum
        try:
            C12avg = round(np.average(Stdtop))
        except ValueError:
            C12avg = 0
        SILscaninfo_sorted = C13scaninfotrim1.sort_values('i', ascending = False)
        SILtop = SILscaninfo_sorted['scan'].values[:int(width)]
        # pd.DataFrame(SILtop).to_csv('Test_{}top.csv'.format(str(Compound)))
        try:
            C13avg = round(np.average(SILtop))
        except ValueError:
            C13avg = 0
        # Adjusts the SIL peak based on the difference between averages
        # print(str(SILtop) + '_' + str(Stdtop))
        if RTAdj == 'y':
            if C13avg > C12avg:
                avgadj = C13avg - C12avg
            elif C12avg > C13avg:
                avgadj = C12avg - C13avg
        elif RTAdj == 'n':
            avgadj = 0
            adjmax = C13avg
        
        if RTAdj == 'y':
            if C13avg > C12avg: 
                C13scaninfo ['scan'] = C13scaninfo['scan'] - float(avgadj)
                adjmax = C13avg - float(avgadj)
            elif C12avg > C13avg:
                C13scaninfo ['scan'] = C13scaninfo['scan'] + float(avgadj)
                adjmax = C13avg + float(avgadj)
            else:
                adjmax = C13avg 
        for scanwindow in scanranges:
            # Extracts all data based on the previously identified peak height and the current scan window
            C12scaninfo2 = C12scaninfo.loc[(C12scaninfo['scan'].values >= (float(C12avg) - float(scanwindow))) & (C12scaninfo['scan'].values <= (float(C12avg) + float(scanwindow)))]
            if RTAdj == 'y':
                C13scaninfo2 = C13scaninfo.loc[(C13scaninfo['scan'].values >= (float(adjmax) - float(scanwindow))) & (C13scaninfo['scan'].values <= (float(adjmax) + float(scanwindow)))]
            else:
                C13scaninfo2 = C13scaninfo.loc[(C13scaninfo['scan'].values >= (float(C13avg) - float(scanwindow))) & (C13scaninfo['scan'].values <= (float(C13avg) + float(scanwindow)))]    
            # Merges data based on scan number, on interesecting scans, then arranges the data for regression
            joined = pd.merge(C12scaninfo2, C13scaninfo2, suffixes=("_Unlabeled","_Labeled"), on = 'scan', how = 'inner')
            joinedlength = len(joined)
            x = np.array(joined['i_Unlabeled']).reshape(-1,1)
            y = np.array(joined['i_Labeled'])
            # Performs regression on the current scan range data and finds the R2 and appends to output, except for an empty dataframe which gives an r2 of 0
            try:
                LRtest = LinearRegression()
                LRtest.fit(x,y)
                r2test = LRtest.score(x,y)
                
                
               
                r2testout = r2test
            except ValueError:
                r2testout = 0
                
            r2results.append(r2testout)
            scanrangeout.append(scanwindow)
            joinlen.append(joinedlength)
            widthout.append(width)
    # Sorts the R2 results and returns the best scan range
    results = pd.DataFrame({'Scans':scanrangeout, 'R2':r2results, 'len':joinlen, 'Width':widthout})
    results['R2Round'] = results['R2'].round(1)
    # results.to_excel('.//Test Lengths/Test_{}.xlsx'.format(featpair))
    results.sort_values(by = ['R2Round','Scans'], ascending = [False, False], ignore_index = True, inplace = True)
    results = results.drop('R2Round', axis = 1)
    # results.sort_values(by = ['R2'], ascending = False, ignore_index = True, inplace = True)
    bestwindow = results['Scans'][0]
    bestwidth = results['Width'][0]
    
    return bestwindow, bestwidth


def findrtminmax(data):
    if len(data) == 0:
        leftrt = 0
        rightrt = 0
        return leftrt, rightrt
    data = data.reset_index()
    data = data.drop([0])
    data = data.reset_index()
    if len(data) == 0:
        leftrt = 0
        rightrt = 0
        return leftrt, rightrt
    
    peak = signal.find_peaks(data['i'], prominence = [0,None])
    try:
        highestpeakindex = np.argmax(peak[1]['prominences'])
    except ValueError:
        leftrt = 0
        rightrt = 0
        return leftrt, rightrt
    left = peak[1]['left_bases'][highestpeakindex]
    right = peak[1]['right_bases'][highestpeakindex]
    leftrt = data['rt'][left]
    rightrt = data['rt'][right]
    return leftrt,rightrt



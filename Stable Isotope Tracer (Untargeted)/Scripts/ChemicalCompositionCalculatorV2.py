# -*- coding: utf-8 -*-
"""
Created on Wed Jul 12 15:51:16 2023

@author: Evan Larson, Hegeman & Cohen Labs - University of Minnesota
"""
import pandas as pd
import re

import math
from p_tqdm import p_umap
import numpy as np


def roundmass(mass, massbin):
    return(round(round(mass / massbin) * massbin, 3))



'Database Searching Below Here'


def searchdatabase(wd, mass, c13mass, polarity, carbons, ppm, full, massbin, preferCHNO, label):
    # Runs Formula Matching Function
    match = matchformuladatabase(mass, c13mass, polarity, carbons, ppm, label, full)
    # Rounds the masses into given bin
    c12round = roundmass(mass, massbin)
    c13round = roundmass(c13mass, massbin)
    # Outputs data to excel if there is more than one match
    if len(match) >= 1:
        match.to_excel('{}/Output/Peaks/Formulas/{}_{}.xlsx'.format(wd,str(c12round), str(c13round)))
    
    if len(match) > 1:
        # Only Prefers to output formulas containing CHNO for main program output
        if preferCHNO == 'Yes':

            matchfilter = match[~match['Formula'].str.contains('P')]
            matchfilter = matchfilter[~matchfilter['Formula'].str.contains(
                'S')]
            matchfilter = matchfilter[~matchfilter['Formula'].str.contains(
                'Cl')]
            matchfilter = matchfilter[~matchfilter['Formula'].str.contains(
                'Br')]
            matchfilter = matchfilter[~matchfilter['Formula'].str.contains(
                'B')]
            matchfilter = matchfilter[~matchfilter['Formula'].str.contains(
                'F')]
            matchfilter.reset_index(inplace=True, drop=True)
            # Sorts the output formulas by the absppm while also dropping that column after
            if len(matchfilter) > 1:
                sortmatch = matchfilter.sort_values(by='absppm')
                try:
                    matchtrunc = sortmatch.truncate(after=0)
                except ValueError:
                    print( sortmatch, matchfilter, match)
                dropmatch = matchtrunc.drop('absppm', axis=1)
            if len(matchfilter) == 1:
                dropmatch = matchfilter.drop('absppm', axis=1)
            if len(matchfilter) == 0:
                matchdict = {'Formula': ['N/A'], 'ppm': ['N/A'], 'Adduct': ['N/A'], 'Labeledppm': ['N/A'], 'Nat13C1': [0], 'TheoMassUnlabeled': ['0']}
                dropmatch = pd.DataFrame.from_dict(matchdict)

            return dropmatch
        else:
            # Sorts the output formulas by the absppm while also dropping that column after
            if len(match) > 1:
               
                sortmatch = match.sort_values(by='absppm')
                matchtrunc = sortmatch.truncate(after=0)
                dropmatch = matchtrunc.drop('absppm', axis=1)
                return dropmatch
    elif len(match) == 0:
        matchdict = {'Formula': ['N/A'], 'ppm': ['N/A'], 'Adduct': ['N/A'], 'Labeledppm': ['N/A'], 'Nat13C1': [0], 'TheoMassUnlabeled': [0]}
        nomatch = pd.DataFrame.from_dict(matchdict)
        return nomatch
    elif len(match) == 1:
        dropmatch = match.drop('absppm', axis=1)
        return dropmatch


def matchformuladatabase(mass, c13mass, polarity, carbons, ppm, label, full):
    # Dataframe to output formula matches
    sortFrame = pd.DataFrame()
    # Reads in formula database
    formdb = pd.read_pickle('.//Database/FullFormulaDataBase.pkl')    
    # calculates the ppm / 10e6 so the calculation only needs to be performed once
    ppmmil = ppm / 10**6
    # Deoending on the type of experiment (Fully labeled vs. Partially Labeled), the polarity, and the label type slices the database dataframe to identify matches
    if full == 'Fully':
        if polarity == 'POSITIVE':
            if label == '13C':
                Hmatches = formdb.loc[(formdb['C'].values == carbons) & ((((formdb['M+H Exact'].values) - (
                    ppmmil*formdb['M+H Exact'].values)) <= mass) & ((formdb['M+H Exact'].values) + (ppmmil*formdb['M+H Exact'].values) >= mass))]
                Namatches = formdb.loc[(formdb['C'].values == carbons) & ((((formdb['M+Na Exact'].values) - (
                    ppmmil*formdb['M+Na Exact'].values)) <= mass) & ((formdb['M+Na Exact'].values) + (ppmmil*formdb['M+Na Exact'].values) >= mass))]
            if label == '15N':
                Hmatches = formdb.loc[(formdb['N'].values == carbons) & ((((formdb['M+H Exact'].values) - (
                    ppmmil*formdb['M+H Exact'].values)) <= mass) & ((formdb['M+H Exact'].values) + (ppmmil*formdb['M+H Exact'].values) >= mass))]
                Namatches = formdb.loc[(formdb['N'].values == carbons) & ((((formdb['M+Na Exact'].values) - (
                    ppmmil*formdb['M+Na Exact'].values)) <= mass) & ((formdb['M+Na Exact'].values) + (ppmmil*formdb['M+Na Exact'].values) >= mass))]
            if label == '2H':
                Hmatches = formdb.loc[(formdb['H'].values == carbons) & ((((formdb['M+H Exact'].values) - (
                    ppmmil*formdb['M+H Exact'].values)) <= mass) & ((formdb['M+H Exact'].values) + (ppmmil*formdb['M+H Exact'].values) >= mass))]
                Namatches = formdb.loc[(formdb['H'].values == carbons) & ((((formdb['M+Na Exact'].values) - (
                    ppmmil*formdb['M+Na Exact'].values)) <= mass) & ((formdb['M+Na Exact'].values) + (ppmmil*formdb['M+Na Exact'].values) >= mass))]
            if label == '18O':
                Hmatches = formdb.loc[(formdb['O'].values == carbons) & ((((formdb['M+H Exact'].values) - (
                    ppmmil*formdb['M+H Exact'].values)) <= mass) & ((formdb['M+H Exact'].values) + (ppmmil*formdb['M+H Exact'].values) >= mass))]
                Namatches = formdb.loc[(formdb['O'].values == carbons) & ((((formdb['M+Na Exact'].values) - (
                    ppmmil*formdb['M+Na Exact'].values)) <= mass) & ((formdb['M+Na Exact'].values) + (ppmmil*formdb['M+Na Exact'].values) >= mass))]
            if label == '34S':
                Hmatches = formdb.loc[(formdb['S'].values == carbons) & ((((formdb['M+H Exact'].values) - (
                    ppmmil*formdb['M+H Exact'].values)) <= mass) & ((formdb['M+H Exact'].values) + (ppmmil*formdb['M+H Exact'].values) >= mass))]
                Namatches = formdb.loc[(formdb['S'].values == carbons) & ((((formdb['M+Na Exact'].values) - (
                    ppmmil*formdb['M+Na Exact'].values)) <= mass) & ((formdb['M+Na Exact'].values) + (ppmmil*formdb['M+Na Exact'].values) >= mass))]
        if polarity == 'NEGATIVE':
            if label == '13C':
                Negmatches = formdb.loc[(formdb['C'].values == carbons) & ((((formdb['M-H Exact'].values) - (
                    ppmmil*formdb['M-H Exact'].values)) <= mass) & ((formdb['M-H Exact'].values) + (ppmmil*formdb['M-H Exact'].values) >= mass))]
            if label == '15N':
                Negmatches = formdb.loc[(formdb['N'].values == carbons) & ((((formdb['M-H Exact'].values) - (
                    ppmmil*formdb['M-H Exact'].values)) <= mass) & ((formdb['M-H Exact'].values) + (ppmmil*formdb['M-H Exact'].values) >= mass))]
            if label == '2H':
                Negmatches = formdb.loc[(formdb['H'].values == carbons) & ((((formdb['M-H Exact'].values) - (
                    ppmmil*formdb['M-H Exact'].values)) <= mass) & ((formdb['M-H Exact'].values) + (ppmmil*formdb['M-H Exact'].values) >= mass))]
            if label == '18O':
                Negmatches = formdb.loc[(formdb['O'].values == carbons) & ((((formdb['M-H Exact'].values) - (
                    ppmmil*formdb['M-H Exact'].values)) <= mass) & ((formdb['M-H Exact'].values) + (ppmmil*formdb['M-H Exact'].values) >= mass))]
            if label == '34S':
                Negmatches = formdb.loc[(formdb['S'].values == carbons) & ((((formdb['M-H Exact'].values) - (
                    ppmmil*formdb['M-H Exact'].values)) <= mass) & ((formdb['M-H Exact'].values) + (ppmmil*formdb['M-H Exact'].values) >= mass))]
    if full == 'Partially':
        if polarity == 'POSITIVE':
            if label == '13C':
                Hmatches = formdb.loc[(formdb['C'].values >= carbons) & ((((formdb['M+H Exact'].values) - (
                    ppmmil*formdb['M+H Exact'].values)) <= mass) & ((formdb['M+H Exact'].values) + (ppmmil*formdb['M+H Exact'].values) >= mass))]
                Namatches = formdb.loc[(formdb['C'].values >= carbons) & ((((formdb['M+Na Exact'].values) - (
                    ppmmil*formdb['M+Na Exact'].values)) <= mass) & ((formdb['M+Na Exact'].values) + (ppmmil*formdb['M+Na Exact'].values) >= mass))]
            if label == '15N':
                Hmatches = formdb.loc[(formdb['N'].values >= carbons) & ((((formdb['M+H Exact'].values) - (
                    ppmmil*formdb['M+H Exact'].values)) <= mass) & ((formdb['M+H Exact'].values) + (ppmmil*formdb['M+H Exact'].values) >= mass))]
                Namatches = formdb.loc[(formdb['N'].values >= carbons) & ((((formdb['M+Na Exact'].values) - (
                    ppmmil*formdb['M+Na Exact'].values)) <= mass) & ((formdb['M+Na Exact'].values) + (ppmmil*formdb['M+Na Exact'].values) >= mass))]
            if label == '2H':
                Hmatches = formdb.loc[(formdb['H'].values >= carbons) & ((((formdb['M+H Exact'].values) - (
                    ppmmil*formdb['M+H Exact'].values)) <= mass) & ((formdb['M+H Exact'].values) + (ppmmil*formdb['M+H Exact'].values) >= mass))]
                Namatches = formdb.loc[(formdb['H'].values >= carbons) & ((((formdb['M+Na Exact'].values) - (
                    ppmmil*formdb['M+Na Exact'].values)) <= mass) & ((formdb['M+Na Exact'].values) + (ppmmil*formdb['M+Na Exact'].values) >= mass))]
            if label == '18O':
                Hmatches = formdb.loc[(formdb['O'].values >= carbons) & ((((formdb['M+H Exact'].values) - (
                    ppmmil*formdb['M+H Exact'].values)) <= mass) & ((formdb['M+H Exact'].values) + (ppmmil*formdb['M+H Exact'].values) >= mass))]
                Namatches = formdb.loc[(formdb['O'].values >= carbons) & ((((formdb['M+Na Exact'].values) - (
                    ppmmil*formdb['M+Na Exact'].values)) <= mass) & ((formdb['M+Na Exact'].values) + (ppmmil*formdb['M+Na Exact'].values) >= mass))]
            if label == '34S':
                Hmatches = formdb.loc[(formdb['S'].values >= carbons) & ((((formdb['M+H Exact'].values) - (
                    ppmmil*formdb['M+H Exact'].values)) <= mass) & ((formdb['M+H Exact'].values) + (ppmmil*formdb['M+H Exact'].values) >= mass))]
                Namatches = formdb.loc[(formdb['S'].values >= carbons) & ((((formdb['M+Na Exact'].values) - (
                    ppmmil*formdb['M+Na Exact'].values)) <= mass) & ((formdb['M+Na Exact'].values) + (ppmmil*formdb['M+Na Exact'].values) >= mass))]
        if polarity == 'NEGATIVE':
            if label == '13C':
                Negmatches = formdb.loc[(formdb['C'].values >= carbons) & ((((formdb['M-H Exact'].values) - (
                    ppmmil*formdb['M-H Exact'].values)) <= mass) & ((formdb['M-H Exact'].values) + (ppmmil*formdb['M-H Exact'].values) >= mass))]
            if label == '15N':
                Negmatches = formdb.loc[(formdb['N'].values >= carbons) & ((((formdb['M-H Exact'].values) - (
                    ppmmil*formdb['M-H Exact'].values)) <= mass) & ((formdb['M-H Exact'].values) + (ppmmil*formdb['M-H Exact'].values) >= mass))]
            if label == '2H':
                Negmatches = formdb.loc[(formdb['H'].values >= carbons) & ((((formdb['M-H Exact'].values) - (
                    ppmmil*formdb['M-H Exact'].values)) <= mass) & ((formdb['M-H Exact'].values) + (ppmmil*formdb['M-H Exact'].values) >= mass))]
            if label == '18O':
                Negmatches = formdb.loc[(formdb['O'].values >= carbons) & ((((formdb['M-H Exact'].values) - (
                    ppmmil*formdb['M-H Exact'].values)) <= mass) & ((formdb['M-H Exact'].values) + (ppmmil*formdb['M-H Exact'].values) >= mass))]
            if label == '34S':
                Negmatches = formdb.loc[(formdb['S'].values >= carbons) & ((((formdb['M-H Exact'].values) - (
                    ppmmil*formdb['M-H Exact'].values)) <= mass) & ((formdb['M-H Exact'].values) + (ppmmil*formdb['M-H Exact'].values) >= mass))]
    # collects data and calculates ppm error for each match, also calculates theoretical 13C intensity which adds all to the output dataframe at the end, also calculates the theoretical labeled mass based on number of labels and uses that to calculate labeled ppm error 
    if polarity == 'POSITIVE':
        Hmatches.reset_index(inplace = True)
        Namatches.reset_index(inplace = True)
        for matches in range(len(Hmatches)):
            addmatch = pd.DataFrame()
            matched = {}
            formula = Hmatches.loc[matches, 'Formula']
            adduct = ['[M+H]']
            theomass = Hmatches.loc[matches, 'M+H Exact']
            ppmerror = (((theomass - mass)/theomass)*(10**6))
            absppm = abs(ppmerror)
            
            if carbons == 1 and label == '13C':
                Cmatch = re.compile('(C\d*)')
                if not re.search(Cmatch, formula) == None:
                  
                    if re.search(Cmatch, formula).group(1) == 'C':
                        Cs = 1
                    elif re.search(Cmatch, formula).group(1) == '':
                        Cs = 0
                    
                    else:
                        Cs = int(re.search(Cmatch, formula).group(1).replace('C',''))
                else:
                    Cs = 0
                Nat13C1 = (Cs * 1.08)
            else:
                Nat13C1 = 0
            if label == '13C':
                theolabel = float(theomass) + (1.003355 * carbons)
            if label == '15N':
                theolabel = float(theomass) + (0.997035 * carbons)
            if label == '2H':
                theolabel = float(theomass) + (1.006276 * carbons)
            if label == '18O':
                theolabel = float(theomass) + (2.004245 * carbons)
            if label == '34S':
                theolabel = float(theomass) + (1.99579 * carbons)
            FinalMassHigh = theolabel + ((ppm / 10**6)*theolabel)
            FinalMassLow = theolabel - ((ppm / 10**6)*theolabel)
            if c13mass >= FinalMassLow and c13mass <= FinalMassHigh:
                quickppm = (((theolabel - c13mass)/theolabel)*(10**6))
            else:
                quickppm = 'OutsideLimit'
            matched = {'Formula': formula, 'Adduct': adduct, 'ppm': ppmerror,
                       'absppm': absppm, 'Nat13C1': Nat13C1, 'Labeledppm':quickppm, 'TheoMassUnlabeled': theomass}
            addmatch = pd.DataFrame(matched, index=[0])
            sortFrame = pd.concat([sortFrame, addmatch],
                                  axis=0, ignore_index=True)
        for matches in range(len(Namatches)):
            addmatch = pd.DataFrame()
            matched = {}
            formula = Namatches.loc[matches, 'Formula']
            adduct = ['[M+Na]']
            theomass = Namatches.loc[matches, 'M+Na Exact']
            ppmerror = (((theomass - mass)/theomass)*(10**6))
            absppm = abs(ppmerror)
            if carbons == 1 and label == '13C':
                Cmatch = re.compile('(C\d*)')
                if not re.search(Cmatch, formula) == None:
                  
                    if re.search(Cmatch, formula).group(1) == 'C':
                        Cs = 1
                    elif re.search(Cmatch, formula).group(1) == '':
                        Cs = 0
                    
                    else:
                        Cs = int(re.search(Cmatch, formula).group(1).replace('C',''))
                else:
                    Cs = 0
                Nat13C1 = (Cs * 1.08)
            else:
                Nat13C1 = 0
            if label == '13C':
                theolabel = float(theomass) + (1.003355 * carbons)
            if label == '15N':
                theolabel = float(theomass) + (0.997035 * carbons)
            if label == '2H':
                theolabel = float(theomass) + (1.006276 * carbons)
            if label == '18O':
                theolabel = float(theomass) + (2.004245 * carbons)
            if label == '34S':
                theolabel = float(theomass) + (1.99579 * carbons)
            FinalMassHigh = theolabel + ((ppm / 10**6)*theolabel)
            FinalMassLow = theolabel - ((ppm / 10**6)*theolabel)
            if c13mass >= FinalMassLow and c13mass <= FinalMassHigh:
                quickppm = (((theolabel - c13mass)/theolabel)*(10**6))
            else:
                quickppm = 'OutsideLimit'
            matched = {'Formula': formula, 'Adduct': adduct, 'ppm': ppmerror,
                       'absppm': absppm, 'Nat13C1': Nat13C1, 'Labeledppm':quickppm, 'TheoMassUnlabeled': theomass}
            addmatch = pd.DataFrame(matched, index=[0])
            sortFrame = pd.concat([sortFrame, addmatch],
                                  axis=0, ignore_index=True)
    
    if polarity == 'NEGATIVE':
        Negmatches.reset_index(inplace = True)
        for matches in range(len(Negmatches)):
            addmatch = pd.DataFrame()
            matched = {}
            formula = Negmatches.loc[matches, 'Formula']
            adduct = ['[M-H]']
            theomass = Negmatches.loc[matches, 'M-H Exact']
            ppmerror = (((theomass - mass)/theomass)*(10**6))
            absppm = abs(ppmerror)
            if carbons == 1 and label == '13C':
                Cmatch = re.compile('(C\d*)')
                if not re.search(Cmatch, formula) == None:
                  
                    if re.search(Cmatch, formula).group(1) == 'C':
                        Cs = 1
                    elif re.search(Cmatch, formula).group(1) == '':
                        Cs = 0
                    
                    else:
                        Cs = int(re.search(Cmatch, formula).group(1).replace('C',''))
                else:
                    Cs = 0
                Nat13C1 = (Cs * 1.08)
            else:
                Nat13C1 = 0
            if label == '13C':
                theolabel = float(theomass) + (1.003355 * carbons)
            if label == '15N':
                theolabel = float(theomass) + (0.997035 * carbons)
            if label == '2H':
                theolabel = float(theomass) + (1.006276 * carbons)
            if label == '18O':
                theolabel = float(theomass) + (2.004245 * carbons)
            if label == '34S':
                theolabel = float(theomass) + (1.99579 * carbons)
            FinalMassHigh = theolabel + ((ppm / 10**6)*theolabel)
            FinalMassLow = theolabel - ((ppm / 10**6)*theolabel)
            if c13mass >= FinalMassLow and c13mass <= FinalMassHigh:
                quickppm = (((theolabel - c13mass)/theolabel)*(10**6))
            else:
                quickppm = 'OutsideLimit'
            matched = {'Formula': formula, 'Adduct': adduct, 'ppm': ppmerror,
                       'absppm': absppm, 'Nat13C1': Nat13C1, 'Labeledppm':quickppm, 'TheoMassUnlabeled': theomass}
            addmatch = pd.DataFrame(matched, index=[0])
            sortFrame = pd.concat([sortFrame, addmatch],
                                  axis=0, ignore_index=True)
    try:
        sortFrame.sort_values(by='absppm', inplace=True, ignore_index=True)
    except:
        pass

    return sortFrame



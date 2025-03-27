# -*- coding: utf-8 -*-
"""
Created on Fri Jan 31 11:47:00 2025

@author: Evan Larson, Hegeman & Cohen Labs - University of Minnesota
"""
import re,os,pandas as pd
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
    
extractfilenames()
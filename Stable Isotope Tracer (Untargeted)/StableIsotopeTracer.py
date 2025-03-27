# -*- coding: utf-8 -*-
"""
Created on Mon Feb 26 14:40:26 2024

@author: Evan Larson, Hegeman & Cohen Labs - University of Minnesota
"""

import sys, re, pandas as pd
from PyQt5.QtWidgets import QMainWindow, QApplication, QDialog, QFileDialog
from Scripts.RunGuiv2 import Ui_MainWindow
from PyQt5.QtGui import QValidator
from Scripts.OverallRun import GUIRUN
from Scripts.SequenceRun import Ui_Dialog
from Scripts.SequencedRunFromGui import SequenceRun
from PyQt5 import QtCore
from Scripts.RunUpdate import Ui_Dialog as Ui_Dialog2
import multiprocessing as mp


class Window(QMainWindow, Ui_MainWindow):
    def __init__(self, parent=None):
        super().__init__()
        self.setupUi(self)
        self.RUN.clicked.connect(self.validateentries)
        self.entries = [self.FDMassBinSize.text(),self.FDMassWindowLabel.text(),self.FDMaxNumLabels.text(),self.FDMinIntensity.text(),self.FDMinNumLabels.text(),self.FDMinScan.text(),self.FDRTBinSize.text(), self.DSPPMError.text(),self.DSTargetFormulaFile.text(),self.RegressMassWindow.text(),self.RegressRTWindow.text(),self.RegressionMinR2.text()]
        self.numericalentries = [self.FDMassBinSize.text(), self.FDMassWindowLabel.text(), self.FDMaxNumLabels.text(), self.FDMinIntensity.text(), self.FDMinNumLabels.text(), self.DSPPMError.text(), self.RegressMassWindow.text(), self.RegressRTWindow.text(), self.RegressionMinR2.text()]
        self.FDMassBinSize.setValidator(NumericalDecimal()) 
        self.FDMassWindowLabel.setValidator(NumericalDecimal()) 
        self.FDMaxNumLabels.setValidator(Numerical())
        self.FDMinIntensity.setValidator(NumericalDecimal()) 
        self.FDMinNumLabels.setValidator(Numerical()) 
        self.DSPPMError.setValidator(NumericalDecimal()) 
        self.RegressMassWindow.setValidator(NumericalDecimal()) 
        self.RegressRTWindow.setValidator(NumericalDecimal()) 
        self.RegressionMinR2.setValidator(NumericalDecimal())
        self.FDMinScan.setValidator(Numerical())
        self.FDRTBinSize.setValidator(NumericalDecimal())
        self.SaveConfig.clicked.connect(self.SaveFile)
        self.SaveConfig_2.clicked.connect(self.LoadFile)
        self.Message.setAlignment(QtCore.Qt.AlignCenter)
    def RUN_clicked(self):
        
        FDMassBinSize = float(self.FDMassBinSize.text())
        FDMassWindowLabel = float(self.FDMassWindowLabel.text())
        FDMaxNumLabels = int(self.FDMaxNumLabels.text())
        FDMinIntensity = float(self.FDMinIntensity.text())
        FDMinNumLabels = int(self.FDMinNumLabels.text())
        FDMinScan = int(self.FDMinScan.text())
        FDRTBinSize = float(self.FDRTBinSize.text())
        Function = str(self.Function.currentText())
        try:
            float(self.Label.currentText())
            Label = float(self.Label.currentText())
        except ValueError:
            Label = str(self.Label.currentText())
        # MzmineLabeled = str(self.MzmineLabeled.text())
        # Mzmineunlabeled = str(self.Mzmineunlabeled.text())
        # SingleFileOrSequence = self.SingleFileOrSequence.currentText()
        DSFullPartial = str(self.DSFullPartial.currentText())
        DSPPMError = float(self.DSPPMError.text())
        DSTargetFormulaFile = str(self.DSTargetFormulaFile.text())
        DSTargetList = str(self.DSTargetList.currentText())
        RegressKeepPlots = str(self.RegressKeepPlots.currentText())
        RegressMassWindow = float(self.RegressMassWindow.text())
        RegressRTWindow = float(self.RegressRTWindow.text())
        RegressionMinR2 = float(self.RegressionMinR2.text())
        PremadeList = str(self.PremadeList.text()).replace('.csv','')
        PremadeCSV = str(self.PremadeList.text())
        PreferCHNO = str(self.DSTargetList_2.currentText())
        UnlabeledControl = str(self.UnlabaledControl.currentText())
        # datafolder = QFileDialog.getExistingDirectory(self, 'Choose Data Folder...', '.\\')
        # print(datafolder)
        GUIRUN(featuretype = Function, premadelist = PremadeList, premadexlsx = PremadeCSV, mincarbons = FDMinNumLabels, maxcarbons = FDMaxNumLabels, masswindow = RegressMassWindow,rtbin=FDRTBinSize, massbin=FDMassBinSize, rtwindow=RegressRTWindow, removeplots=RegressKeepPlots, C12C13Window=FDMassWindowLabel, formppm=DSPPMError, fully = DSFullPartial, minR2=RegressionMinR2,minint=FDMinIntensity,label=Label,targetforms=DSTargetList,target=DSTargetFormulaFile,minscan=FDMinScan,preferCHNO=PreferCHNO, UnlabeledControl = UnlabeledControl)
    
    def validateentries(self):
        self.updateentries()
        self.Message.setText('Starting Analysis')
        if self.SingleFileOrSequence.currentText() == 'Single Set':    
            if self.Function.currentText() == 'Work From Mzmine Feature List':
                if any(items == '' for items in self.entries) and self.DSTargetList.currentText() == 'Yes':
                   self.Message.setText('Missing Required Field')
                   return
                elif any('.csv' not in items for items in self.csventries):
                        self.Message.setText('Mzmine Files must be CSV')
                        return
                elif self.DSTargetList.currentText() == 'Yes' and '.xlsx' not in self.DSTargetFormulaFile.text():
                    self.Message.setText('Target File must be XLSX')
                    return
                elif checklabelinput(self.Label.currentText()) == False:
                    self.Message.setText('Incorrect Label Entry')
                    return
                else: 
                    self.Message.setText('Check Python Terminal')
                    self.runwindow()
                    # self.RUN_clicked()
                    # self.Message.setText('Analysis Finished! Check Output Folder')
                    return
            elif self.Function.currentText() == 'Work From Premade Feature List':
                if any(items == '' for items in self.entries) and self.DSTargetList.currentText() == 'Yes' and self.PremadeList.text() == '':
                   self.Message.setText('Missing Required Field')
                   return
                elif '.csv' not in self.PremadeList.text():
                        self.Message.setText('Feature Lists must be CSV')
                        return
                elif self.DSTargetList.currentText() == 'Yes' and '.xlsx' not in self.DSTargetFormulaFile.text():
                    self.Message.setText('Target File must be XLSX')
                    return
                elif checklabelinput(self.Label.currentText())== False:
                    self.Message.setText('Incorrect Label Entry')
                    return
                else:    
                    self.Message.setText('Check Python Terminal')
                    self.runwindow()
                    # self.RUN_clicked()
                    # self.Message.setText('Analysis Finished! Check Output Folder')
                    return
            elif any(items == '' for items in self.nonmzmineentries) and self.DSTargetList.currentText() == 'Yes':
                self.Message.setText('Missing Required Field')
                return
            elif self.DSTargetList.currentText() == 'Yes' and '.xlsx' not in self.DSTargetFormulaFile.text():
                self.Message.setText('Target File must be XLSX')
                return
            elif checklabelinput(self.Label.currentText()) == False:
                self.Message.setText('Incorrect Label Entry')
                return
            else:    
                self.Message.setText('Check Python Terminal')
                self.runwindow()
                # self.RUN_clicked()
                # self.Message.setText('Analysis Finished! Check Output Folder')
                return
                
        else:
            self.Message.setText('Check Python Terminal')
            self.opensequencewindow()
            return
    
    def runwindow(self):
        runwindowopen = RunWindow()
        runwindowopen.exec()
        return
    def opensequencewindow(self):
        seqwindow = SequenceWindow()
        seqwindow.exec()
    def updateentries(self):
        self.entries = [self.FDMassBinSize.text(),self.FDMassWindowLabel.text(),self.FDMaxNumLabels.text(),self.FDMinIntensity.text(),self.FDMinNumLabels.text(),self.FDMinScan.text(),self.FDRTBinSize.text(),self.DSPPMError.text(),self.DSTargetFormulaFile.text(),self.RegressMassWindow.text(),self.RegressRTWindow.text(),self.RegressionMinR2.text()]
        self.nonmzmineentries = [self.FDMassBinSize.text(),self.FDMassWindowLabel.text(),self.FDMaxNumLabels.text(),self.FDMinIntensity.text(),self.FDMinNumLabels.text(),self.FDMinScan.text(),self.FDRTBinSize.text(),self.DSPPMError.text(),self.DSTargetFormulaFile.text(),self.RegressMassWindow.text(),self.RegressRTWindow.text(),self.RegressionMinR2.text()]
    def SaveFile(self):
        self.Message.setText('Saving...')
        filename = QFileDialog.getSaveFileName(self, 'Save Config File', '.\\',"CSV File (*.csv)")
        if filename[0] == '':
            return
        else:
            SaveDict = {'SingeDataSet or Sequence': [self.SingleFileOrSequence.currentText()],'Unlabeled Control':[self.UnlabaledControl.currentText()], 'Program Function':[self.Function.currentText()], 'Select Isotopic Label':[self.Label.currentText()],'PremadeList':[self.PremadeList.text()],'PreferCHNO':[self.DSTargetList_2.currentText()], 'Minimum Number of Labels':[self.FDMinNumLabels.text()],'Maximum Number of Labels':[self.FDMaxNumLabels.text()],'Mass Window for Label':[self.FDMassWindowLabel.text()],'Retention Time Bin Size':[self.FDRTBinSize.text()],'Mass Bin Size':[self.FDMassBinSize.text()],'Minimum Intensity for Features':[self.FDMinIntensity.text()],'Minimum Binned Scans for Features':[self.FDMinScan.text()],'Retention Time Window For Peak Height':[self.RegressRTWindow.text()],'Regression Mass Window':[self.RegressMassWindow.text()],'Minimum R2':[self.RegressionMinR2.text()],'RemovePlots':[self.RegressKeepPlots.currentText()], 'PPM Error':[self.DSPPMError.text()],'FullyorPartially':[self.DSFullPartial.currentText()],'TargetFormulas':[self.DSTargetFormulaFile.text()],'FormulaListFile':[self.DSTargetList.currentText()]  }
            pd.DataFrame.from_dict(SaveDict, orient = 'index', columns = ['Parameters']).to_csv(filename[0])
        self.Message.setText('Saved')
        
    def LoadFile(self):
        self.Message.setText('Loading...')
        filename = QFileDialog.getOpenFileName(self, 'Load Config File', '.\\', "CSV File (*.csv)")
        if filename[0] == '':
            return
        else:
            file = pd.read_csv(filename[0], index_col = 0)
    
            sfosindex = self.SingleFileOrSequence.findText(file.loc['SingeDataSet or Sequence','Parameters'])
            self.SingleFileOrSequence.setCurrentIndex(sfosindex)
            funindex = self.Function.findText(file.loc['Program Function','Parameters'])
            self.Function.setCurrentIndex(funindex)
            labindex = self.Label.findText(file.loc['Select Isotopic Label','Parameters'])
            self.Label.setCurrentIndex(labindex)
            unlabeledindex = self.UnlabaledControl.findText(file.loc['Unlabeled Control', 'Parameters'])
            self.UnlabaledControl.setCurrentIndex(unlabeledindex)
            self.FDMinNumLabels.setText(str(file.loc['Minimum Number of Labels','Parameters']))
            self.FDMaxNumLabels.setText(str(file.loc['Maximum Number of Labels','Parameters']))
            self.FDMassWindowLabel.setText(str(file.loc['Mass Window for Label','Parameters']))
            self.FDRTBinSize.setText(str(file.loc['Retention Time Bin Size','Parameters']))
            self.FDMassBinSize.setText(str(file.loc['Mass Bin Size','Parameters']))
            self.FDMinIntensity.setText(str(file.loc['Minimum Intensity for Features','Parameters']))
            self.FDMinScan.setText(str(file.loc['Minimum Binned Scans for Features','Parameters']))
            self.RegressRTWindow.setText(str(file.loc['Retention Time Window For Peak Height','Parameters']))
            self.RegressMassWindow.setText(str(file.loc['Regression Mass Window','Parameters']))
            self.RegressionMinR2.setText(str(file.loc['Minimum R2','Parameters']))
            self.DSPPMError.setText(str(file.loc['PPM Error','Parameters']))
            self.DSTargetFormulaFile.setText(str(file.loc['TargetFormulas','Parameters']))
            removeindex = self.RegressKeepPlots.findText(file.loc['RemovePlots','Parameters'])
            self.RegressKeepPlots.setCurrentIndex(removeindex)
            fupaindex = self.DSFullPartial.findText(file.loc['FullyorPartially','Parameters'])
            self.DSFullPartial.setCurrentIndex(fupaindex)
            tarlistindex = self.DSTargetList.findText(file.loc['FormulaListFile','Parameters'])
            self.DSTargetList.setCurrentIndex(tarlistindex)
            self.PremadeList.setText(str(file.loc['PremadeList','Parameters']))
            preferindex = self.DSTargetList_2.findText(file.loc['PreferCHNO','Parameters'])
            self.DSTargetList_2.setCurrentIndex(preferindex)
            self.Message.setText('Loaded')
            
class SequenceWindow(QDialog, Ui_Dialog):
    def __init__(self, parent=None):
        super().__init__()
        self.setupUi(self)
        self.StartSequence.setText('Confirm')
        self.StartSequence.clicked.connect(self.update)
    def update(self):
        self.label.setText('Click to Start\nCheck Python Console for Updates\nPlease Leave Window Open')
        self.label.setAlignment(QtCore.Qt.AlignCenter)
        self.StartSequence.setText('Start Sequence')
        self.StartSequence.clicked.connect(self.start)
    def start(self):
        SequenceRun()
        
        self.label.setText('Analysis Finished! Check Output Folder')
        self.StartSequence.setText('Close Dialog')
        self.StartSequence.clicked.disconnect()
        self.StartSequence.clicked.connect(self.close)
        return

class RunWindow(QDialog, Ui_Dialog2):
    def __init__(self, parent=None):
        super().__init__()
        self.setupUi(self)
        self.Close.clicked.connect(self.RUN_clickedWindow)
    def RUN_clickedWindow(self):
        
        win.RUN_clicked()
        
        win.Message.setText('Analysis Finished! Check Output Folder')
        self.close()
        return
        
class NumericalDecimal(QValidator):
    def validate(self, string, index):
        numbers = re.compile('[\d\.]+')
        
        if string == '':
            return QValidator.State.Acceptable, string, index
        if numbers.fullmatch(string):
            return QValidator.State.Acceptable, string, index
        else:
            return QValidator.State.Invalid, string, index
            
class Numerical(QValidator):
    def validate(self, string, index):
        numbers = re.compile('[\d]+')
        
        if string == '':
            return QValidator.State.Acceptable, string, index
        if numbers.fullmatch(string):
            return QValidator.State.Acceptable, string, index
        else:
            return QValidator.State.Invalid, string, index       
        
def checklabelinput(num):
    try:
        float(num)
        return True
    except ValueError:
        if num == '13C' or num == '15N' or num == '2H' or num == '18O' or num == '34S':
            return True
        else:
            return False
        
if __name__ == "__main__":
    mp.freeze_support()
    app = QApplication(sys.argv)
    win = Window()
    win.show()
    sys.exit(app.exec())
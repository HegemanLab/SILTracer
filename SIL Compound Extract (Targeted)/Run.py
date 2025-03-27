# -*- coding: utf-8 -*-
"""
Created on Tue Jan 21 13:30:48 2025

@author: Evan Larson, Hegeman & Cohen Labs - University of Minnesota
"""
import sys, pandas as pd, time
from PyQt5.QtWidgets import QMainWindow, QApplication, QDialog, QFileDialog
from Scripts.GUI import Ui_MainWindow
from Scripts.SILRegressionV11 import SILCmpdExtract, removepkls, makepkls


class Window (QMainWindow, Ui_MainWindow):
    def __init__(self, parent=None):
        super().__init__()
        self.setupUi(self)
        self.SelectFile.clicked.connect(self.selectfile)
        self.Run.clicked.connect(self.RunAnalysis)
    def selectfile(self):
        filename = QFileDialog.getOpenFileName(self, 'Load Masslist File', '.\\', "Excel File (*.xlsx)")
        self.Masslist.setText(filename[0])
    def RunAnalysis(self):
        masslistlink = str(self.Masslist.text())
        try:
            masslist = pd.read_excel(str(masslistlink))
        except FileNotFoundError:
            self.Message.setText('Incorrect Masslist Path')
            return
        try:
            for tup in masslist.itertuples():
                regressstats ={}
                scantype = tup.ScanType
                Compound = tup.Cmpds
                StdMass = tup.Std
                StdMass2 = tup.Std2
                SILMass = tup.SIL
                SILMass2 = tup.SIL2
                SILpre = tup.SILPre
                Stdpre = tup.StdPre
                RTMin = tup.StdRtMin
                RTMax = tup.StdRtMax
                SILRTMin = tup.SILRtMin
                SILRTMax = tup.SILRtMax
                masswindow = tup.MassWindow
                if '_' in Compound:
                   cmpdshort = Compound.split('_',1)[0]
                else:
                    cmpdshort = Compound
                RTAdj = tup.RTAdj
        except AttributeError:
            self.Message.setText('Incorrect Masslist Format')
            return
        if self.Replicates.currentText() == 'Yes':
            replicates = True
        else:
            replicates = False
        start = time.time()
        self.Message.setText('Running Analysis, Please Wait')
        self.Message.repaint()
        makepkls()
        self.Message.setText('Running Analysis, Please Wait')
        self.Message.repaint()
        count = SILCmpdExtract(masslistlink = masslistlink, replicates = replicates)
        self.Message.setText('Removing Pickle Files')
        self.Message.repaint()
        removepkls()
        end = time.time()
        self.Message.setText("\n\nTook {} minutes to analyze {} files.".format(str(round(((end-start)/60),1)), str(count)))
        self.Message.repaint()
        return
        
            


if __name__ == "__main__":
    app = QApplication(sys.argv)
    win = Window()
    win.show()
    sys.exit(app.exec())
#Attributions:
#Integration Matplotlib to QTDesigner: https://yapayzekalabs.blogspot.com/2018/11/pyqt5-gui-qt-designer-matplotlib.html
#PyQt6 and QTDesigner tutotiral: https://www.pythonguis.com/pyqt6/
#Icons: the japanese guy

#General modules
import sys
import datetime
import numpy as np
import traceback

#Module to deal with files from the system
import os

#Modules for GUI
from PyQt6 import QtWidgets, uic
from PyQt6.QtCore import Qt, QRunnable, QThreadPool, pyqtSlot, pyqtSignal, QObject
from PyQt6 import QtWidgets
from PyQt6.QtWidgets import QApplication, QWidget, QCheckBox, QFileDialog, QDialog
from PyQt6.QtGui import QDoubleValidator, QIntValidator

#Modules for plotting
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.figure import Figure
from matplotlib.ticker import MultipleLocator
import matplotlib.pyplot as plt
plt.ioff()
import matplotlib
matplotlib.use('agg')

import warnings
#warnings.simplefilter("ignore", UserWarning)

#Modules for modal analysis and EMM-AMR processing
import CESSIPy_modEMMARM as SSI 
from MRPy import MRPy #Library with modal analysis functions
import auxFunctions_postProcessEMMARM as auxEMMARM #Library with auxiliary functions for EMM-ARM post processing

#Temporary variable
debugActivatedForDevelopment=False

#Variable to make all addresses relative
basedir = os.path.dirname(__file__)

#Make windows recognize this app in the correct group when it is distributed in an installer
try:
    from ctypes import windll  # Only exists on Windows.
    myappid = 'com.EMMARMpostprocess.v.0.1'
    windll.shell32.SetCurrentProcessExplicitAppUserModelID(myappid)
except ImportError:
    pass

class eModulus_mainWindow(QtWidgets.QMainWindow):

    def __init__(self, *args, **kwargs):
        super(eModulus_mainWindow, self).__init__(*args, **kwargs)

        #Load the UI Page
        uic.loadUi(os.path.join(basedir,'UI/eModulus_mainwindow.ui'), self)
        
        #Define global properties with default values
        #These values may be used as reference for validating user inputs
        self.pathForFile=r''
        self.testDescription="Test"
        self.ageAtBeginningOfTest=30
        self.massOfEmptyTube=67.4/1000
        self.internalTubeDiameter=0.016814
        self.externalTubeDiameter=0.02014
        self.totalLengthOfTube=0.55
        self.freeCantileverLength=0.45
        self.massOfFilledTube=0.28632
        self.massAtTip=14.9/1000
        self.freeCantileverLengthEmptyTube=0.45
        self.massAtTipEmptyTube=7.8/1000
        self.frequencyEmptyTube=25.88473
        self.dormantAgeThreshold=None

        #Computed variables
        self.openModalApp=None
        self.modulusElasticity=None
        self.ages=None
        self.tubeMomentOfInertia=None
        self.materialMomentOfInertia=None
        self.tubeEmptyLinearMass=None
        self.tubeEmptyMassFreeLength=None
        self.tubeModulusInitialGuess=None
        self.tubeModulus=None
        self.tubeFullLinearMass=None
        self.tubeFullMassFreeLength=None
        self.calculationDetails={
                'Mould section inertia (m4)':None,
                'Material section inertia (m4)':None, 
                'Mould linear mass (kg/m)': None, 
                'Total mass of empty mould free length (kg)': None,
                'Initial guess of mould modulus (MPa)': None,
                'Mould modulus (MPa)': None, 
                'Frequency of empty tube (Hz)': None,  
                'Filled mould linear masss (kg/m)': None, 
                'Total mass of filled mould free length (kg/m)': None}

        #Define threadpool for multithreading
        self.threadpool = QThreadPool()

        #Define the runtime configuration
        self.currentAnalysisConfiguration={}
        self.updateCurrentAnalysisConfiguration()

        #Update display values in the GUI with the current values of the parameters
        self.updateGUIDisplayValues()

        #Define validators
        self.lineEdit_ageAtBeginningOfTest.setValidator(QDoubleValidator())
        self.lineEdit_massEmptyMould.setValidator(QDoubleValidator())
        self.lineEdit_internalMouldDiameter.setValidator(QDoubleValidator())
        self.lineEdit_externalMouldDiameter.setValidator(QDoubleValidator())
        self.lineEdit_totalLengthMould.setValidator(QDoubleValidator())
        self.lineEdit_freeLengthAtTest.setValidator(QDoubleValidator())
        self.lineEdit_totalMassFilledMould.setValidator(QDoubleValidator())
        self.lineEdit_massAtTip.setValidator(QDoubleValidator())
        self.lineEdit_freeLengthEmptyTest.setValidator(QDoubleValidator())
        self.lineEdit_massTipEmptyTest.setValidator(QDoubleValidator())
        self.lineEdit_averageFrequencyEmptyTest.setValidator(QDoubleValidator())
        self.lineEdit_dormantThreshold.setValidator(QDoubleValidator())
        
        #Define signals
        # Sends the current index (position) of the selected item.
        self.pushButton_openFilePathDialog.clicked.connect(self.openFilePathDialog)
        self.lineEdit_filePath.textEdited.connect(self.update_pathForFile)
        self.lineEdit_testDescription.textEdited.connect(self.update_testDescription)
        self.lineEdit_ageAtBeginningOfTest.textChanged.connect(self.update_ageAtBeginningOfTest)
        self.lineEdit_massEmptyMould.textChanged.connect(self.update_massOfEmptyTube)
        self.lineEdit_internalMouldDiameter.textChanged.connect(self.update_internalTubeDiameter)
        self.lineEdit_externalMouldDiameter.textChanged.connect(self.update_externalTubeDiameter)
        self.lineEdit_totalLengthMould.textChanged.connect(self.update_totalLengthOfTube)
        self.lineEdit_freeLengthAtTest.textChanged.connect(self.update_freeCantileverLength)
        self.lineEdit_totalMassFilledMould.textChanged.connect(self.update_massOfFilledTube)
        self.lineEdit_massAtTip.textChanged.connect(self.update_massAtTip)
        self.lineEdit_freeLengthEmptyTest.textChanged.connect(self.update_freeCantileverLengthEmptyTube)
        self.lineEdit_massTipEmptyTest.textChanged.connect(self.update_massAtTipEmptyTube)
        self.lineEdit_averageFrequencyEmptyTest.textChanged.connect(self.update_frequencyEmptyTube)
        self.lineEdit_dormantThreshold.textChanged.connect(self.update_dormantAgeThreshold)
        self.checkBox_mouldEModulusEstimation.stateChanged.connect(self.update_useDormantAge)

        ##Run analysis
        self.pushButton_runAnalysis.clicked.connect(self.runAnalysis)

        ##Frequency evolution tab
        self.pushButton_openFilePathModulusEvolution.clicked.connect(self.openModulusEvolutionPathDialog)
        self.lineEdit_filePathModulusEvolution.textChanged.connect(self.updateModulusEvolutionPath)
        self.toolButton_exportModulusEvolution.clicked.connect(self.exportModulusEvolution)

        ##Menu bar
        self.action_exportAnalysisConfiguration.triggered.connect(self.exportAnalysisConfiguration)
        self.action_importAnalysisConfiguration.triggered.connect(self.importAnalysisConfiguration)
        self.action_openModalAnalysisApp.triggered.connect(self.openModalAnalysisApp)
        self.action_close.triggered.connect(self.closeApp)
        self.action_close.triggered.connect(self.exportAnalysisConfiguration)
        self.actionAbout.triggered.connect(self.aboutWindow)

    #Definition of methods to perform with the signals
    ##Run analysis
    def runAnalysis(self):
        worker = Worker(self)
        worker.signals.updateGraph.connect(self.updateGraphs)
        worker.signals.updateProgressBar.connect(self.updateProgressBar)
        worker.signals.updateLogMessage.connect(self.updateLogMessage)
        # Execute
        self.threadpool.start(worker)

    ##Menu panel
    def update_pathForFile(self, new_value):
        #Only updates if new_value is different than current value
        if self.pathForFile != new_value:
            self.pathForFile = new_value
            print(self.pathForFile) if debugActivatedForDevelopment else None 
    def openFilePathDialog(self):
        #Create a QFileDialog instance
        file_dialog = QFileDialog(self)
        # Show the file dialog and get the selected file
        selected_file_path, selected_filter = file_dialog.getOpenFileName(self, "Open File", "",  "Text files (*.txt);;All Files (*)")
        #Verify the user has selected an appropriate file extension
        if selected_file_path and selected_file_path.endswith("txt"):
            self.pathForFile = selected_file_path
            self.lineEdit_filePath.setText(selected_file_path)
            print("Selected file:", selected_file_path) if debugActivatedForDevelopment else None   
        else:
            #TODO: implement textBrowser_log messages instead of print
            message = ">{0}: Error opening file: please select a .txt file.".format(datetime.datetime.now())
            self.textBrowser_log.append(message)
            self.pathForFile = None
    def update_testDescription(self, new_value):
        try:
            self.testDescription=new_value
        except Exception as e:
            #Expection to deal with when the field is
            self.testDescription=None
        print(self.testDescription) if debugActivatedForDevelopment else None 
    def update_ageAtBeginningOfTest(self, new_value):
        try:
            self.ageAtBeginningOfTest=self.try_float(new_value)
        except Exception as e:
            #Expection to deal with when the field is
            self.ageAtBeginningOfTest=None
        print(self.ageAtBeginningOfTest) if debugActivatedForDevelopment else None   
    def update_massOfEmptyTube(self, new_value):
        try:
            self.massOfEmptyTube=self.try_float(new_value)
        except Exception as e:
            #Expection to deal with when the field is
            self.massOfEmptyTube=None
        print(self.massOfEmptyTube) if debugActivatedForDevelopment else None   
    def update_internalTubeDiameter(self, new_value):
        try:
            self.internalTubeDiameter=self.try_float(new_value)
        except Exception as e:
            #Expection to deal with when the field is
            self.internalTubeDiameter=None
        print(self.internalTubeDiameter) if debugActivatedForDevelopment else None   
    def update_externalTubeDiameter(self, new_value):
        try:
            self.externalTubeDiameter=self.try_float(new_value)
        except Exception as e:
            #Expection to deal with when the field is
            self.externalTubeDiameter=None
        print(self.externalTubeDiameter) if debugActivatedForDevelopment else None   
    def update_totalLengthOfTube(self, new_value):
        try:
            self.totalLengthOfTube=self.try_float(new_value)
        except Exception as e:
            #Expection to deal with when the field is
            self.totalLengthOfTube=None
        print(self.totalLengthOfTube) if debugActivatedForDevelopment else None   
    def update_freeCantileverLength(self, new_value):
        try:
            self.freeCantileverLength=self.try_float(new_value)
        except Exception as e:
            #Expection to deal with when the field is
            self.freeCantileverLength=None
        print(self.freeCantileverLength) if debugActivatedForDevelopment else None   
    def update_massOfFilledTube(self, new_value):
        try:
            self.massOfFilledTube=self.try_float(new_value)
        except Exception as e:
            #Expection to deal with when the field is
            self.massOfFilledTube=None
        print(self.massOfFilledTube) if debugActivatedForDevelopment else None   
    def update_massAtTip(self, new_value):
        try:
            self.massAtTip=self.try_float(new_value)
        except Exception as e:
            #Expection to deal with when the field is
            self.massAtTip=None
        print(self.massAtTip) if debugActivatedForDevelopment else None 
    def update_freeCantileverLengthEmptyTube(self, new_value):
        try:
            self.freeCantileverLengthEmptyTube=self.try_float(new_value)
        except Exception as e:
            #Expection to deal with when the field is
            self.freeCantileverLengthEmptyTube=None
        print(self.freeCantileverLengthEmptyTube) if debugActivatedForDevelopment else None 
    def update_massAtTipEmptyTube(self, new_value):
        try:
            self.massAtTipEmptyTube=self.try_float(new_value)
        except Exception as e:
            #Expection to deal with when the field is
            self.massAtTipEmptyTube=None
        print(self.massAtTipEmptyTube) if debugActivatedForDevelopment else None 
    def update_frequencyEmptyTube(self, new_value):
        try:
            self.frequencyEmptyTube=self.try_float(new_value)
        except Exception as e:
            #Expection to deal with when the field is
            self.frequencyEmptyTube=None
        print(self.frequencyEmptyTube) if debugActivatedForDevelopment else None 
    def update_dormantAgeThreshold(self, new_value):
        try:
            self.dormantAgeThreshold=self.try_float(new_value)
        except Exception as e:
            #Expection to deal with when the field is
            self.dormantAgeThreshold=None
        print(self.dormantAgeThreshold) if debugActivatedForDevelopment else None 
    def update_useDormantAge(self, state):
        if state == Qt.CheckState.Checked.value:
            self.lineEdit_dormantThreshold.setEnabled(True)
            self.label_dormantThreshold.setEnabled(True)
            self.label_dormantThresholdValue.setEnabled(True)
        else:
            self.dormantAgeThreshold=None
            self.lineEdit_dormantThreshold.setEnabled(False)
            self.label_dormantThreshold.setEnabled(False)
            self.label_dormantThresholdValue.setEnabled(False)
    def updateProgressBar(self, list):
        progress, typeOfInformation = list
        if typeOfInformation == 'percentage':
            self.progressBar.setFormat('%p%')
            self.progressBar.setValue(progress)
        else:
            self.progressBar.setFormat(typeOfInformation)

    ##Modulus evolution panel
    def openModulusEvolutionPathDialog(self):
        #Create a QFileDialog instance
        file_dialog = QFileDialog(self)
        selected_folder_path = file_dialog.getExistingDirectory(self, "Select Folder", "")
        # Process the selected folder (e.g., print the folder path)
        if selected_folder_path:
            self.lineEdit_filePathModulusEvolution.setText(selected_folder_path)
            print(selected_folder_path) if debugActivatedForDevelopment else None 
    def updateModulusEvolutionPath(self, new_value):
        if os.path.exists(new_value):
            self.toolButton_exportModulusEvolution.setEnabled(True)
        else:
            message = ">{0}: Folder does not exist. Please select an existing folder.".format(datetime.datetime.now())
            self.textBrowser_log.append(message)
            self.toolButton_exportModulusEvolution.setEnabled(False)
    def exportModulusEvolution(self):
        saveFilePrefix=self.lineEdit_filePathModulusEvolution.text()+'/'+self.lineEdit_fileNameModulusEvolution.text()
        np.savetxt(saveFilePrefix+'.txt', np.vstack((self.ages, self.modulusElasticity.T)).T, delimiter='\t', fmt='%f', header="=======================\nTest description\n"+self.testDescription+"\n=======================\nAGE(SECONDS)\tELASTIC MODULUS(PA)")
        logMessage = ">{0}: E-modulus file saved succesfuly.".format(datetime.datetime.now())
        self.updateLogMessage(logMessage)
        
    ##Toolbar
    def exportAnalysisConfiguration(self,s):
        #Create a QFileDialog instance
        file_dialog = QFileDialog(self)
        selected_file_path, selected_filter = file_dialog.getSaveFileName(self, "Enter name", "", "Text file (*.txt)")
        self.updateCurrentAnalysisConfiguration()
        import json
        with open(selected_file_path, 'w') as file:
            json.dump(self.currentAnalysisConfiguration, file)  # encode dict into JSON
        logMessage = ">{0}: Export analysis configuration succesful.".format(datetime.datetime.now())
        self.updateLogMessage(logMessage)
    def importAnalysisConfiguration(self,s):
        #Create a QFileDialog instance
        file_dialog = QFileDialog(self)
        selected_file_path, selected_filter = file_dialog.getOpenFileName(self, "Enter name", "", "Text file (*.txt)")
        import json
        with open(selected_file_path, 'r') as file:
            self.currentAnalysisConfiguration=json.load(file)  # encode dict into JSON
        self.readCurrentAnalysisConfiguration()
        self.updateGUIDisplayValues()
        logMessage = ">{0}: Import analysis configuration succesful.".format(datetime.datetime.now())
        self.updateLogMessage(logMessage)
    def openModalAnalysisApp(self,s):
        self.openModalApp=True
        self.close()
    def closeApp(self,s):
        self.close()
    def aboutWindow(self,s):
        dlg_about = aboutDlg(self)
        dlg_about.exec()

    ##Plot graph
    def updateGraphs(self, plotToUpdate=None):
        #Plot to update is an int from 0...N
        if plotToUpdate == 0:
            try:
                self.figModulus = Figure()
                self.ax = self.figModulus.add_subplot(111)
                self.ax.plot(self.ages/(60*60), self.modulusElasticity/(1e9), label="E-modulus evolution")
                self.ax.set_ylabel("E-modulus (GPa)")
                self.ax.set_xlabel("Age (hours)")
                self.ax.set_xscale('log')
                #plt.xlim([1e-2,22])
                self.ax.set_ylim((-0.2, None))
                self.ax.grid(which='both', axis='both', linestyle='-', color='whitesmoke') 
                self.ax.set_axisbelow(True)
                self.ax.legend(loc='lower right')
                self.figModulus.set_tight_layout(True)
                self.graph_modulusEvolution.update_figure(self.figModulus)
            except Exception as e:
                logMessage = ">{0}: An exception occurred. This may be due to bad input. Exception description: \n {1}".format(datetime.datetime.now(),traceback.format_exc())
                self.updateLogMessage(logMessage)
    ##Complementary functions
    def updateLogMessage(self, message):
        self.textBrowser_log.append(message)
    def try_float(self,v):
        #Convert a string to float and handle the case the string is empty, converting it to None
        try:
            return float(v)
        except Exception:
            return None
    def updateCurrentAnalysisConfiguration(self):
        self.currentAnalysisConfiguration['pathForFile']=self.pathForFile
        self.currentAnalysisConfiguration['testDescription']=self.testDescription
        self.currentAnalysisConfiguration['ageAtBeginningOfTest']=self.ageAtBeginningOfTest
        self.currentAnalysisConfiguration['massOfEmptyTube']=self.massOfEmptyTube
        self.currentAnalysisConfiguration['internalTubeDiameter']=self.internalTubeDiameter        
        self.currentAnalysisConfiguration['externalTubeDiameter']=self.externalTubeDiameter
        self.currentAnalysisConfiguration['totalLengthOfTube']=self.totalLengthOfTube
        self.currentAnalysisConfiguration['freeCantileverLength']=self.freeCantileverLength
        self.currentAnalysisConfiguration['massOfFilledTube']=self.massOfFilledTube
        self.currentAnalysisConfiguration['massAtTip']=self.massAtTip
        self.currentAnalysisConfiguration['freeCantileverLengthEmptyTube']=self.freeCantileverLengthEmptyTube
        self.currentAnalysisConfiguration['massAtTipEmptyTube']=self.massAtTipEmptyTube
        self.currentAnalysisConfiguration['frequencyEmptyTube']=self.frequencyEmptyTube
        self.currentAnalysisConfiguration['dormantAgeThreshold']=self.dormantAgeThreshold
    def readCurrentAnalysisConfiguration(self):
        self.pathForFile=self.currentAnalysisConfiguration['pathForFile']
        self.testDescription=self.currentAnalysisConfiguration['testDescription']
        self.ageAtBeginningOfTest=self.currentAnalysisConfiguration['ageAtBeginningOfTest']
        self.massOfEmptyTube=self.currentAnalysisConfiguration['massOfEmptyTube']
        self.internalTubeDiameter=self.currentAnalysisConfiguration['internalTubeDiameter']
        self.externalTubeDiameter=self.currentAnalysisConfiguration['externalTubeDiameter']
        self.totalLengthOfTube=self.currentAnalysisConfiguration['totalLengthOfTube']
        self.freeCantileverLength=self.currentAnalysisConfiguration['freeCantileverLength']
        self.massOfFilledTube=self.currentAnalysisConfiguration['massOfFilledTube']
        self.massAtTip=self.currentAnalysisConfiguration['massAtTip']
        self.freeCantileverLengthEmptyTube=self.currentAnalysisConfiguration['freeCantileverLengthEmptyTube']
        self.massAtTipEmptyTube=self.currentAnalysisConfiguration['massAtTipEmptyTube']
        self.frequencyEmptyTube=self.currentAnalysisConfiguration['frequencyEmptyTube']
        self.dormantAgeThreshold=self.currentAnalysisConfiguration['dormantAgeThreshold']
    def updateGUIDisplayValues(self): 
        self.lineEdit_filePath.setText(str(self.pathForFile))
        self.lineEdit_testDescription.setText(str(self.testDescription))
        self.lineEdit_ageAtBeginningOfTest.setText(str(self.ageAtBeginningOfTest))
        self.lineEdit_massEmptyMould.setText(str(self.massOfEmptyTube))
        self.lineEdit_internalMouldDiameter.setText(str(self.internalTubeDiameter))
        self.lineEdit_externalMouldDiameter.setText(str(self.externalTubeDiameter))
        self.lineEdit_totalLengthMould.setText(str(self.totalLengthOfTube))
        self.lineEdit_freeLengthAtTest.setText(str(self.freeCantileverLength))
        self.lineEdit_totalMassFilledMould.setText(str(self.massOfFilledTube))
        self.lineEdit_massAtTip.setText(str(self.massAtTip))
        self.lineEdit_freeLengthEmptyTest.setText(str(self.freeCantileverLengthEmptyTube))
        self.lineEdit_massTipEmptyTest.setText(str(self.massAtTipEmptyTube))
        self.lineEdit_averageFrequencyEmptyTest.setText(str(self.frequencyEmptyTube))
        self.lineEdit_dormantThreshold.setText(str(self.dormantAgeThreshold))

class aboutDlg(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        #Load the dialog's GUI
        uic.loadUi(os.path.join(basedir,'UI/aboutDialog.ui'), self)

class WorkerSignals(QObject):
    updateGraph = pyqtSignal(int)
    updateProgressBar = pyqtSignal(list)
    updateLogMessage = pyqtSignal(str)

class Worker(QRunnable):
    '''
    Worker thread

    Inherits from QRunnable to handler worker thread setup, signals and wrap-up.

    :param callback: The function callback to run on this worker thread. Supplied args and
                     kwargs will be passed through to the runner.
    :type callback: function
    :param args: Arguments to pass to the callback function
    :param kwargs: Keywords to pass to the callback function
    '''

    def __init__(self, analysisObject):
        super(Worker, self).__init__()
        '''
        updatePlotInBatchAnalysis: bool
            This variable is used in the context of the "Update plot" button in Batch Analysis
        '''
        # Store constructor arguments (re-used for processing)
        self.analysisObject = analysisObject
        self.signals = WorkerSignals()

    @pyqtSlot()
    def run(self):       
        #Perform the competent analysis
        try:
            self.signals.updateProgressBar.emit([0,'percentage'])
            #Open file with vibration frequnecies of the EMM-ARM test
            modalIdentificationData = np.loadtxt(self.analysisObject.pathForFile, dtype='float', delimiter='\t')
            #Extract only the first mode vibration frequencies
            vibrationFrequencies = modalIdentificationData[:,1]
            #Create the vector that will store the predicted composite flexural stiffness
            self.analysisObject.modulusElasticity = np.zeros(len(vibrationFrequencies))
            #Correct ages accordingly to the instant of beggining of the test.
            #Result file is always in seconds
            self.analysisObject.ages = modalIdentificationData[:,0]+self.analysisObject.ageAtBeginningOfTest*60 

            #Compute geometrical and mechanical properties of the mould tube
            self.analysisObject.tubeMomentOfInertia = np.pi*((self.analysisObject.externalTubeDiameter**4)-(self.analysisObject.internalTubeDiameter**4))/64
            self.analysisObject.materialMomentOfInertia = np.pi*((self.analysisObject.internalTubeDiameter**4))/64
            self.analysisObject.tubeEmptyLinearMass = self.analysisObject.massOfEmptyTube/self.analysisObject.totalLengthOfTube
            self.analysisObject.tubeEmptyMassFreeLength = self.analysisObject.freeCantileverLengthEmptyTube*self.analysisObject.tubeEmptyLinearMass
            self.signals.updateProgressBar.emit([5,'percentage'])
            #Estimate the E-modulus of the mould tube from the transcental equation
            if str(self.analysisObject.comboBox_typeOfTest.currentText())=="Cantilever":
                self.analysisObject.tubeModulusInitialGuess = (self.analysisObject.freeCantileverLengthEmptyTube**3)*(self.analysisObject.massAtTipEmptyTube+0.24*self.analysisObject.tubeEmptyMassFreeLength)*((self.analysisObject.frequencyEmptyTube*2*np.pi)**2)/(3*self.analysisObject.tubeMomentOfInertia) #Modulus given in Pa
                self.analysisObject.tubeFlexuralStiffness = auxEMMARM.solveTranscendentalEquation(self.analysisObject.tubeModulusInitialGuess*self.analysisObject.tubeMomentOfInertia, self.analysisObject.frequencyEmptyTube, self.analysisObject.tubeEmptyLinearMass, self.analysisObject.freeCantileverLengthEmptyTube, self.analysisObject.massAtTipEmptyTube)
                self.analysisObject.tubeModulus = self.analysisObject.tubeFlexuralStiffness/self.analysisObject.tubeMomentOfInertia
            else:
                #It is simply supported, as the value comes from a combobox with only 2 options
                self.analysisObject.tubeModulusInitialGuess = (self.analysisObject.freeCantileverLengthEmptyTube**3)*(self.analysisObject.massAtTipEmptyTube+0.49*self.analysisObject.tubeEmptyMassFreeLength)*((self.analysisObject.frequencyEmptyTube*np.pi)**2)/(12*self.analysisObject.tubeMomentOfInertia) #Modulus given in Pa
                self.analysisObject.tubeFlexuralStiffness = auxEMMARM.solveTranscendentalEquation(self.analysisObject.tubeModulusInitialGuess*self.analysisObject.tubeMomentOfInertia, self.analysisObject.frequencyEmptyTube, self.analysisObject.tubeEmptyLinearMass, self.analysisObject.freeCantileverLengthEmptyTube, self.analysisObject.massAtTipEmptyTube,typeOfEquation="Simply supported" )
                self.analysisObject.tubeModulus = self.analysisObject.tubeFlexuralStiffness/self.analysisObject.tubeMomentOfInertia
            #Properties of filled tube
            self.analysisObject.tubeFullLinearMass = self.analysisObject.massOfFilledTube/self.analysisObject.totalLengthOfTube
            self.analysisObject.tubeFullMassFreeLength = self.analysisObject.freeCantileverLength*self.analysisObject.tubeFullLinearMass

            #Store all parameters calculated for the current tube to save in a checking file for later validation of the calculations performed
            self.analysisObject.calculationDetails={
                'Mould section inertia (m4)':self.analysisObject.tubeMomentOfInertia,
                'Material section inertia (m4)':self.analysisObject.materialMomentOfInertia, 
                'Mould linear mass (kg/m)': self.analysisObject.tubeEmptyLinearMass, 
                'Total mass of empty mould free length (kg)': self.analysisObject.tubeEmptyMassFreeLength,
                'Initial guess of mould modulus (GPa)':self.analysisObject.tubeModulusInitialGuess/1e9,
                'Mould modulus (GPa)':self.analysisObject.tubeModulus.tolist()[0]/1e9, 
                'Frequency of empty tube (Hz)': self.analysisObject.frequencyEmptyTube,  
                'Filled mould linear masss (kg/m)': self.analysisObject.tubeFullLinearMass, 
                'Total mass of filled mould free length (kg/m)': self.analysisObject.tubeFullMassFreeLength}

            self.signals.updateProgressBar.emit([10,'percentage'])
            progressStep = np.linspace(10+100/len(vibrationFrequencies),90,len(vibrationFrequencies))
            lastProgress = 10
            #Estimate flexural stiffness associated to every vibration frequency identified
            firstPoint=True
            for i, frequency in enumerate(vibrationFrequencies):
                currentProgressStep = int(progressStep[i])-lastProgress
                if frequency != 0:
                    if (firstPoint is True) and (self.analysisObject.dormantAgeThreshold is not None):
                        #If it is the first iteration and user set the parameter 'dormantAgeThreshold' it is because
                        #we need to find the proper EI that will make the first dormant times, obtained by visual inspection, equal to zero
                        #Module to solve the transcendental function to estimate E-modulus
                        self.signals.updateLogMessage.emit("Dormant age activated")
                        tubeFlexuralStiffnessFirstGuess = self.analysisObject.tubeFlexuralStiffness
                        self.signals.updateProgressBar.emit([int(lastProgress+currentProgressStep/3),'percentage'])
                        self.analysisObject.tubeFlexuralStiffness=auxEMMARM.findTubeFlexuralStiffnessFromDormantAge(tubeFlexuralStiffnessFirstGuess, self.analysisObject.freeCantileverLength, self.analysisObject.massAtTip, 
                                                                                                                    self.analysisObject.tubeFullMassFreeLength, self.analysisObject.tubeFullLinearMass, self.analysisObject.ages, 
                                                                                                                    vibrationFrequencies, self.analysisObject.dormantAgeThreshold, typeOfEquation=str(self.analysisObject.comboBox_typeOfTest.currentText()))
                        #Estimate an initial guess based on a simplified equation
                        if str(self.analysisObject.comboBox_typeOfTest.currentText())=="Cantilever":
                            #Initial guess from Blevins book, pag 158
                            compositeFlexuralStiffnessInitialGuess = ((self.analysisObject.freeCantileverLength)**3)*(self.analysisObject.massAtTip+0.24*self.analysisObject.tubeFullMassFreeLength)*((frequency*2*np.pi)**2)/(3) #Modulus given in Pa
                            self.signals.updateProgressBar.emit([int(lastProgress+currentProgressStep/2),'percentage'])
                            compositeFlexuralStiffness = auxEMMARM.solveTranscendentalEquation(compositeFlexuralStiffnessInitialGuess, frequency, self.analysisObject.tubeFullLinearMass, self.analysisObject.freeCantileverLength, self.analysisObject.massAtTip)
                        else:
                            #Initial guess from Blevins book, pag 158
                            compositeFlexuralStiffnessInitialGuess = ((self.analysisObject.freeCantileverLength)**3)*(self.analysisObject.massAtTip+0.49*self.analysisObject.tubeFullMassFreeLength)*((frequency*np.pi)**2)/(12) #Modulus given in Pa
                            self.signals.updateProgressBar.emit([int(lastProgress+currentProgressStep/2),'percentage'])
                            compositeFlexuralStiffness = auxEMMARM.solveTranscendentalEquation(compositeFlexuralStiffnessInitialGuess, frequency, self.analysisObject.tubeFullLinearMass, self.analysisObject.freeCantileverLength, self.analysisObject.massAtTip, typeOfEquation="Simply supported")
                        self.analysisObject.modulusElasticity[i] = (compositeFlexuralStiffness-self.analysisObject.tubeFlexuralStiffness)/self.analysisObject.materialMomentOfInertia
                        firstPoint = False
                    else:
                        #If the above is not true, then we simply need to use the tube flexural stiffness estimated from the empty tests
                        #Estimate an initial guess based on a simplified equation
                        #Initial guess from Blevins book, pag 158
                        if str(self.analysisObject.comboBox_typeOfTest.currentText())=="Cantilever":
                            compositeFlexuralStiffnessInitialGuess = ((self.analysisObject.freeCantileverLength)**3)*(self.analysisObject.massAtTip+0.24*self.analysisObject.tubeFullMassFreeLength)*((frequency*2*np.pi)**2)/(3) #Modulus given in Pa
                            self.signals.updateProgressBar.emit([int(lastProgress+currentProgressStep/2),'percentage'])
                            compositeFlexuralStiffness = auxEMMARM.solveTranscendentalEquation(compositeFlexuralStiffnessInitialGuess, frequency, self.analysisObject.tubeFullLinearMass, self.analysisObject.freeCantileverLength, self.analysisObject.massAtTip)
                        else:
                            compositeFlexuralStiffnessInitialGuess = ((self.analysisObject.freeCantileverLength)**3)*(self.analysisObject.massAtTip+0.49*self.analysisObject.tubeFullMassFreeLength)*((frequency*np.pi)**2)/(12) #Modulus given in Pa
                            self.signals.updateProgressBar.emit([int(lastProgress+currentProgressStep/2),'percentage'])
                            compositeFlexuralStiffness = auxEMMARM.solveTranscendentalEquation(compositeFlexuralStiffnessInitialGuess, frequency, self.analysisObject.tubeFullLinearMass, self.analysisObject.freeCantileverLength, self.analysisObject.massAtTip, typeOfEquation="Simply supported")
                        self.analysisObject.modulusElasticity[i] = (compositeFlexuralStiffness-self.analysisObject.tubeFlexuralStiffness)/self.analysisObject.materialMomentOfInertia
                else:
                    self.analysisObject.modulusElasticity[i] = np.nan
                    self.analysisObject.ages[i] = np.nan
                self.signals.updateProgressBar.emit([int(lastProgress+currentProgressStep),'percentage'])
                lastProgress = lastProgress+currentProgressStep
            
            self.analysisObject.modulusElasticity = np.array([elem for elem in self.analysisObject.modulusElasticity if not(np.isnan(elem))])
            self.analysisObject.ages = np.array([elem for elem in self.analysisObject.ages if not(np.isnan(elem))])
            self.signals.updateGraph.emit(0)
            self.signals.updateProgressBar.emit([100,'percentage'])
            self.signals.updateProgressBar.emit([100,'Analysis complete'])
            cDK = [key for key in self.analysisObject.calculationDetails]
            logMessage = ">{0}: Calculation details: \n {1}: {2} \n  {3}: {4}\n  {5}: {6} \n  {7}: {8} \n  {9}: {10} \n  {11}: {12} \n  {13}: {14} \n  {15}: {16} \n  {17}: {18}  ".format(datetime.datetime.now(),
                                                                                                                                                                                            cDK[0],self.analysisObject.calculationDetails[cDK[0]],
                                                                                                                                                                                            cDK[1],self.analysisObject.calculationDetails[cDK[1]],
                                                                                                                                                                                            cDK[2],self.analysisObject.calculationDetails[cDK[2]],
                                                                                                                                                                                            cDK[3],self.analysisObject.calculationDetails[cDK[3]],
                                                                                                                                                                                            cDK[4],self.analysisObject.calculationDetails[cDK[4]],
                                                                                                                                                                                            cDK[5],self.analysisObject.calculationDetails[cDK[5]],
                                                                                                                                                                                            cDK[6],self.analysisObject.calculationDetails[cDK[6]],
                                                                                                                                                                                            cDK[7],self.analysisObject.calculationDetails[cDK[7]],
                                                                                                                                                                                            cDK[8],self.analysisObject.calculationDetails[cDK[8]])
            self.signals.updateLogMessage.emit(logMessage)
        except Exception as e:
            # Print the exception error message
            self.signals.updateProgressBar.emit([100,'percentage'])
            self.signals.updateProgressBar.emit([0,'Error has occurred'])
            logMessage = ">{0}: An exception occurred. This may be due to bad input. Exception description: \n {1}".format(datetime.datetime.now(),traceback.format_exc())
            self.signals.updateLogMessage.emit(logMessage)

def main():

    app = QtWidgets.QApplication(sys.argv)
    main = eModulus_mainWindow()
    main.show()
    sys.exit(app.exec())

if __name__ == '__main__':
    main()
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
from PyQt6.QtCore import Qt, QRunnable, QThreadPool, pyqtSlot, pyqtSignal, QObject, QSignalBlocker
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

class modalAnalysis_mainWindow(QtWidgets.QMainWindow):

    def __init__(self, *args, **kwargs):
        super(modalAnalysis_mainWindow, self).__init__(*args, **kwargs)

        #Load the UI Page
        uic.loadUi(os.path.join(basedir,'UI/modalAnalysis_mainwindow.ui'), self)
        
        #Define global properties with default values
        #These values may be used as reference for validating user inputs
        self.pathForFile=r''
        self.fileVisualizedPath = None
        self.filesToRead = None #List of files to read, in case of a batch analysis
        self.selectedFileType = ('uEMMARM files', '*.emm')
        self.selectedFileExtension = self.selectedFileType[1][-5:]
        
        self.selectedSystem='uEMMARM'
        self.selectedProcessing='singleFile'
        self.desiredChannel = 1
        self.calibrationFactor = 1.0/2400.0
        self.samplingFrequencyOriginal = 860
        self.filterConfiguration = []
        self.nps = 4096
        self.modalIdentificationMethodToPerform = {'peak-picking':False, 
                                            'BFD':False, 
                                            'EFDD':False, 
                                            'SSI-COV':False}
        #Peak-picking properties
        self.intervalForAveragingHz=5
        #BFD properties
        self.fint_BFD = [0.8,1.2]
        #EFDD properties
        self.fint_EFDD = [0.8, 1.2]
        self.tint = np.array([0, 0.2])
        #SSI-COV properties
        self.numModesToBeConsidered = 1
        self.i = 60
        self.startingOrderNumber = 2
        self.endOrderNumber = 60
        self.incrementBetweenOrder = 1
        self.refs = [0]
        self.tol =  np.array((
                    [0.01,8,150],     # [allowed_variation, lower_bound, higher_bound]: limits for eigenfrequencies
                    [0.1,0.001,0.025],     # [allowed_variation, lower_bound, higher_bound]: limits for damping ratios
                    [0.02,0,1]))    # [allowed_variation, lower_bound, higher_bound]: limits for MAC
        
        #Define and initialize default values in plot configuration variable
        self.openModulusApp = None
        self.plotGeneral = None
        self.plotPSD = None
        self.plotNPSD = None
        self.plotPeakPicking=None
        self.plotBFD=None
        self.plotEFDD=None
        self.plotSSI=None
        self.plotConfiguration = None
        self.populate_plotConfiguration(firstInitialization=True)

        #Define threadpool for multithreading
        self.threadpool = QThreadPool()

        #Define the runtime configuration
        self.currentAnalysisConfiguration={}
        self.updateCurrentAnalysisConfiguration()

        #Update display values in the GUI with the current values of the parameters
        self.updateGUIDisplayValues()

        #Disable all graphical tabs but the first (time series and PSD)
        for tabIndex in np.arange(1,self.tabWidget_graphicalPanel.count()):
            self.tabWidget_graphicalPanel.setTabVisible(tabIndex,False)

        #Define validators
        self.lineEdit_accelConvFactor.setValidator(QDoubleValidator())
        self.lineEdit_samplingFrequency.setValidator(QIntValidator())
        self.lineEdit_butterworthLowPassValue.setValidator(QDoubleValidator())
        self.lineEdit_butterworthHighPassValue.setValidator(QDoubleValidator())
        self.lineEdit_windowLength_PSD.setValidator(QIntValidator())
        self.lineEdit_averagingInterval_PP.setValidator(QDoubleValidator())
        
        self.lineEdit_curveFittingFrequencyLowerBound_BFD.setValidator(QDoubleValidator())
        self.lineEdit_curveFittingFrequencyUpperBound_BFD.setValidator(QDoubleValidator())

        self.lineEdit_curveFittingFrequencyLowerBound_EFDD.setValidator(QDoubleValidator())
        self.lineEdit_curveFittingFrequencyUpperBound_EFDD.setValidator(QDoubleValidator())
        self.lineEdit_curveFittingTimeLowerBound_EFDD.setValidator(QDoubleValidator())
        self.lineEdit_curveFittingTimeUpperBound_EFDD.setValidator(QDoubleValidator())

        self.lineEdit_stabTolFrequency_allowedVariation_SSI.setValidator(QDoubleValidator())
        self.lineEdit_stabTolFrequency_lowerBound_SSI.setValidator(QDoubleValidator())
        self.lineEdit_stabTolFrequency_upperBound_SSI.setValidator(QDoubleValidator()) 
        self.lineEdit_stabTolDamping_allowedVariation_SSI.setValidator(QDoubleValidator()) 
        self.lineEdit_stabTolDamping_lowerBound_SSI.setValidator(QDoubleValidator()) 
        self.lineEdit_stabTolDamping_upperBound_SSI.setValidator(QDoubleValidator()) 
        
        #Define signals
        # Sends the current index (position) of the selected item.

        ## General info panel
        self.comboBox_typeOfSystem.currentIndexChanged.connect(self.update_selectedSystem)
        self.radioButton_singleFile.toggled.connect(self.update_selectedProcessing)
        self.radioButton_multipleFile.toggled.connect(self.update_selectedProcessing)
        self.spinBox_channelSelection.valueChanged.connect(self.update_desiredChannel)
        self.lineEdit_accelConvFactor.textChanged.connect(self.update_calibrationFactor)
        self.lineEdit_samplingFrequency.textChanged.connect(self.update_samplingFrequencyOriginal)
        self.pushButton_openFilePathDialog.clicked.connect(self.openFilePathDialog)
        #self.lineEdit_filePath.textChanged.connect(self.update_pathForFile)

        ## File visualized panel
        self.pushButton_fileBackwards.clicked.connect(self.updateFileBackwards)
        self.pushButton_fileForwards.clicked.connect(self.updateFileForwards)
        self.comboBox_fileVisualized.currentTextChanged.connect(self.update_fileVisualized)
        self.pushButton_updatePlot.clicked.connect(self.update_individualPlotInBatchAnalysis)

        ## Filtering panel
        #Detrend groupbox
        self.checkBox_detrend.stateChanged.connect(self.update_activateDetrendFilter)
        self.comboBox_detrendType.currentIndexChanged.connect(self.update_detrendType)
        #Butterworth groupbox
        self.checkBox_butterworthActivate.stateChanged.connect(self.update_activateButterworthFilter)
        self.spinBox_butterworthOrder.valueChanged.connect(self.update_butterworthOrder)
        self.lineEdit_butterworthLowPassValue.textChanged.connect(self.update_butterworthLowPassFrequency)
        self.lineEdit_butterworthHighPassValue.textChanged.connect(self.update_butterworthHighPassFrequency)
        #Decimation groupbox
        self.checkBox_decimationActivate.stateChanged.connect(self.update_activateDecimationFilter)
        self.spinBox_decimationFactor.valueChanged.connect(self.update_decimationFactor)

        ## Modal analysis methods panel
        #PSD groupbox
        self.lineEdit_windowLength_PSD.textChanged.connect(self.update_windowLengthPSD)
        self.pushButton_graphicalConfig_PSD.clicked.connect(lambda checked: self.openPlotConfigurationDialog(tab=0))
        #PP groupbox
        self.checkBox_activate_PP.stateChanged.connect(self.update_activatePP)
        self.lineEdit_averagingInterval_PP.textChanged.connect(self.update_intervalForAveragingHz)
        self.pushButton_graphicalConfig_PP.clicked.connect(lambda checked: self.openPlotConfigurationDialog(tab=1))
        #BFD groupbox
        self.checkBox_activate_BFD.stateChanged.connect(self.update_activateBFD)
        self.lineEdit_curveFittingFrequencyLowerBound_BFD.textChanged.connect(self.update_lower_fint_BFD)
        self.lineEdit_curveFittingFrequencyUpperBound_BFD.textChanged.connect(self.update_upper_fint_BFD)
        self.pushButton_graphicalConfig_BFD.clicked.connect(lambda checked: self.openPlotConfigurationDialog(tab=2))
        #EFDD groupbox
        self.checkBox_activate_EFDD.stateChanged.connect(self.update_activateEFDD)
        self.lineEdit_curveFittingFrequencyLowerBound_EFDD.textChanged.connect(self.update_lower_fint_EFDD)
        self.lineEdit_curveFittingFrequencyUpperBound_EFDD.textChanged.connect(self.update_upper_fint_EFDD)
        self.lineEdit_curveFittingTimeLowerBound_EFDD.textChanged.connect(self.update_lower_tint)
        self.lineEdit_curveFittingTimeUpperBound_EFDD.textChanged.connect(self.update_upper_tint)
        self.pushButton_graphicalConfig_EFDD.clicked.connect(lambda checked: self.openPlotConfigurationDialog(tab=3))
        #SSI-COV groupbox
        self.checkBox_activate_SSI.stateChanged.connect(self.update_activateSSI)
        self.spinBox_numModes_SSI.valueChanged.connect(self.update_numModesToBeConsidered)
        self.spinBox_numTimeLags_SSI.valueChanged.connect(self.update_i)
        self.spinBox_startingModelOrder_SSI.valueChanged.connect(self.update_startingOrderNumber)
        self.spinBox_endingModelOrder_SSI.valueChanged.connect(self.update_endOrderNumber)
        self.spinBox_orderIncrement_SSI.valueChanged.connect(self.update_incrementBetweenOrder)
        self.lineEdit_stabTolFrequency_allowedVariation_SSI.textChanged.connect(self.update_tol_frequencyVariation)
        self.lineEdit_stabTolFrequency_lowerBound_SSI.textChanged.connect(self.update_tol_frequencyLowerBound)
        self.lineEdit_stabTolFrequency_upperBound_SSI.textChanged.connect(self.update_tol_frequencyUpperBound)
        self.lineEdit_stabTolDamping_allowedVariation_SSI.textChanged.connect(self.update_tol_dampingUpperBound)
        self.lineEdit_stabTolDamping_lowerBound_SSI.textChanged.connect(self.update_tol_dampingUpperBound)
        self.lineEdit_stabTolDamping_upperBound_SSI.textChanged.connect(self.update_tol_dampingUpperBound)
        self.pushButton_graphicalConfig_SSI.clicked.connect(lambda checked: self.openPlotConfigurationDialog(tab=4))

        ##Run analysis
        self.pushButton_runAnalysis.clicked.connect(self.runAnalysis)

        ##Frequency evolution tab
        self.pushButton_openFilePathFrequencyEvolution.clicked.connect(self.openFrequencyEvolutionPathDialog)
        self.lineEdit_filePathFrequencyEvolution.textChanged.connect(self.updateFrequencyEvolutionPath)
        self.toolButton_exportFrequencyEvolution.clicked.connect(self.exportFrequencyEvolution)

        ##Damping evolution tab
        self.pushButton_openFilePathDampingEvolution.clicked.connect(self.openDampingEvolutionPathDialog)
        self.lineEdit_filePathDampingEvolution.textChanged.connect(self.updateDampingEvolutionPath)
        self.toolButton_exportDampingEvolution.clicked.connect(self.exportDampingEvolution)

        ##Menu bar
        self.action_exportAnalysisConfiguration.triggered.connect(self.exportAnalysisConfiguration)
        self.action_importAnalysisConfiguration.triggered.connect(self.importAnalysisConfiguration)
        self.action_openEModulusEstimationApp.triggered.connect(self.openModulusEstimationApp)
        self.action_close.triggered.connect(self.closeApp)
        self.actionAbout.triggered.connect(self.aboutWindow)

    #Definition of methods to perform with the signals
    ##Menu par
    def exportAnalysisConfiguration(self,s):
        #Create a QFileDialog instance
        file_dialog = QFileDialog(self)
        selected_file_path, selected_filter = file_dialog.getSaveFileName(self, "Enter name", "", "Text file (*.txt)")
        self.updateCurrentAnalysisConfiguration()
        import json
        with open(selected_file_path, 'w') as file:
            json.dump(self.currentAnalysisConfiguration, file)  # encode dict into JSON
    def importAnalysisConfiguration(self,s):
        #Create a QFileDialog instance
        try:
            file_dialog = QFileDialog(self)
            selected_file_path, selected_filter = file_dialog.getOpenFileName(self, "Enter name", "", "Text file (*.txt)")
            import json
            self.resetAnalysisConfiguration()
            with open(selected_file_path, 'r') as file:
                self.currentAnalysisConfiguration=json.load(file)  # encode dict into JSON
            self.readCurrentAnalysisConfiguration()
            self.updateGUIDisplayValues()
        except Exception as e:
            # Print the exception error message
            logMessage = ">{0}: An error has occurred. This may be due to bad input file. Exception description: \n {1}".format(datetime.datetime.now(),traceback.format_exc())
            self.textBrowser_log.append(logMessage)
    def openModulusEstimationApp(self,s):
        self.openModulusApp=True
        self.close()
    def closeApp(self,s):
        self.close()
    def aboutWindow(self,s):
        dlg_about = aboutDlg(self)
        dlg_about.exec()
    ##General info panel
    '''
    def update_pathForFile(self, new_value):
        #Only updates if new_value is different than current value
        if self.pathForFile != new_value:
            if self.selectedProcessing == 'singleFile':
                if self.pathForFile is None:
                    pass
                else:
                    self.pathForFile = new_value
            elif self.selectedProcessing == 'batchProcessing':
                self.pathForFile = new_value
                if self.pathForFile is None:
                    pass
                else:
                    try:
                        self.filesToRead = [x for x in os.listdir(self.pathForFile) if x[-4:] == self.selectedFileExtension[-4:]] #Get array with file names in the selected folder
                        if len(self.filesToRead) < 1:
                            raise #Raise an exception because no valid files could be found
                        self.comboBox_fileVisualized.clear()
                        self.comboBox_fileVisualized.addItems(self.filesToRead)
                        self.comboBox_fileVisualized.setCurrentIndex(0)
                    except:
                        message = ">{0}: Path informed does not exist or there are no EMM-ARM files there.".format(datetime.datetime.now())
                        self.textBrowser_log.append(message)
            print(self.pathForFile) if debugActivatedForDevelopment else None   
    '''
    def update_selectedSystem(self, currentIndex): # i is an int
        self.selectedSystem=self.comboBox_typeOfSystem.currentText()
        #Verify which fileType is to be allowed base on the informed system
        fileTypes = [('uEMMARM files', '*.emm'), ('Text files', '*.txt'), ('CSV files', '*.csv')]
        if self.selectedSystem == "uEMMARM":
            self.selectedFileType = fileTypes[0]
            self.selectedFileExtension = self.selectedFileType[1][-5:]
            self.lineEdit_samplingFrequency.setEnabled(False)
        elif self.selectedSystem == "National":
            self.selectedFileType = fileTypes[1]
            self.selectedFileExtension = self.selectedFileType[1][-5:]
            self.lineEdit_samplingFrequency.setEnabled(True)
        elif self.selectedSystem == "old_uEMMARM":
            self.selectedFileType = fileTypes[1]
            self.selectedFileExtension = self.selectedFileType[1][-5:]
            self.lineEdit_samplingFrequency.setEnabled(True)
        elif self.selectedSystem == "RaspberryPi":
            self.selectedFileType = fileTypes[2]
            self.selectedFileExtension = self.selectedFileType[1][-5:]
            self.lineEdit_samplingFrequency.setEnabled(True)
        print(self.selectedSystem) if debugActivatedForDevelopment else None
    def update_selectedProcessing(self):
        sender = self.sender()  # Get the sender of the signal
        if sender.isChecked():
            if (sender.text() == 'Single file'):
                self.selectedProcessing = 'singleFile'
                self.pathForFile = None
                self.lineEdit_filePath.setText("")
                self.comboBox_fileVisualized.setEnabled(False)
                self.pushButton_updatePlot.setEnabled(False)
                self.pushButton_fileBackwards.setEnabled(False)
                self.pushButton_fileForwards.setEnabled(False)
                self.tabWidget_graphicalPanel.setTabVisible(5,False)
                self.tabWidget_graphicalPanel.setTabVisible(6,True)
            elif sender.text() == 'Multiple files':
                self.selectedProcessing = 'batchProcessing'
                self.pathForFile = None
                self.lineEdit_filePath.setText("")
                self.comboBox_fileVisualized.setEnabled(True)
                self.pushButton_updatePlot.setEnabled(True)
                self.pushButton_fileBackwards.setEnabled(True)
                self.pushButton_fileForwards.setEnabled(True)
                self.tabWidget_graphicalPanel.setTabVisible(5,True)
                self.tabWidget_graphicalPanel.setTabVisible(6,True)
            else:
                #State not to be reached.
                QApplication.quit()
            print(self.selectedProcessing) if debugActivatedForDevelopment else None
    def update_desiredChannel(self, new_value):
        self.desiredChannel=int(new_value)
        print(self.desiredChannel) if debugActivatedForDevelopment else None
    def update_calibrationFactor(self, new_value):
        try:
            self.calibrationFactor=1.0/float(new_value)
        except Exception as e:
            None
        print(self.calibrationFactor) if debugActivatedForDevelopment else None   
    def update_samplingFrequencyOriginal(self, new_value):
        try:
            self.samplingFrequencyOriginal=int(new_value)
        except Exception as e:
            #Expection to deal with when the field is
            self.samplingFrequencyOriginal=None
        print(self.samplingFrequencyOriginal) if debugActivatedForDevelopment else None   
    def openFilePathDialog(self):
        #Create a QFileDialog instance
        file_dialog = QFileDialog(self)
        if self.selectedProcessing == 'singleFile':
            # Show the file dialog and get the selected file
            file_types = "".join((self.selectedFileType[0]," ","(",self.selectedFileExtension,");;All Files (*)"))
            selected_file_path, selected_filter = file_dialog.getOpenFileName(self, "Open File", "", file_types)
            #Verify the user has selected an appropriate file extension
            if selected_file_path and selected_file_path.endswith(self.selectedFileExtension[-3:]):
                self.pathForFile = selected_file_path
                self.lineEdit_filePath.setText(selected_file_path)
                self.comboBox_fileVisualized.clear()
                self.comboBox_fileVisualized.addItems([selected_file_path.split('/')[-1]])
                print(self.pathForFile) if debugActivatedForDevelopment else None   
            else:
                message = ">{0}: Error opening file: please select a .{1} file.".format(datetime.datetime.now(),self.selectedFileExtension[-3:])
                self.textBrowser_log.append(message)
                self.pathForFile = None
        elif self.selectedProcessing == 'batchProcessing':
            # Show the folder dialog and get the selected folder path
            selected_folder_path = file_dialog.getExistingDirectory(self, "Select Folder", "")
            # Process the selected folder (e.g., print the folder path)
            if selected_folder_path:
                self.pathForFile = selected_folder_path
                self.lineEdit_filePath.setText(selected_folder_path)
                print(self.pathForFile) if debugActivatedForDevelopment else None 
                #Retrieve list of files to read
                self.filesToRead = [x for x in os.listdir(self.pathForFile) if x[-4:] == self.selectedFileExtension[-4:]] #Get array with file names in the selected folder
                self.comboBox_fileVisualized.clear()
                self.comboBox_fileVisualized.addItems(self.filesToRead)
                self.comboBox_fileVisualized.setCurrentIndex(0)
            else:
                message = ">{0}: Please select a folder containing .{1} files.".format(datetime.datetime.now(),self.selectedFileExtension[-3:])
                self.textBrowser_log.append(message)
                self.pathForFile = None
    def updateCurrentAnalysisConfiguration(self):
        self.currentAnalysisConfiguration['pathForFile']=self.pathForFile
        self.currentAnalysisConfiguration['fileVisualizedPath']=self.fileVisualizedPath
        self.currentAnalysisConfiguration['filesToRead']=self.filesToRead
        self.currentAnalysisConfiguration['selectedFileType']=self.selectedFileType
        self.currentAnalysisConfiguration['selectedFileExtension']=self.selectedFileExtension        
        self.currentAnalysisConfiguration['selectedSystem']=self.selectedSystem
        self.currentAnalysisConfiguration['selectedProcessing']=self.selectedProcessing
        self.currentAnalysisConfiguration['desiredChannel']=self.desiredChannel
        self.currentAnalysisConfiguration['calibrationFactor']=self.calibrationFactor
        self.currentAnalysisConfiguration['samplingFrequencyOriginal']=self.samplingFrequencyOriginal
        self.currentAnalysisConfiguration['filterConfiguration']=self.filterConfiguration
        self.currentAnalysisConfiguration['nps']=self.nps
        self.currentAnalysisConfiguration['modalIdentificationMethodToPerform']=self.modalIdentificationMethodToPerform
        self.currentAnalysisConfiguration['intervalForAveragingHz']=self.intervalForAveragingHz
        self.currentAnalysisConfiguration['fint_BFD']=self.fint_BFD
        self.currentAnalysisConfiguration['fint_EFDD']=self.fint_EFDD
        self.currentAnalysisConfiguration['tint']=self.tint.tolist()
        self.currentAnalysisConfiguration['numModesToBeConsidered']=self.numModesToBeConsidered
        self.currentAnalysisConfiguration['i']=self.i
        self.currentAnalysisConfiguration['startingOrderNumber']=self.startingOrderNumber
        self.currentAnalysisConfiguration['endOrderNumber']=self.endOrderNumber
        self.currentAnalysisConfiguration['incrementBetweenOrder']=self.incrementBetweenOrder
        self.currentAnalysisConfiguration['refs']=self.refs
        self.currentAnalysisConfiguration['tol']=self.tol.tolist()
        self.currentAnalysisConfiguration['plotConfiguration']=self.plotConfiguration
    def readCurrentAnalysisConfiguration(self):
        self.pathForFile=self.currentAnalysisConfiguration['pathForFile']
        self.fileVisualizedPath=self.currentAnalysisConfiguration['fileVisualizedPath']
        self.filesToRead=self.currentAnalysisConfiguration['filesToRead']
        self.selectedFileType=self.currentAnalysisConfiguration['selectedFileType']
        self.selectedFileExtension=self.currentAnalysisConfiguration['selectedFileExtension']
        self.selectedSystem=self.currentAnalysisConfiguration['selectedSystem']
        self.selectedProcessing=self.currentAnalysisConfiguration['selectedProcessing']
        self.desiredChannel=self.currentAnalysisConfiguration['desiredChannel']
        self.calibrationFactor=self.currentAnalysisConfiguration['calibrationFactor']
        self.samplingFrequencyOriginal=self.currentAnalysisConfiguration['samplingFrequencyOriginal']
        self.filterConfiguration=self.currentAnalysisConfiguration['filterConfiguration']
        self.nps=self.currentAnalysisConfiguration['nps']
        self.modalIdentificationMethodToPerform=self.currentAnalysisConfiguration['modalIdentificationMethodToPerform']
        self.intervalForAveragingHz=self.currentAnalysisConfiguration['intervalForAveragingHz']
        self.fint_BFD=self.currentAnalysisConfiguration['fint_BFD']
        self.fint_EFDD=self.currentAnalysisConfiguration['fint_EFDD']
        self.tint=np.array(self.currentAnalysisConfiguration['tint'])
        self.numModesToBeConsidered=self.currentAnalysisConfiguration['numModesToBeConsidered']
        self.i=self.currentAnalysisConfiguration['i']
        self.startingOrderNumber=self.currentAnalysisConfiguration['startingOrderNumber']
        self.endOrderNumber=self.currentAnalysisConfiguration['endOrderNumber']
        self.incrementBetweenOrder=self.currentAnalysisConfiguration['incrementBetweenOrder']
        self.refs=self.currentAnalysisConfiguration['refs']
        self.tol=np.array(self.currentAnalysisConfiguration['tol'])
        self.plotConfiguration=self.currentAnalysisConfiguration['plotConfiguration']
    def resetAnalysisConfiguration(self):
        #Reset currentAnalysisConfiguration variable
        for element in self.currentAnalysisConfiguration:
            self.currentAnalysisConfiguration[element] = None
        #TODO:  These configuration variables should be put in a list or dict, as the currentAnalysisConfiguration is. 
        #       Or they should be merged in a single variable
        #       I don't recall why it was created two distinct variables actually
        #Because self.currentAnalysisConfiguration was set to None above, this will also reset the configuration variables
        self.readCurrentAnalysisConfiguration()

    ##Frequency and Damping panel
    def openFrequencyEvolutionPathDialog(self):
        #Create a QFileDialog instance
        file_dialog = QFileDialog(self)
        selected_folder_path = file_dialog.getExistingDirectory(self, "Select Folder", "")
        # Process the selected folder (e.g., print the folder path)
        if selected_folder_path:
            self.lineEdit_filePathFrequencyEvolution.setText(selected_folder_path)
            print(selected_folder_path) if debugActivatedForDevelopment else None 
    def updateFrequencyEvolutionPath(self, new_value):
        if os.path.exists(new_value):
            self.toolButton_exportFrequencyEvolution.setEnabled(True)
        else:
            message = ">{0}: Folder does not exist. Please select an existing folder.".format(datetime.datetime.now())
            self.textBrowser_log.append(message)
            self.toolButton_exportFrequencyEvolution.setEnabled(False)
    def exportFrequencyEvolution(self):
        saveFilePrefix=self.lineEdit_filePathFrequencyEvolution.text()+'/'+self.lineEdit_fileNameFrequencyEvolution.text()
        #Flag to inform the user whether a modal analysis was set
        #If no modal analysis was set to be done, then the exporting will do nothing and the user might get confused
        exportingDone = False
        if self.modalIdentificationMethodToPerform['peak-picking']==True:
            np.savetxt(saveFilePrefix+'_Frequency_PP.txt', np.vstack((self.agesOfMeasurementOriginal, self.FPP)).T, delimiter='\t', fmt='%f', header="=======================\nAGE(DAYS)\tFREQUENCY(HZ)")
            exportingDone = True
        if self.modalIdentificationMethodToPerform['BFD']==True:
            np.savetxt(saveFilePrefix+'_Frequency_BFD.txt', np.vstack((self.agesOfMeasurementOriginal, self.FBFD)).T, delimiter='\t', fmt='%f', header="=======================\nAGE(DAYS)\tFREQUENCY(HZ)")
            exportingDone = True
        if self.modalIdentificationMethodToPerform['EFDD']==True:
            np.savetxt(saveFilePrefix+'_Frequency_EFFD.txt', np.vstack((self.agesOfMeasurementOriginal, self.FEFDD)).T, delimiter='\t', fmt='%f', header="=======================\nAGE(DAYS)\tFREQUENCY(HZ)")
            exportingDone = True
        if self.modalIdentificationMethodToPerform['SSI-COV']==True:
            np.savetxt(saveFilePrefix+'_Frequency_SSI.txt', np.vstack((self.agesOfMeasurementOriginal, self.FSSI_CLUSTER.T)).T, delimiter='\t', fmt='%f', header="=======================\nAGE(DAYS)\tFREQUENCIES(Hz)")
            np.savetxt(saveFilePrefix+'_numPoles_SSI.txt', np.vstack((self.agesOfMeasurementOriginal, self.numStablePoles_CLUSTER.T)).T, delimiter='\t', fmt='%f', header="=======================\nAGE(DAYS)\tNUMBER POLES(counts)")
            exportingDone = True
        #Log message to inform status of exporting
        if exportingDone is False:
            message = ">{0}: No modal analysis method was activated, so no file was exported.".format(datetime.datetime.now())
            self.textBrowser_log.append(message)
        else:
            message = ">{0}: File exported successfully.".format(datetime.datetime.now())
            self.textBrowser_log.append(message)
    
    def openDampingEvolutionPathDialog(self):
        #Create a QFileDialog instance
        file_dialog = QFileDialog(self)
        selected_folder_path = file_dialog.getExistingDirectory(self, "Select Folder", "")
        # Process the selected folder (e.g., print the folder path)
        if selected_folder_path:
            self.lineEdit_filePathDampingEvolution.setText(selected_folder_path)
            print(selected_folder_path) if debugActivatedForDevelopment else None 
    def updateDampingEvolutionPath(self, new_value):
        if os.path.exists(new_value):
            self.toolButton_exportDampingEvolution.setEnabled(True)
        else:
            message = ">{0}: Folder does not exist. Please select an existing folder.".format(datetime.datetime.now())
            self.textBrowser_log.append(message)
            self.toolButton_exportDampingEvolution.setEnabled(False)
    def exportDampingEvolution(self):
        saveFilePrefix=self.lineEdit_filePathDampingEvolution.text()+'/'+self.lineEdit_fileNameDampingEvolution.text()
        #Flag to inform the user whether a modal analysis was set
        #If no modal analysis was set to be done, then the exporting will do nothing and the user might get confused
        exportingDone = False
        if self.modalIdentificationMethodToPerform['peak-picking']==True:
            np.savetxt(saveFilePrefix+'_Damping_PP.txt', np.vstack((self.agesOfMeasurementOriginal, self.ZPP_HP)).T, delimiter='\t', fmt='%f', header="=======================\nAGE(DAYS)\tDAMPING(%)")
            exportingDone = True
        if self.modalIdentificationMethodToPerform['BFD']==True:
            np.savetxt(saveFilePrefix+'_Damping_BFD.txt', np.vstack((self.agesOfMeasurementOriginal, self.ZBFD_FT)).T, delimiter='\t', fmt='%f', header="=======================\nAGE(DAYS)\tDAMPING(%)")
            exportingDone = True
        if self.modalIdentificationMethodToPerform['EFDD']==True:
            np.savetxt(saveFilePrefix+'_Damping_EFFD.txt', np.vstack((self.agesOfMeasurementOriginal, self.ZEFDD)).T, delimiter='\t', fmt='%f', header="=======================\nAGE(DAYS)\tDAMPING(%)")
            exportingDone = True
        if self.modalIdentificationMethodToPerform['SSI-COV']==True:
            np.savetxt(saveFilePrefix+'_Damping_SSI.txt', np.vstack((self.agesOfMeasurementOriginal, self.ZSSI_CLUSTER.T)).T, delimiter='\t', fmt='%f', header="=======================\nAGE(DAYS)\tDAMPING(%)")
            exportingDone = True
        #Log message to inform status of exporting
        if exportingDone is False:
            message = ">{0}: No modal analysis method was activated, so no file was exported.".format(datetime.datetime.now())
            self.textBrowser_log.append(message)
        else:
            message = ">{0}: File exported successfully.".format(datetime.datetime.now())
            self.textBrowser_log.append(message)

    ##File visualized panel
    def updateFileBackwards(self):
        if (self.comboBox_fileVisualized.currentIndex()-1) < 0:
            self.comboBox_fileVisualized.setCurrentIndex(0) 
        else:
            self.comboBox_fileVisualized.setCurrentIndex(self.comboBox_fileVisualized.currentIndex()-1) 
    def updateFileForwards(self):
        if (self.comboBox_fileVisualized.currentIndex()+1) >= self.comboBox_fileVisualized.count():
            self.comboBox_fileVisualized.setCurrentIndex(self.comboBox_fileVisualized.count()-1)
        else:
            self.comboBox_fileVisualized.setCurrentIndex(self.comboBox_fileVisualized.currentIndex()+1)  
    def update_fileVisualized(self, new_value):
        self.fileVisualizedPath = self.pathForFile+"/"+self.comboBox_fileVisualized.currentText()
    def update_individualPlotInBatchAnalysis(self):
        self.runAnalysis(updatePlotInBatchAnalysis=True)

    ##Filtering panel
    #Detrend groupbox
    def update_activateDetrendFilter(self, state):
        #if state == Qt.CheckState.Checked.value:
        if self.checkBox_detrend.isChecked():
            #If the checkBox is checked, user wants to include detrending filter
            #Always insert detrend in the beginning of the list
            self.filterConfiguration.insert(0,{'filter': 'detrend','type': self.comboBox_detrendType.currentText()})
            #Activate the panel's control
            self.comboBox_detrendType.setEnabled(True)
            self.label_detrendType.setEnabled(True)
            print("Detrend is checked") if debugActivatedForDevelopment else None
        else:
            #If the checkBox is unchecked, then we need to erase the entry related to detrend from the self.filterConfiguration variable
            self.filterConfiguration = [d for d in self.filterConfiguration if not ('filter' in d and d['filter'] == 'detrend')]
            #And we need to deactivate the control
            self.comboBox_detrendType.setEnabled(False)
            self.label_detrendType.setEnabled(False)
            print("Detrend is unchecked") if debugActivatedForDevelopment else None
        print(self.filterConfiguration) if debugActivatedForDevelopment else None
    def update_detrendType(self, currentIndex): # i is an int
        #Change the 'type' key in the 'detrend' dictionary contained in self.filterConfiguration variable
        # Iterate through the list of dictionaries
        for filter in self.filterConfiguration:
            if filter.get('filter') == 'detrend':
                filter['type'] = self.comboBox_detrendType.currentText()
        print(self.filterConfiguration) if debugActivatedForDevelopment else None
    #Butterworth groupbox
    def update_activateButterworthFilter(self, state):
        if state == Qt.CheckState.Checked.value:
            #If the checkBox is checked, user wants to include a Butterworth filter
            #First we need to activate all controls in this Group
            for widget in self.groupBox_butterworth.findChildren(QWidget):
                if isinstance(widget, QCheckBox):
                    if widget.text() != 'Activate':
                        widget.setEnabled(True)
                else:
                    widget.setEnabled(True)        
            #Always insert detrend in the beginning of the list
            self.filterConfiguration.append({'filter': 'butterworth','order': self.spinBox_butterworthOrder.value(),'frequencies': [self.try_float(self.lineEdit_butterworthHighPassValue.text()),self.try_float(self.lineEdit_butterworthLowPassValue.text())]})
            print("Butterworth is checked") if debugActivatedForDevelopment else None
        else:
            #If the checkBox is unchecked, then we need to erase the entry related to detrend from the self.filterConfiguration variable
            self.filterConfiguration = [d for d in self.filterConfiguration if not ('filter' in d and d['filter'] == 'butterworth')]
            #And we need to deactivate the control
            for widget in self.groupBox_butterworth.findChildren(QWidget):
                if isinstance(widget, QCheckBox):
                    if widget.text() != 'Activate':
                        widget.setEnabled(False)
                else:
                    widget.setEnabled(False)
            print("Butterworth is unchecked") if debugActivatedForDevelopment else None
        print(self.filterConfiguration) if debugActivatedForDevelopment else None
    def update_butterworthOrder(self, new_value):
        for filter in self.filterConfiguration:
            if filter.get('filter') == 'butterworth':
                filter['order'] = int(new_value)
        print(self.filterConfiguration) if debugActivatedForDevelopment else None
    def update_butterworthLowPassFrequency(self, new_value):
        for filter in self.filterConfiguration:
            if filter.get('filter') == 'butterworth':
                try:
                    filter['frequencies'][1] = self.try_float(new_value)
                except Exception as e:
                    None
        print(self.filterConfiguration) if debugActivatedForDevelopment else None   
    def update_butterworthHighPassFrequency(self, new_value):
        for filter in self.filterConfiguration:
            if filter.get('filter') == 'butterworth':
                try:
                    filter['frequencies'][0] = self.try_float(new_value)
                except Exception as e:
                    None
        print(self.filterConfiguration) if debugActivatedForDevelopment else None
    #Decimation groupbox
    def update_activateDecimationFilter(self, state):
        if state == Qt.CheckState.Checked.value:
            #If the checkBox is checked, user wants to include decimation filter
            #Always insert detrend in the beginning of the list
            self.filterConfiguration.append({'filter': 'decimation','decimationFactor': self.spinBox_decimationFactor.value()})
            #Activate the panel's control
            self.label_decimationFactor.setEnabled(True)
            self.spinBox_decimationFactor.setEnabled(True)
            print("Decimation is checked") if debugActivatedForDevelopment else None
        else:
            #If the checkBox is unchecked, then we need to erase the entry related to detrend from the self.filterConfiguration variable
            self.filterConfiguration = [d for d in self.filterConfiguration if not ('filter' in d and d['filter'] == 'decimation')]
            #And we need to deactivate the control
            self.label_decimationFactor.setEnabled(False)
            self.spinBox_decimationFactor.setEnabled(False)
            print("Decimation is unchecked") if debugActivatedForDevelopment else None
        print(self.filterConfiguration) if debugActivatedForDevelopment else None
    def update_decimationFactor(self,new_value):
        for filter in self.filterConfiguration:
            if filter.get('filter') == 'decimation':
                filter['decimationFactor'] = int(new_value)
        print(self.filterConfiguration) if debugActivatedForDevelopment else None
    
    ##Modal analysis methods panel
    #PSD groupbox
    def update_windowLengthPSD(self,new_value):
        try:
            self.nps=int(new_value)
        except Exception as e:
            #Expection to deal with when the field is empty
            self.nps=None
        print(self.nps) if debugActivatedForDevelopment else None
    #PP groupbox
    def update_activatePP(self, state):
        if state == Qt.CheckState.Checked.value:
            #If the checkBox is checked, user wants to include a Peak-Picking method
            #First we need to activate all controls in this Group
            for widget in self.groupBox_PP.findChildren(QWidget):
                if isinstance(widget, QCheckBox):
                    if widget.text() != 'Activate':
                        widget.setEnabled(True)
                else:
                    widget.setEnabled(True)
            self.tabWidget_graphicalPanel.setTabVisible(1,True)
            #Then, we need to enable the corresponding graphic tab
            #Make peak-picking active in the self.modalIdentificationMethodToPerform dict
            self.modalIdentificationMethodToPerform['peak-picking']=True
            self.intervalForAveragingHz=self.try_float(self.lineEdit_averagingInterval_PP.text())
            print("Peak-Picking is checked") if debugActivatedForDevelopment else None
        else:
            #If the checkBox is unchecked, then we need to update the self.modalIdentificationMethodToPerform dict
            self.modalIdentificationMethodToPerform['peak-picking']=False
            self.intervalForAveragingHz=None
            #And we need to deactivate the control
            for widget in self.groupBox_PP.findChildren(QWidget):
                if isinstance(widget, QCheckBox):
                    if widget.text() != 'Activate':
                        widget.setEnabled(False)
                else:
                    widget.setEnabled(False)
            self.tabWidget_graphicalPanel.setTabVisible(1,False)
            print("Peak-Picking is unchecked") if debugActivatedForDevelopment else None
        print("Modal identification dict:",self.modalIdentificationMethodToPerform,"\nPeak-picking variables: ",self.intervalForAveragingHz) if debugActivatedForDevelopment else None
    def update_intervalForAveragingHz(self,new_value):
        try:
            self.intervalForAveragingHz=self.try_float(new_value)
        except Exception as e:
            #Expection to deal with when the field is empty
            self.intervalForAveragingHz=None
        print(self.intervalForAveragingHz) if debugActivatedForDevelopment else None
    #BFD groupbox
    def update_activateBFD(self, state):
        if state == Qt.CheckState.Checked.value:
            #If the checkBox is checked, user wants to include a BFD method
            #First we need to activate all controls in this Group
            for widget in self.groupBox_BFD.findChildren(QWidget):
                if isinstance(widget, QCheckBox):
                    if widget.text() != 'Activate':
                        widget.setEnabled(True)
                else:
                    widget.setEnabled(True)        
            self.tabWidget_graphicalPanel.setTabVisible(2,True)
            #Make BFD active in the self.modalIdentificationMethodToPerform dict
            self.modalIdentificationMethodToPerform['BFD']=True
            self.fint_BFD=[self.try_float(self.lineEdit_curveFittingFrequencyLowerBound_BFD.text())/100.0,self.try_float(self.lineEdit_curveFittingFrequencyUpperBound_BFD.text())/100.0]
            print("BFD is checked") if debugActivatedForDevelopment else None
        else:
            #If the checkBox is unchecked, then we need to update the self.modalIdentificationMethodToPerform dict
            self.modalIdentificationMethodToPerform['BFD']=False
            self.fint_BFD=[None,None]
            #And we need to deactivate the control
            for widget in self.groupBox_BFD.findChildren(QWidget):
                if isinstance(widget, QCheckBox):
                    if widget.text() != 'Activate':
                        widget.setEnabled(False)
                else:
                    widget.setEnabled(False)
            self.tabWidget_graphicalPanel.setTabVisible(2,False)
            print("BFD is unchecked") if debugActivatedForDevelopment else None
        print("Modal identification dict:",self.modalIdentificationMethodToPerform,"\nBFD variables: ",self.fint_BFD) if debugActivatedForDevelopment else None
    def update_lower_fint_BFD(self,new_value):
        try:
            self.fint_BFD[0]=self.try_float(new_value)/100.0
        except Exception as e:
            #Expection to deal with when the field is empty
            self.fint_BFD[0]=None
        print("fint_BFD:",self.fint_BFD) if debugActivatedForDevelopment else None
    def update_upper_fint_BFD(self,new_value):
        try:
            self.fint_BFD[1]=self.try_float(new_value)/100.0
        except Exception as e:
            #Expection to deal with when the field is empty
            self.fint_BFD[1]=None
        print("fint_BFD:",self.fint_BFD) if debugActivatedForDevelopment else None
    #EFDD groupbox
    def update_activateEFDD(self, state):
        if state == Qt.CheckState.Checked.value:
            #If the checkBox is checked, user wants to include a EFDD method
            #First we need to activate all controls in this Group
            for widget in self.groupBox_EFDD.findChildren(QWidget):
                if isinstance(widget, QCheckBox):
                    if widget.text() != 'Activate':
                        widget.setEnabled(True)
                else:
                    widget.setEnabled(True)        
            self.tabWidget_graphicalPanel.setTabVisible(3,True)
            #Make EFDD active in the self.modalIdentificationMethodToPerform dict
            self.modalIdentificationMethodToPerform['EFDD']=True
            #Fetch parameters
            self.fint_EFDD=[self.try_float(self.lineEdit_curveFittingFrequencyLowerBound_EFDD.text())/100.0,self.try_float(self.lineEdit_curveFittingFrequencyUpperBound_EFDD.text())/100.0]
            self.tint= np.array([self.try_float(self.lineEdit_curveFittingTimeLowerBound_EFDD.text()),self.try_float(self.lineEdit_curveFittingTimeUpperBound_EFDD.text())])
            print("EFDD is checked") if debugActivatedForDevelopment else None
        else:
            #If the checkBox is unchecked, then we need to update the self.modalIdentificationMethodToPerform dict
            self.modalIdentificationMethodToPerform['EFDD']=False
            self.fint_EFDD=[None,None]
            self.tint= np.array([None,None])
            #And we need to deactivate the control
            for widget in self.groupBox_EFDD.findChildren(QWidget):
                if isinstance(widget, QCheckBox):
                    if widget.text() != 'Activate':
                        widget.setEnabled(False)
                else:
                    widget.setEnabled(False)
            self.tabWidget_graphicalPanel.setTabVisible(3,False)
            print("EFDD is unchecked") if debugActivatedForDevelopment else None
        print("Modal identification dict:",self.modalIdentificationMethodToPerform,"\nEFDD variables: fint",self.fint_EFDD, ", tint", self.tint) if debugActivatedForDevelopment else None
    def update_lower_fint_EFDD(self,new_value):
        try:
            self.fint_EFDD[0]=self.try_float(new_value)/100.0
        except Exception as e:
            #Expection to deal with when the field is empty
            self.fint_EFDD[0]=None
        print("fint_EFDD:",self.fint_EFDD) if debugActivatedForDevelopment else None
    def update_upper_fint_EFDD(self,new_value):
        try:
            self.fint_EFDD[1]=self.try_float(new_value)/100.0
        except Exception as e:
            #Expection to deal with when the field is empty
            self.fint_EFDD[1]=None
        print("fint_EFDD:",self.fint_EFDD) if debugActivatedForDevelopment else None
    def update_lower_tint(self,new_value):
        try:
            self.tint[0]=self.try_float(new_value)
        except Exception as e:
            #Expection to deal with when the field is empty
            self.tint[0]=None
        print("tint:",self.tint) if debugActivatedForDevelopment else None
    def update_upper_tint(self,new_value):
        try:
            self.tint[1]=self.try_float(new_value)
        except Exception as e:
            #Expection to deal with when the field is empty
            self.tint[1]=None
        print("tint:",self.tint) if debugActivatedForDevelopment else None
    #SSI-COV groupbox
    def update_activateSSI(self, state):
        if state == Qt.CheckState.Checked.value:
            #If the checkBox is checked, user wants to include a SSI method
            #First we need to activate all controls in this Group
            for widget in self.groupBox_SSI.findChildren(QWidget):
                if isinstance(widget, QCheckBox):
                    if widget.text() != 'Activate':
                        widget.setEnabled(True)
                else:
                    widget.setEnabled(True)        
            self.tabWidget_graphicalPanel.setTabVisible(4,True)
            #Make SSI active in the self.modalIdentificationMethodToPerform dict
            self.modalIdentificationMethodToPerform['SSI-COV']=True
            #Fetch parameters
            self.numModesToBeConsidered = self.spinBox_numModes_SSI.value()
            self.i = self.spinBox_numTimeLags_SSI.value()
            self.startingOrderNumber = self.spinBox_startingModelOrder_SSI.value()
            self.endOrderNumber = self.spinBox_endingModelOrder_SSI.value()
            self.incrementBetweenOrder = self.spinBox_orderIncrement_SSI.value()
            self.tol[0][0] = self.try_float(self.lineEdit_stabTolFrequency_allowedVariation_SSI.text())/100.0 
            self.tol[0][1] = self.try_float(self.lineEdit_stabTolFrequency_lowerBound_SSI.text()) 
            self.tol[0][2] = self.try_float(self.lineEdit_stabTolFrequency_upperBound_SSI.text()) 
            self.tol[1][0] = self.try_float(self.lineEdit_stabTolDamping_allowedVariation_SSI.text())/100.0 
            self.tol[1][1] = self.try_float(self.lineEdit_stabTolDamping_lowerBound_SSI.text())/100.0 
            self.tol[1][2] = self.try_float(self.lineEdit_stabTolDamping_upperBound_SSI.text())/100.0 
            print("SSI-COV is checked") if debugActivatedForDevelopment else None
        else:
            #If the checkBox is unchecked, then we need to update the self.modalIdentificationMethodToPerform dict
            self.modalIdentificationMethodToPerform['SSI-COV']=False
            self.numModesToBeConsidered = None
            self.i = None
            self.startingOrderNumber = None
            self.endOrderNumber = None
            self.incrementBetweenOrder = None
            self.refs = [0]
            self.tol =  np.array(
                        ([None,None,None],     # [allowed_variation, lower_bound, higher_bound]: limits for eigenfrequencies
                        [None,None,None],     # [allowed_variation, lower_bound, higher_bound]: limits for damping ratios
                        [None,None,None]))    # [allowed_variation, lower_bound, higher_bound]: limits for MAC
            #And we need to deactivate the control
            for widget in self.groupBox_SSI.findChildren(QWidget):
                if isinstance(widget, QCheckBox):
                    if widget.text() != 'Activate':
                        widget.setEnabled(False)
                else:
                    widget.setEnabled(False)
            self.tabWidget_graphicalPanel.setTabVisible(4,False)
            print("SSI-COV is unchecked") if debugActivatedForDevelopment else None
        print("Modal identification dict:",self.modalIdentificationMethodToPerform,
              "\nSSI-COV variables:",
              "\nnumModesToBeConsidered",self.numModesToBeConsidered, 
              "\nnumber of time lags (i)", self.i,
              "\nstartingOrderNumber:", self.startingOrderNumber,
              "\nendOrderNumber:",self.endOrderNumber,
              "\nincrementBetweenOrder:",self.incrementBetweenOrder,
              "\ntol:",self.tol) if debugActivatedForDevelopment else None
    def update_numModesToBeConsidered(self, new_value):
        self.numModesToBeConsidered= int(new_value)
        print("numModesToBeConsidered:",self.numModesToBeConsidered) if debugActivatedForDevelopment else None
    def update_i(self, new_value):
        self.i= int(new_value)
        print("i:",self.i) if debugActivatedForDevelopment else None
    def update_startingOrderNumber(self, new_value):
        self.startingOrderNumber= int(new_value)
        print("startingOrderNumber:",self.startingOrderNumber) if debugActivatedForDevelopment else None
    def update_endOrderNumber(self, new_value):
        self.endOrderNumber= int(new_value)
        print("update_endOrderNumber:",self.endOrderNumber) if debugActivatedForDevelopment else None
    def update_incrementBetweenOrder(self, new_value):
        self.incrementBetweenOrder= int(new_value)
        print("update_incrementBetweenOrder:",self.incrementBetweenOrder) if debugActivatedForDevelopment else None  
    def update_tol_frequencyVariation(self,new_value):
        try:
            self.tol[0][0]=self.try_float(new_value)/100.0
        except Exception as e:
            #Expection to deal with when the field is empty
            self.tol[0][0]=None
        print("tol:",self.tol) if debugActivatedForDevelopment else None
    def update_tol_frequencyLowerBound(self,new_value):
        try:
            self.tol[0][1]=self.try_float(new_value)
        except Exception as e:
            #Expection to deal with when the field is empty
            self.tol[0][1]=None
        print("tol:",self.tol) if debugActivatedForDevelopment else None
    def update_tol_frequencyUpperBound(self,new_value):
        try:
            self.tol[0][2]=self.try_float(new_value)
        except Exception as e:
            #Expection to deal with when the field is empty
            self.tol[0][2]=None
        print("tol:",self.tol) if debugActivatedForDevelopment else None
    def update_tol_dampingVariation(self,new_value):
        try:
            self.tol[1][0]=self.try_float(new_value)/100.0
        except Exception as e:
            #Expection to deal with when the field is empty
            self.tol[1][0]=None
        print("tol:",self.tol) if debugActivatedForDevelopment else None
    def update_tol_dampingLowerBound(self,new_value):
        try:
            self.tol[1][1]=self.try_float(new_value)/100.0
        except Exception as e:
            #Expection to deal with when the field is empty
            self.tol[1][1]=None
        print("tol:",self.tol) if debugActivatedForDevelopment else None
    def update_tol_dampingUpperBound(self,new_value):
        try:
            self.tol[1][2]=self.try_float(new_value)/100.0
        except Exception as e:
            #Expection to deal with when the field is empty
            self.tol[1][2]=None
        print("tol:",self.tol) if debugActivatedForDevelopment else None
    
    ##Run analysis
    def runAnalysis(self, updatePlotInBatchAnalysis=False):
        worker = Worker(updatePlotInBatchAnalysis,self)
        worker.signals.updateGraph.connect(self.updateGraphs)
        worker.signals.updateProgressBar.connect(self.updateProgressBar)
        worker.signals.updateLogMessage.connect(self.updateLogMessage)
        # Execute
        self.threadpool.start(worker)

    def runAnalysis_legacy(self, updatePlotInBatchAnalysis=False):
        '''
        This function is here just for legacy. It was used when multithreading was not deployed yet.
        updatePlotInBatchAnalysis: bool
            This variable is used in the context of the "Update plot" button in Batch Analysis
        '''
        try:
            #Perform the competent analysis
            if self.selectedProcessing == 'singleFile' or (updatePlotInBatchAnalysis is True): 
                #Initial processing common to any modal analysis method selected
                #Check if this call is in a single file analysis, or the update of graphs in batch analysis
                if updatePlotInBatchAnalysis is False:
                    accelerationDigital = auxEMMARM.readSingleFile(self.pathForFile, self.selectedSystem,self. desiredChannel)
                    acceleration = auxEMMARM.convertToG(accelerationDigital,self.calibrationFactor, self.selectedSystem)
                    if self.selectedSystem == 'uEMMARM':
                        folderPath = "/".join(self.pathForFile.split('/')[:-1])+"/"
                        files = self.pathForFile.split('/')[-1]
                        self.samplingFrequencyOriginal = auxEMMARM.getSamplingFrequency_uEMMARM(folderPath, files, len(acceleration))
                    else:
                        #Sampling frequency has been informed by the user
                        pass
                else:
                    accelerationDigital = auxEMMARM.readSingleFile(self.fileVisualizedPath, self.selectedSystem,self. desiredChannel)
                    acceleration = auxEMMARM.convertToG(accelerationDigital,self.calibrationFactor, self.selectedSystem)
                    if self.selectedSystem == 'uEMMARM':
                        folderPath = "/".join(self.fileVisualizedPath.split('/')[:-1])+"/"
                        files = self.fileVisualizedPath.split('/')[-1]
                        self.samplingFrequencyOriginal = auxEMMARM.getSamplingFrequency_uEMMARM(folderPath, files, len(acceleration))
                    else:
                        #Sampling frequency has been informed by the user
                        pass
                accelerationFiltered, samplingFrequencyFiltered  = auxEMMARM.filtering(acceleration, self.samplingFrequencyOriginal, self.filterConfiguration)
                yk = MRPy(accelerationFiltered,fs=samplingFrequencyFiltered)
                self.figurePSD, PSD = SSI.SDM(yk, nperseg=self.nps, plot=self.plotConfiguration, window='hann', nfft=2*self.nps)
                self.progressBar.setFormat('%p%')
                self.progressBar.setValue(15)

                #Peak picking method
                if self.modalIdentificationMethodToPerform['peak-picking']==True:
                    self.figurePP, self.textualResult_PP, FPP, ZPP_HP, PSDAveragedFrequency, PSDAveragedPeakIndex, yMaxPeakIndex=auxEMMARM.averagedPeakPickingMethod(PSD, self.intervalForAveragingHz, plot=self.plotConfiguration, verbose=False, textualResults=True)
                else:
                    figIgnore, FPP, ZPP_HP, ignoreThis, PSDAveragedPeakIndex, ignoreThis=auxEMMARM.averagedPeakPickingMethod(PSD, self.intervalForAveragingHz, plot=self.plotConfiguration, verbose=False)
                self.progressBar.setValue(30)

                #Build PSD object manually
                PSD.pki  = np.array([PSDAveragedPeakIndex], dtype=int)
                PSD.MGi  = np.array([0], dtype=int)
                PSD.svi  = np.array([0], dtype=int)
                PSD.tint = self.tint
                self.figureANPSD, PSD = SSI.ANPSD_from_SDM(PSD,plot=self.plotConfiguration, mode="batch")
                self.progressBar.setValue(45)

                #BFD method
                if self.modalIdentificationMethodToPerform['BFD']==True:
                    PSD.fint = np.array([self.fint_BFD[0]*FPP, self.fint_BFD[1]*FPP])
                    self.figureBFD, self.textualResult_BFD, FBFD, ZBFD_FT, VBFD, PSD_BFD = SSI.BFD(yk, PSD, plot=self.plotConfiguration, mode='batch', verbose=False, textualResults=True)
                self.progressBar.setValue(60)

                #EFDD method
                if self.modalIdentificationMethodToPerform['EFDD']==True:
                    PSD.fint = np.array([self.fint_EFDD[0]*FPP, self.fint_EFDD[1]*FPP])
                    self.figureEFDD, self.figureAutc, self.textualResult_EFDD, FEFDD, ZEFDD, VEFDD, PSD_EFDD = SSI.EFDD(yk, PSD, plot=self.plotConfiguration, mode='batch', verbose=False, textualResults=True)
                self.progressBar.setValue(75)

                #SSI method
                if self.modalIdentificationMethodToPerform['SSI-COV']==True:
                    yk = SSI.rearrange_data(yk,self.refs) 

                    FSSI_MODEL, ZSSI_MODEL, VSSI_MODEL = SSI.SSI_COV_iterator(yk,self.i,self.startingOrderNumber, self.endOrderNumber, self.incrementBetweenOrder, plot=False)
                    
                    self.figureSSI, stableModes = SSI.stabilization_diagram(FSSI_MODEL,ZSSI_MODEL,VSSI_MODEL, tol=self.tol, plot=self.plotConfiguration, PSD=PSD, verbose=False)

                    self.textualResult_SSI, FSSI, ZSSI, VSSI, numStablePoles = SSI.stable_modes(FSSI_MODEL, ZSSI_MODEL, VSSI_MODEL, stableModes, tol=0.01, spo=10, verbose=False, textualResults=True, numStablePoles=self.numModesToBeConsidered)
                    FSSI_CLUSTER=np.zeros(self.numModesToBeConsidered)
                    ZSSI_CLUSTER=np.zeros(self.numModesToBeConsidered)
                    numStablePoles_CLUSTER=np.zeros(self.numModesToBeConsidered)
                    eigenfrequenciesIndices = np.flip(np.argsort(numStablePoles))
                    
                    if FSSI.size == 0:
                        print('No sufficiently large frequency clusters could be identified') if debugActivatedForDevelopment else None
                    else:
                        for r, l in enumerate(np.take_along_axis(FSSI, eigenfrequenciesIndices, 0)[0:self.numModesToBeConsidered]): FSSI_CLUSTER[r]=l

                    if ZSSI.size == 0:
                        print('No sufficiently large damping clusters could be identified') if debugActivatedForDevelopment else None
                    else:
                        for r, l in enumerate(np.take_along_axis(ZSSI, eigenfrequenciesIndices, 0)[0:self.numModesToBeConsidered]): ZSSI_CLUSTER[r]=100*l
                    
                    if numStablePoles.size == 0:
                        print('No sufficiently large clusters could be identified') if debugActivatedForDevelopment else None 
                    else:
                        for r, l in enumerate(np.take_along_axis(numStablePoles, eigenfrequenciesIndices, 0)[0:self.numModesToBeConsidered]): numStablePoles_CLUSTER[r]=int(l)
                self.progressBar.setValue(90)

                #Display time series and PSD
                self.figureTimeSeriesRaw=auxEMMARM.plotAccelerationTimeSeries([[acceleration,self.samplingFrequencyOriginal,'Original'],], plot=self.plotConfiguration)
                self.graph_timeSeriesRaw.update_figure(self.figureTimeSeriesRaw)
                self.figureTimeSeriesFiltered=auxEMMARM.plotAccelerationTimeSeries([[accelerationFiltered,samplingFrequencyFiltered,'Filtered'],], plot=self.plotConfiguration)
                self.graph_timeSeriesFiltered.update_figure(self.figureTimeSeriesFiltered)
                self.graph_PSD.update_figure(self.figurePSD)

                #Display peak-picking method
                if (self.plotConfiguration['typeForPeakPicking']==True) and (self.modalIdentificationMethodToPerform['peak-picking']==True):
                    self.graph_PP.update_figure(self.figurePP)
                    self.graph_NPSD.update_figure(self.figureANPSD)
                    self.textBrowser_results_PP.setText(self.textualResult_PP)
                
                #Display BFD method
                if (self.plotConfiguration['typeForBFD']==True) and (self.modalIdentificationMethodToPerform['BFD']==True):
                    self.graph_BFD.update_figure(self.figureBFD)
                    self.textBrowser_results_BFD.setText(self.textualResult_BFD)

                self.progressBar.setValue(95)

                #Display EFDD method
                if (self.plotConfiguration['typeForEFDD']=="Autocorrelation-SVD") and (self.modalIdentificationMethodToPerform['EFDD']==True):
                    self.graph_spectraEFDD.update_figure(self.figureEFDD)
                    self.graph_autcEFDD.update_figure(self.figureAutc)
                    self.textBrowser_results_EFDD.setText(self.textualResult_EFDD)
                
                #Display SSI method
                if (self.plotConfiguration['typeForEFDD-AutocorrelationFitting'] is not False) and (self.modalIdentificationMethodToPerform['SSI-COV']==True):
                    self.graph_SSI.update_figure(self.figureSSI)
                    self.textBrowser_results_SSI.setText(self.textualResult_SSI)

                self.progressBar.setValue(100)
                if updatePlotInBatchAnalysis is True:
                    self.progressBar.setFormat('Update plot complete')
                else:
                    self.progressBar.setFormat('Analysis complete')
            elif self.selectedProcessing == 'batchProcessing':

                self.agesOfMeasurementOriginal = np.zeros(len(self.filesToRead)) 

                if self.modalIdentificationMethodToPerform['peak-picking']==True:
                    self.FPP = np.zeros(len(self.filesToRead))
                    self.ZPP_HP = np.zeros(len(self.filesToRead)) #To hold half-power bandwitdh
                if self.modalIdentificationMethodToPerform['BFD']==True:
                    self.FBFD = np.zeros(len(self.filesToRead))
                    self.ZBFD_FT = np.zeros(len(self.filesToRead)) #To hold curve fitting
                if self.modalIdentificationMethodToPerform['EFDD']==True:
                    self.FEFDD = np.zeros(len(self.filesToRead))
                    self.ZEFDD = np.zeros(len(self.filesToRead))
                if self.modalIdentificationMethodToPerform['SSI-COV']==True:
                    self.numStablePoles_CLUSTER=np.zeros((len(self.filesToRead),self.numModesToBeConsidered))
                    self.FSSI_CLUSTER=np.zeros((len(self.filesToRead),self.numModesToBeConsidered))
                    self.ZSSI_CLUSTER=np.zeros((len(self.filesToRead),self.numModesToBeConsidered))

                progressStep = np.linspace(100/len(self.filesToRead),100,len(self.filesToRead))
                self.progressBar.setFormat('%p%')
                for iteration, files in enumerate(self.filesToRead):

                    accelerationDigital = auxEMMARM.readBatchFile(self.pathForFile, files, self.selectedSystem, self.desiredChannel)
                    self.agesOfMeasurementOriginal[iteration] = auxEMMARM.getAgeAtMeasurementBatchFile(self.pathForFile, files, self.filesToRead[0], self.selectedSystem)
                    
                    acceleration = auxEMMARM.convertToG(accelerationDigital,self.calibrationFactor, self.selectedSystem)
                    if self.selectedSystem == 'uEMMARM':
                        samplingFrequency = auxEMMARM.getSamplingFrequency_uEMMARM(self.pathForFile, files, len(acceleration))
                    else:
                        samplingFrequency = self.samplingFrequencyOriginal
                    
                    accelerationFiltered, samplingFrequencyFiltered  = auxEMMARM.filtering(acceleration, samplingFrequency, self.filterConfiguration)
                    
                    yk = MRPy(accelerationFiltered,fs=samplingFrequencyFiltered)
                    ignoreThis, PSD = SSI.SDM(yk, nperseg=self.nps, plot=self.plotConfiguration, window='hann', nfft=2*self.nps)

                    #1) PEAK-PICKING METHOD
                    if self.modalIdentificationMethodToPerform['peak-picking']==True:
                        ignoreThis, self.FPP[iteration], self.ZPP_HP[iteration], PSDAveragedFrequency, PSDAveragedPeakIndex, yMaxPeakIndex=auxEMMARM.averagedPeakPickingMethod(PSD, self.intervalForAveragingHz, verbose=False)
                    else:
                        ignoreThis, self.FPP[iteration], self.ZPP_HP[iteration], ignoreThis, PSDAveragedPeakIndex, ignoreThis=auxEMMARM.averagedPeakPickingMethod(PSD, 0.1, verbose=False)

                    #2) MAKE ALL SORTS OF COMPUTATIONS TO PREPARE FOR FURTHER COMPUTATIONS
                    PSD.pki  = np.array([PSDAveragedPeakIndex], dtype=int)
                    PSD.MGi  = np.array([0], dtype=int)
                    PSD.svi  = np.array([0], dtype=int)
                    PSD.tint = self.tint
                    ignoreThis, PSD = SSI.ANPSD_from_SDM(PSD,plot=self.plotConfiguration, mode="batch")

                    #3) PERFORM BFD METHOD
                    if self.modalIdentificationMethodToPerform['BFD']==True:
                        PSD.fint = np.array([self.fint_BFD[0]*self.FPP[iteration], self.fint_BFD[1]*self.FPP[iteration]])
                        ignoreThis, self.FBFD[iteration], self.ZBFD_FT[iteration], VBFD, PSD_BFD = SSI.BFD(yk, PSD, plot=self.plotConfiguration, mode='batch', verbose=False)
                        #ZBFD[iteration,:]
                    
                    #4) EFDD METHOD
                    if self.modalIdentificationMethodToPerform['EFDD']==True:
                        PSD.fint = np.array([self.fint_EFDD[0]*self.FPP[iteration], self.fint_EFDD[1]*self.FPP[iteration]])
                        ignoreThis, ignoreThis, self.FEFDD[iteration], self.ZEFDD[iteration], VEFDD, PSD_EFDD = SSI.EFDD(yk, PSD, plot=self.plotConfiguration, mode='batch', verbose=False)

                    #5) PERFORM SSI-COV METHOD
                    if self.modalIdentificationMethodToPerform['SSI-COV']==True:
                        yk = SSI.rearrange_data(yk,self.refs) 
                        FSSI_MODEL, ZSSI_MODEL, VSSI_MODEL = SSI.SSI_COV_iterator(yk,self.i,self.startingOrderNumber, self.endOrderNumber, self.incrementBetweenOrder, plot=False)
                        ignoreThis, stableModes = SSI.stabilization_diagram(FSSI_MODEL,ZSSI_MODEL,VSSI_MODEL, tol=self.tol, PSD=PSD, verbose=False)
                        FSSI, ZSSI, VSSI, numStablePoles = SSI.stable_modes(FSSI_MODEL, ZSSI_MODEL, VSSI_MODEL, stableModes, tol=0.01, spo=10, verbose=False,numStablePoles=self.numModesToBeConsidered)
                        eigenfrequenciesIndices = np.flip(np.argsort(numStablePoles))
                        for r, l in enumerate(np.take_along_axis(FSSI, eigenfrequenciesIndices, 0)[0:self.numModesToBeConsidered]): self.FSSI_CLUSTER[iteration][r]=l
                        for r, l in enumerate(np.take_along_axis(ZSSI, eigenfrequenciesIndices, 0)[0:self.numModesToBeConsidered]): self.ZSSI_CLUSTER[iteration][r]=100*l
                        for r, l in enumerate(np.take_along_axis(numStablePoles, eigenfrequenciesIndices, 0)[0:self.numModesToBeConsidered]): self.numStablePoles_CLUSTER[iteration][r]=int(l)

                    #6) POPULATE HEAT MAP VARIABLE
                    '''
                    if heatMap_BatchAnalysis['save'] == True:
                        if iteration == 0:
                            #The heatMap variable will be a 2D numpy array to store all PSDs in a organized way to plot a heat map of the EMM-ARM test
                            #The first row will contain the ages associated to each PSD
                            #The first column will contain the frequency bins of each PSD (they will be the same as all PSDs are processed in the same way)
                            #The number of frequency bins will depend on the specified frequencies of interest.
                            #The first element will contain a np.nan, as it will not contain anything meaningful (it is just the crossing of the row with ages and column with frequency bins)

                            #First, check if it is the first iteration, populate the first column (with the frequency bins)
                            indicesOfInterest=np.where(np.logical_and(PSD.f>=heatMap_BatchAnalysis['frequenciesOfInterest'][0], PSD.f<=heatMap_BatchAnalysis['frequenciesOfInterest'][1]))
                            heatMap=np.zeros((len(indicesOfInterest[0])+1,len(filesToRead)+1))
                            heatMap[0,0]=np.nan
                            heatMap[1:,0]=PSD.f[indicesOfInterest] #Remember that first row is dedicated for age associated to the PSDs
                        #Now, populate the next column with PSD
                        #If iteration = n, then the respective colum in heatMapa is n+1, as the first column is already populated by frequency bins)
                        heatMap[1:, iteration+1]=PSD.ANPSD[indicesOfInterest]
                        #Populate the first element of the current column being populated with the associated age of the PSD
                        heatMap[0, iteration+1]=agesOfMeasurementOriginal[iteration]
                    '''
                    #7) PROGRESS BAR
                    #Update progress bar
                    self.progressBar.setValue(int(progressStep[iteration]))
                
                #8) PLOT RESULTS IN FREQUENCY EVOLUTION TAB
                figureFrequencyEvolution = None
                figureDampingEvolution = None
                figureFrequencyEvolution = self.buildFrequencyEvolutionFigure()
                figureDampingEvolution = self.buildDampingEvolutionFigure()
                self.graph_frequencyEvolution.update_figure(figureFrequencyEvolution)
                self.graph_dampingEvolution.update_figure(figureDampingEvolution)
                #Generate result files with the modal identification
                '''
                if self.modalIdentificationMethodToPerform['peak-picking']==True:
                    np.savetxt(resultFilePreffix+'_Frequency_PP.txt', np.vstack((agesOfMeasurementOriginal, FPP)).T, delimiter='\t', fmt='%f', header=headerResultFiles+"=======================\nAGE(DAYS)\tFREQUENCY(HZ)")
                    np.savetxt(resultFilePreffix+'_Damping_HP_PP.txt', np.vstack((agesOfMeasurementOriginal, ZPP_HP)).T, delimiter='\t', fmt='%f', header=headerResultFiles+"=======================\nAGE(DAYS)\tHALF-POWER DAMPING RATIO(%)")
                if self.modalIdentificationMethodToPerform['BFD']==True:
                    np.savetxt(resultFilePreffix+'_Frequency_BFD.txt', np.vstack((agesOfMeasurementOriginal, FBFD)).T, delimiter='\t', fmt='%f', header=headerResultFiles+"=======================\nAGE(DAYS)\tFREQUENCY(HZ)")
                    np.savetxt(resultFilePreffix+'_Damping_FT_BFD.txt', np.vstack((agesOfMeasurementOriginal, ZBFD_FT)).T, delimiter='\t', fmt='%f', header=headerResultFiles+"=======================\nAGE(DAYS)\tFITTING DAMPING RATIO(%)")
                if self.modalIdentificationMethodToPerform['EFDD']==True:
                    np.savetxt(resultFilePreffix+'_Frequency_EFFD.txt', np.vstack((agesOfMeasurementOriginal, FEFDD)).T, delimiter='\t', fmt='%f', header=headerResultFiles+"=======================\nAGE(DAYS)\tFREQUENCY(HZ)")
                    np.savetxt(resultFilePreffix+'_Damping_EFFD.txt', np.vstack((agesOfMeasurementOriginal, ZEFDD)).T, delimiter='\t', fmt='%f', header=headerResultFiles+"=======================\nAGE(DAYS)\tDAMPING RATIO(%)")
                if self.modalIdentificationMethodToPerform['SSI-COV']==True:
                    np.savetxt(resultFilePreffix+'_Frequency_SSI.txt', np.vstack((agesOfMeasurementOriginal, FSSI_CLUSTER.T)).T, delimiter='\t', fmt='%f', header=headerResultFiles+"=======================\nAGE(DAYS)\tFREQUENCIES(Hz)")
                    np.savetxt(resultFilePreffix+'_Damping_SSI.txt', np.vstack((agesOfMeasurementOriginal, ZSSI_CLUSTER.T)).T, delimiter='\t', fmt='%f', header=headerResultFiles+"=======================\nAGE(DAYS)\tDAMPING RATIOS(%)")
                    np.savetxt(resultFilePreffix+'_numPoles_SSI.txt', np.vstack((agesOfMeasurementOriginal, numStablePoles_CLUSTER.T)).T, delimiter='\t', fmt='%f', header=headerResultFiles+"=======================\nAGE(DAYS)\tNUMBER POLES(counts)")  
                '''
                '''
                if  heatMap_BatchAnalysis['save'] == True:
                    savez_compressed(resultFilePreffix+'_'+headerResultFiles[0:3]+'_heatMap.npz', heatMap)     
                '''
                #Set progress bar message
                self.progressBar.setFormat('Analysis complete')
                self.runAnalysis(updatePlotInBatchAnalysis=True)
        except Exception as e:
            # Print the exception error message
            logMessage = ">{0}: An exception occurred. This may be due to bad input. Exception description: \n {1}".format(datetime.datetime.now(),traceback.format_exc())
            self.textBrowser_log.append(logMessage)

    #Function to validate inputs
    def updateProgressBar (self, list):
        progress, typeOfInformation = list
        if typeOfInformation == 'percentage':
            self.progressBar.setFormat('%p%')
            self.progressBar.setValue(progress)
        else:
            self.progressBar.setFormat(typeOfInformation)
    def updateGraphs (self, plotToUpdate=None):
        #Plot to update is an int from 0...N
        if plotToUpdate == 0:
            self.graph_timeSeriesRaw.update_figure(self.figureTimeSeriesRaw)
            self.graph_timeSeriesFiltered.update_figure(self.figureTimeSeriesFiltered)
            self.graph_PSD.update_figure(self.figurePSD)
        elif plotToUpdate == 1:
            self.graph_PP.update_figure(self.figurePP)
            self.graph_NPSD.update_figure(self.figureANPSD)
            self.textBrowser_results_PP.setText(self.textualResult_PP)
        elif plotToUpdate == 2:
            self.graph_BFD.update_figure(self.figureBFD)
            self.textBrowser_results_BFD.setText(self.textualResult_BFD)
        elif plotToUpdate == 3:
            self.graph_spectraEFDD.update_figure(self.figureEFDD)
            self.graph_autcEFDD.update_figure(self.figureAutc)
            self.textBrowser_results_EFDD.setText(self.textualResult_EFDD)
        elif plotToUpdate == 4:
            self.graph_SSI.update_figure(self.figureSSI)
            self.textBrowser_results_SSI.setText(self.textualResult_SSI)
        elif plotToUpdate == 5:
            figureFrequencyEvolution = self.buildFrequencyEvolutionFigure()
            figureDampingEvolution = self.buildDampingEvolutionFigure()
            self.graph_frequencyEvolution.update_figure(figureFrequencyEvolution)
            self.graph_dampingEvolution.update_figure(figureDampingEvolution)
    def updateLogMessage (self, message):
        self.textBrowser_log.append(message)

    def validateInputs (self):
        #Validate an acceleration conversion factor was inserted
        #Validate a sampling frequency was informed
        #Validate the inputs of butterworth filter
        #Validate a set of modal analysis methods were informed
        None

    #Complementary functions
    ##Plot configuration dialog
    def openPlotConfigurationDialog(self, tab=0):
        print(tab) if debugActivatedForDevelopment else None
        dlg_plotConfiguration = plotConfigDlg(self)
        dlg_plotConfiguration.tabWidget_plotConfiguration.setCurrentIndex(tab)
        dlg_plotConfiguration.exec()

        self.plotConfiguration=dlg_plotConfiguration.plotConfiguration
        print("Plot configuration dialog exit.\n",self.plotConfiguration) if debugActivatedForDevelopment else None
    def populate_plotConfiguration (self, firstInitialization=False):
        if firstInitialization is True:
            self.plotGeneral = {'frequencyBandOfInterest': [0,250], 'fontSize': 10, 'fontName':'Times New Roman', 'figSize': (5,2), 'dpi': 150,'lowerYFactorPlotPSD': 0.8, 'upperYFactorPlotPSD': 1.2}
            self.plotPSD = {'typeForPSD':'Single_PSD'} #The only option currently supported
            self.plotNPSD = {'typeForANPSD':'only_ANPSD','figSizeANPSD': self.plotGeneral['figSize']} #The only option currently supported
            self.plotPeakPicking={'typeForPeakPicking': True,'figSizePeakPicking': (5,5)}
            self.plotBFD={'typeForBFD':True, 'figSizeBFD': (5,5)}
            self.plotEFDD={'typeForEFDD-AutocorrelationFitting':True,'typeForEFDD': 'Autocorrelation-SVD','figSizeEFDD': (5,5)} #'typeForEFDD': 'Autocorrelation-SVD' is the only option currently supported option
            self.plotSSI={'typeForStabilizationDiagram': 'StabilizationPSD', 'figSizeStabilization': (5,7)}
        #Rebuilt self.plotConfiguration
        self.plotConfiguration = {}
        self.plotConfiguration.update(self.plotGeneral)
        self.plotConfiguration.update(self.plotPSD)
        self.plotConfiguration.update(self.plotNPSD)
        self.plotConfiguration.update(self.plotPeakPicking)
        self.plotConfiguration.update(self.plotBFD)
        self.plotConfiguration.update(self.plotEFDD)
        self.plotConfiguration.update(self.plotSSI)

        print(self.plotConfiguration) if debugActivatedForDevelopment else None
    def try_float(self,v):
        #Convert a string to float and handle the case the string is empty, converting it to None
        try:
            return float(v)
        except Exception:
            return None
    def buildFrequencyEvolutionFigure(self):
        fig = Figure(figsize=self.plotConfiguration['figSize'], dpi=self.plotConfiguration['dpi'])
        ax = fig.add_subplot(111)
        ax.set_xlabel('Time since start of test (hours)',size=self.plotConfiguration['fontSize'], fontname=self.plotConfiguration['fontName'])
        ax.set_ylabel('Frequency (Hz)',size=self.plotConfiguration['fontSize'], fontname=self.plotConfiguration['fontName'])
        ax.grid(which='both', axis='both', linestyle='-', color='whitesmoke') 
        ax.xaxis.set_minor_locator(MultipleLocator(24))
        if self.modalIdentificationMethodToPerform['peak-picking'] is True:
            ax.plot(self.agesOfMeasurementOriginal/(60*60),self.FPP,label='peak picking')
        if self.modalIdentificationMethodToPerform['BFD'] is True:
            ax.plot(self.agesOfMeasurementOriginal/(60*60),self.FBFD,label='BFD')
        if self.modalIdentificationMethodToPerform['EFDD'] is True:
            ax.plot(self.agesOfMeasurementOriginal/(60*60),self.FEFDD,label='EFDD')
        if self.modalIdentificationMethodToPerform['SSI-COV'] is True:
            ax.plot(self.agesOfMeasurementOriginal/(60*60),self.FSSI_CLUSTER[:,0],label='SSI-COV')
        ax.legend(loc="lower right", fontsize=self.plotConfiguration['fontSize'])
        fig.set_tight_layout(True)  
        return fig
    def buildDampingEvolutionFigure(self):
        fig = Figure(figsize=self.plotConfiguration['figSize'], dpi=self.plotConfiguration['dpi'])
        ax = fig.add_subplot(111)
        ax.set_xlabel('Time since start of test (hours)',size=self.plotConfiguration['fontSize'], fontname=self.plotConfiguration['fontName'])
        ax.set_ylabel('Damping (%)',size=self.plotConfiguration['fontSize'], fontname=self.plotConfiguration['fontName'])
        ax.grid(which='both', axis='both', linestyle='-', color='whitesmoke') 
        ax.xaxis.set_minor_locator(MultipleLocator(24))
        if self.modalIdentificationMethodToPerform['peak-picking'] is True:
            ax.plot(self.agesOfMeasurementOriginal/(60*60),self.ZPP_HP,label='peak picking')
        if self.modalIdentificationMethodToPerform['BFD'] is True:
            ax.plot(self.agesOfMeasurementOriginal/(60*60),self.ZBFD_FT,label='BFD')
        if self.modalIdentificationMethodToPerform['EFDD'] is True:
            ax.plot(self.agesOfMeasurementOriginal/(60*60),self.ZEFDD,label='EFDD')
        if self.modalIdentificationMethodToPerform['SSI-COV'] is True:
            ax.plot(self.agesOfMeasurementOriginal/(60*60),self.ZSSI_CLUSTER[:,0],label='SSI-COV')
        ax.legend(loc="lower right", fontsize=self.plotConfiguration['fontSize'])
        fig.set_tight_layout(True)    
        return fig
    
    def updateGUIDisplayValues(self): 
        #TODO: Updating this parameters programatically created conflicts with the signals
        #as the updating through a configuration file was included later in the development
        #the signals were implemented to update whenever the values of the QObjects changed
        #This causes this to fire signals that changes the configuration itself. 
        #Perhaps a refactoring on how configuration is handled cna simplify things and avoid
        #unwanted behaviours
        
        self.comboBox_typeOfSystem.setCurrentText(self.selectedSystem)
        self.spinBox_channelSelection.setValue(self.desiredChannel)
        self.lineEdit_accelConvFactor.setText(str(1.0/self.calibrationFactor))
        self.lineEdit_samplingFrequency.setText(str(self.samplingFrequencyOriginal))
        #TODO: The way the signal of the radio buttons above are implemented now makes them erase
        #the self.pathForFile before it is set here. So we have to repopulate the variable here
        #This is bad implementation and should fixed in the future
        if self.selectedProcessing == 'singleFile':
            self.radioButton_multipleFile.setChecked(False)
            self.radioButton_singleFile.setChecked(True)
            self.pathForFile = self.currentAnalysisConfiguration['pathForFile']
            self.lineEdit_filePath.setText(self.pathForFile)
        elif self.selectedProcessing == 'batchProcessing':
            self.radioButton_singleFile.setChecked(False)
            self.radioButton_multipleFile.setChecked(True)
            self.pathForFile = self.currentAnalysisConfiguration['pathForFile']
            self.lineEdit_filePath.setText(self.pathForFile)
            #Retrieve list of files to read
            self.filesToRead = [x for x in os.listdir(self.pathForFile) if x[-4:] == self.selectedFileExtension[-4:]] #Get array with file names in the selected folder
            self.comboBox_fileVisualized.clear()
            self.comboBox_fileVisualized.addItems(self.filesToRead)
            self.comboBox_fileVisualized.setCurrentIndex(0)

        #self.filterConfiguration=self.currentAnalysisConfiguration['filterConfiguration']
        #Variables to control if filters are present in the configuration
        detrendPresent=False
        butterworthPresent=False
        decimationPresent=False
        for filter in self.filterConfiguration:
            if filter['filter']=='detrend':
                self.comboBox_detrendType.setCurrentText(filter['type'])
                self.checkBox_detrend.setChecked(True)
                detrendPresent=True
            if filter['filter']=='butterworth':
                self.spinBox_butterworthOrder.setValue(filter['order'])
                self.lineEdit_butterworthLowPassValue.setText(str(filter['frequencies'][1]))
                self.lineEdit_butterworthHighPassValue.setText(str(filter['frequencies'][0]))
                self.checkBox_butterworthActivate.setChecked(True)
                butterworthPresent=True
            if filter['filter']=='decimation':
                self.spinBox_decimationFactor.setValue(filter['decimationFactor'])
                self.checkBox_decimationActivate.setChecked(True)
                decimationPresent=True
        #Now, deactivate those filters that were not present in the configuration
        if detrendPresent is False:
            self.checkBox_detrend.setChecked(False)
        if butterworthPresent is False:
            self.checkBox_butterworthActivate.setChecked(False)
        if decimationPresent is False:
            self.checkBox_decimationActivate.setChecked(False)

        self.lineEdit_windowLength_PSD.setText(str(self.nps))
        
        #self.modalIdentificationMethodToPerform=self.currentAnalysisConfiguration['modalIdentificationMethodToPerform']
        
        #Variables to control if methods are present in the configuration
        PP_Present=False
        BFD_Present=False
        EFDD_Present=False
        SSI_Present=False
        if self.modalIdentificationMethodToPerform['peak-picking'] is True:
            self.lineEdit_averagingInterval_PP.setText(str(self.intervalForAveragingHz))
            self.checkBox_activate_PP.setChecked(True)
            PP_Present=True
        if self.modalIdentificationMethodToPerform['BFD'] is True:
            self.lineEdit_curveFittingFrequencyLowerBound_BFD.setText(str(self.fint_BFD[0]*100))
            self.lineEdit_curveFittingFrequencyUpperBound_BFD.setText(str(self.fint_BFD[1]*100))
            self.checkBox_activate_BFD.setChecked(True)
            BFD_Present=True
        if self.modalIdentificationMethodToPerform['EFDD'] is True:
            self.lineEdit_curveFittingFrequencyLowerBound_EFDD.setText(str(self.fint_EFDD[0]*100))
            self.lineEdit_curveFittingFrequencyUpperBound_EFDD.setText(str(self.fint_EFDD[1]*100))
            self.lineEdit_curveFittingTimeLowerBound_EFDD.setText(str(self.tint[0]))
            self.lineEdit_curveFittingTimeUpperBound_EFDD.setText(str(self.tint[1]))            
            self.checkBox_activate_EFDD.setChecked(True)
            EFDD_Present=True
        if self.modalIdentificationMethodToPerform['SSI-COV'] is True:
            self.spinBox_numModes_SSI.setValue(self.numModesToBeConsidered)
            self.spinBox_numTimeLags_SSI.setValue(self.i)
            self.spinBox_startingModelOrder_SSI.setValue(self.startingOrderNumber)
            self.spinBox_endingModelOrder_SSI.setValue(self.endOrderNumber)
            self.spinBox_orderIncrement_SSI.setValue(self.incrementBetweenOrder)
            self.lineEdit_stabTolFrequency_allowedVariation_SSI.setText(str(self.tol[0][0]))
            self.lineEdit_stabTolFrequency_lowerBound_SSI.setText(str(self.tol[0][1]))
            self.lineEdit_stabTolFrequency_upperBound_SSI.setText(str(self.tol[0][2]))
            self.lineEdit_stabTolDamping_allowedVariation_SSI.setText(str(self.tol[1][0]))
            self.lineEdit_stabTolDamping_lowerBound_SSI.setText(str(self.tol[1][1]))
            self.lineEdit_stabTolDamping_upperBound_SSI.setText(str(self.tol[1][2]))
            self.checkBox_activate_SSI.setChecked(True)
            SSI_Present=True
        #Deactivate filters that are not present in the configuration
        if PP_Present is False:
            self.checkBox_activate_PP.setChecked(False)
        if BFD_Present is False:
            self.checkBox_activate_BFD.setChecked(False)
        if EFDD_Present is False:
            self.checkBox_activate_EFDD.setChecked(False)
        if SSI_Present is False:
            self.checkBox_activate_SSI.setChecked(False)

        #self.plotConfiguration=self.currentAnalysisConfiguration['plotConfiguration']

class aboutDlg(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        #Load the dialog's GUI
        uic.loadUi(os.path.join(basedir,'UI/aboutDialog.ui'), self)

class plotConfigDlg(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        #Load the dialog's GUI
        uic.loadUi(os.path.join(basedir,'UI/dialog_plotConfiguration.ui'), self)
        self.populate_plotConfiguration(firstInitialization=True)
        
        #Define validator
        self.lineEdit_fontSize_general.setValidator(QIntValidator())
        self.lineEdit_figureDPI_general.setValidator(QIntValidator())

        #Define signals
        self.lineEdit_frequencyBand_general.textEdited.connect(lambda new_value: self.updatePlotConfiguration(new_value, parameter=0))
        self.lineEdit_fontSize_general.textEdited.connect(lambda new_value: self.updatePlotConfiguration(new_value, parameter=1))
        self.fontComboBox_fontName_general.currentFontChanged.connect(lambda new_value: self.updatePlotConfiguration(new_value, parameter=2))
        self.lineEdit_figureSize_general.textEdited.connect(lambda new_value: self.updatePlotConfiguration(new_value, parameter=3))
        self.lineEdit_figureDPI_general.textEdited.connect(lambda new_value: self.updatePlotConfiguration(new_value, parameter=4))
        self.lineEdit_yFactorPSD_general.textEdited.connect(lambda new_value: self.updatePlotConfiguration(new_value, parameter=5))
        self.comboBox_plotType_PSD.currentIndexChanged.connect(lambda new_value: self.updatePlotConfiguration(new_value, parameter=6))
        self.comboBox_plotType_PP.currentIndexChanged.connect(lambda new_value: self.updatePlotConfiguration(new_value, parameter=7))
        self.lineEdit_figureSize_PP.textEdited.connect(lambda new_value: self.updatePlotConfiguration(new_value, parameter=8))
        self.comboBox_plotType_BFD.currentIndexChanged.connect(lambda new_value: self.updatePlotConfiguration(new_value, parameter=9))
        self.lineEdit_figureSize_BFD.textEdited.connect(lambda new_value: self.updatePlotConfiguration(new_value, parameter=10))
        self.comboBox_plotType_EFDD.currentIndexChanged.connect(lambda new_value: self.updatePlotConfiguration(new_value, parameter=11))
        self.comboBox_autocorrelationPlot_EFDD.currentIndexChanged.connect(lambda new_value: self.updatePlotConfiguration(new_value, parameter=12))
        self.lineEdit_figureSize_EFDD.textEdited.connect(lambda new_value: self.updatePlotConfiguration(new_value, parameter=13))
        self.comboBox_plotType_SSI.currentIndexChanged.connect(lambda new_value: self.updatePlotConfiguration(new_value, parameter=14))
        self.lineEdit_figureSize_SSI.textEdited.connect(lambda new_value: self.updatePlotConfiguration(new_value, parameter=15))

    def updatePlotConfiguration(self, new_value, parameter=None):
        #Check which control was updated and perform the corresponding updates in the variable
        if parameter == 0:
            try:
                new_value=new_value.split(";")
                self.plotGeneral['frequencyBandOfInterest']=[self.try_float(new_value[0]), self.try_float(new_value[1])]
            except:
                None
        elif parameter == 1:
            try:
                self.plotGeneral['fontSize']=int(new_value)
            except:
                None
        elif parameter == 2:
            try:
                self.plotGeneral['fontName']=new_value.family()
            except:
                None
        elif parameter == 3:
            try:
                new_value=new_value.split(";")
                self.plotGeneral['figSize']=[self.try_float(new_value[0]), self.try_float(new_value[1])]
            except:
                None
        elif parameter == 4:
            try:
                self.plotGeneral['dpi']=int(new_value)
            except:
                None
        elif parameter == 5:
            try:
                new_value=new_value.split(";")
                self.plotGeneral['lowerYFactorPlotPSD']=self.try_float(new_value[0])
                self.plotGeneral['upperYFactorPlotPSD']=self.try_float(new_value[1])
            except:
                None
        elif parameter == 6:
            try:
                valueSelected = self.comboBox_plotType_PSD.currentText()
                if valueSelected == 'Single_PSD':
                    self.plotPSD['typeForPSD'] = valueSelected
                elif valueSelected == 'No plot':
                    self.plotPSD['typeForPSD'] = False
            except:
                None
        elif parameter == 7:
            try:
                valueSelected = self.comboBox_plotType_PP.currentText()
                if valueSelected == 'Standard':
                    self.plotPeakPicking['typeForPeakPicking'] = True
                    self.plotNPSD['typeForANPSD'] = 'only_ANPSD'

                elif valueSelected == 'No plot':
                    self.plotPeakPicking['typeForPeakPicking'] = False
                    self.plotNPSD['typeForANPSD'] = False
            except:
                None
        elif parameter == 8:
            try:
                new_value=new_value.split(";")
                self.plotPeakPicking['figSizePeakPicking']=[self.try_float(new_value[0]), self.try_float(new_value[1])]
                self.plotPeakPicking['figSizeANPSD']=self.plotPeakPicking['figSizePeakPicking']
            except:
                None
        elif parameter == 9:
            try:
                valueSelected = self.comboBox_plotType_BFD.currentText()
                if valueSelected == 'Standard':
                    self.plotBFD['typeForBFD'] = True
                elif valueSelected == 'No plot':
                    self.plotBFD['typeForBFD'] = False
            except:
                None
        elif parameter == 10:
            try:
                new_value=new_value.split(";")
                self.plotBFD['figSizeBFD']=[self.try_float(new_value[0]), self.try_float(new_value[1])]
            except:
                None
        elif parameter == 11:
            try:
                valueSelected = self.comboBox_plotType_EFDD.currentText()
                if valueSelected == 'Autocorrelation-SVD':
                    self.plotEFDD['typeForEFDD'] = 'Autocorrelation-SVD'
                elif valueSelected == 'No plot':
                    self.plotEFDD['typeForEFDD'] = False
            except:
                None
        elif parameter == 12:
            try:
                valueSelected = self.comboBox_autocorrelationPlot_EFDD.currentText()
                if valueSelected == 'True':
                    self.plotEFDD['typeForEFDD-AutocorrelationFitting'] = True
                elif valueSelected == 'False':
                    self.plotEFDD['typeForEFDD-AutocorrelationFitting'] = False
            except:
                None
        elif parameter == 13:
            try:
                new_value=new_value.split(";")
                self.plotEFDD['figSizeEFDD']=[self.try_float(new_value[0]), self.try_float(new_value[1])]
            except:
                None
        elif parameter == 14:
            try:
                valueSelected = self.comboBox_plotType_SSI.currentText()
                if valueSelected == 'StabilizationPSD':
                    self.plotSSI['typeForStabilizationDiagram'] = 'StabilizationPSD'
                if valueSelected == 'StabilizationOnly':
                    self.plotSSI['typeForStabilizationDiagram'] = 'StabilizationOnly'
                elif valueSelected == 'No plot':
                    self.plotSSI['typeForStabilizationDiagram'] = False
            except:
                None
        elif parameter == 15:
            try:
                new_value=new_value.split(";")
                self.plotSSI['figSizeStabilization']=[self.try_float(new_value[0]), self.try_float(new_value[1])]
            except:
                None
        
        self.populate_plotConfiguration()
        print(self.plotConfiguration) if debugActivatedForDevelopment else None

    def populate_plotConfiguration (self, firstInitialization=False):
        if firstInitialization is True:
            self.plotGeneral = {'frequencyBandOfInterest':[0,250], 'fontSize':10, 'fontName':'Segoe UI', 'figSize':[5,2], 'dpi':150,'lowerYFactorPlotPSD':0.8, 'upperYFactorPlotPSD':1.2}
            self.plotPSD = {'typeForPSD': 'Single_PSD'} #The only option currently supported
            self.plotNPSD = {'typeForANPSD':'only_ANPSD','figSizeANPSD': self.plotGeneral['figSize']} #The only option currently supported
            self.plotPeakPicking={'typeForPeakPicking':True,'figSizePeakPicking':[5,5]}
            self.plotBFD={'typeForBFD':True, 'figSizeBFD':[5,5]}
            self.plotEFDD={'typeForEFDD-AutocorrelationFitting': True,'typeForEFDD':'Autocorrelation-SVD','figSizeEFDD':[5,5]} #'typeForEFDD': 'Autocorrelation-SVD' is the only option currently supported option
            self.plotSSI={'typeForStabilizationDiagram':'StabilizationPSD', 'figSizeStabilization':[5,5]}
            self.plotConfiguration = None
        #Rebuilt self.plotConfiguration
        self.plotConfiguration = {}
        self.plotConfiguration.update(self.plotGeneral)
        self.plotConfiguration.update(self.plotPSD)
        self.plotConfiguration.update(self.plotNPSD)
        self.plotConfiguration.update(self.plotPeakPicking)
        self.plotConfiguration.update(self.plotBFD)
        self.plotConfiguration.update(self.plotEFDD)
        self.plotConfiguration.update(self.plotSSI)

        print(self.plotConfiguration) if debugActivatedForDevelopment else None

    def try_float(self,v):
        #Convert a string to float and handle the case the string is empty, converting it to None
        try:
            return float(v)
        except Exception:
            return None

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

    def __init__(self, updatePlotInBatchAnalysis, analysisObject):
        super(Worker, self).__init__()
        '''
        updatePlotInBatchAnalysis: bool
            This variable is used in the context of the "Update plot" button in Batch Analysis
        '''
        # Store constructor arguments (re-used for processing)
        self.analysisObject = analysisObject
        self.updatePlotInBatchAnalysis = updatePlotInBatchAnalysis
        self.signals = WorkerSignals()

    @pyqtSlot()
    def run(self):       
        #Perform the competent analysis
        try:
            if self.analysisObject.selectedProcessing == 'singleFile' or (self.updatePlotInBatchAnalysis is True): 
                #Initial processing common to any modal analysis method selected
                #Check if this call is in a single file analysis, or the update of graphs in batch analysis
                self.signals.updateProgressBar.emit([0,'percentage'])
                if self.updatePlotInBatchAnalysis is False:
                    accelerationDigital = auxEMMARM.readSingleFile(self.analysisObject.pathForFile, self.analysisObject.selectedSystem,self.analysisObject. desiredChannel)
                    acceleration = auxEMMARM.convertToG(accelerationDigital,self.analysisObject.calibrationFactor, self.analysisObject.selectedSystem)
                    if self.analysisObject.selectedSystem == 'uEMMARM':
                        folderPath = "/".join(self.analysisObject.pathForFile.split('/')[:-1])+"/"
                        files = self.analysisObject.pathForFile.split('/')[-1]
                        self.analysisObject.samplingFrequencyOriginal = auxEMMARM.getSamplingFrequency_uEMMARM(folderPath, files, len(acceleration))
                    else:
                        #Sampling frequency has been informed by the user
                        pass
                else:
                    accelerationDigital = auxEMMARM.readSingleFile(self.analysisObject.fileVisualizedPath, self.analysisObject.selectedSystem,self.analysisObject. desiredChannel)
                    acceleration = auxEMMARM.convertToG(accelerationDigital,self.analysisObject.calibrationFactor, self.analysisObject.selectedSystem)
                    if self.analysisObject.selectedSystem == 'uEMMARM':
                        folderPath = "/".join(self.analysisObject.fileVisualizedPath.split('/')[:-1])+"/"
                        files = self.analysisObject.fileVisualizedPath.split('/')[-1]
                        self.analysisObject.samplingFrequencyOriginal = auxEMMARM.getSamplingFrequency_uEMMARM(folderPath, files, len(acceleration))
                    else:
                        #Sampling frequency has been informed by the user
                        pass
                accelerationFiltered, samplingFrequencyFiltered  = auxEMMARM.filtering(acceleration, self.analysisObject.samplingFrequencyOriginal, self.analysisObject.filterConfiguration)
                yk = MRPy(accelerationFiltered,fs=samplingFrequencyFiltered)
                self.analysisObject.figurePSD, PSD = SSI.SDM(yk, nperseg=self.analysisObject.nps, plot=self.analysisObject.plotConfiguration, window='hann', nfft=2*self.analysisObject.nps)
                self.signals.updateProgressBar.emit([10,'percentage'])

                #Display time series and PSD
                self.analysisObject.figureTimeSeriesRaw=auxEMMARM.plotAccelerationTimeSeries([[acceleration,self.analysisObject.samplingFrequencyOriginal,'Original'],], plot=self.analysisObject.plotConfiguration)
                self.analysisObject.graph_timeSeriesRaw.update_figure(self.analysisObject.figureTimeSeriesRaw)
                self.analysisObject.figureTimeSeriesFiltered=auxEMMARM.plotAccelerationTimeSeries([[accelerationFiltered,samplingFrequencyFiltered,'Filtered'],], plot=self.analysisObject.plotConfiguration)
                self.analysisObject.graph_timeSeriesFiltered.update_figure(self.analysisObject.figureTimeSeriesFiltered)
                self.analysisObject.graph_PSD.update_figure(self.analysisObject.figurePSD)
                self.signals.updateGraph.emit(0)
                self.signals.updateProgressBar.emit([20,'percentage'])

                #Peak picking method
                if self.analysisObject.modalIdentificationMethodToPerform['peak-picking']==True:
                    self.analysisObject.figurePP, self.analysisObject.textualResult_PP, FPP, ZPP_HP, PSDAveragedFrequency, PSDAveragedPeakIndex, yMaxPeakIndex=auxEMMARM.averagedPeakPickingMethod(PSD, self.analysisObject.intervalForAveragingHz, plot=self.analysisObject.plotConfiguration, verbose=False, textualResults=True)
                    #Display peak-picking method
                    self.signals.updateProgressBar.emit([30,'percentage'])
                    if (self.analysisObject.plotConfiguration['typeForPeakPicking']==True):
                        self.signals.updateGraph.emit(1)
                        self.signals.updateProgressBar.emit([40,'percentage'])
                else:
                    figIgnore, FPP, ZPP_HP, ignoreThis, PSDAveragedPeakIndex, ignoreThis=auxEMMARM.averagedPeakPickingMethod(PSD, self.analysisObject.intervalForAveragingHz, plot=self.analysisObject.plotConfiguration, verbose=False)
                    self.signals.updateProgressBar.emit([30,'percentage'])
                    self.signals.updateProgressBar.emit([40,'percentage'])


                #Build PSD object manually
                PSD.pki  = np.array([PSDAveragedPeakIndex], dtype=int)
                PSD.MGi  = np.array([0], dtype=int)
                PSD.svi  = np.array([0], dtype=int)
                PSD.tint = self.analysisObject.tint
                self.analysisObject.figureANPSD, PSD = SSI.ANPSD_from_SDM(PSD,plot=self.analysisObject.plotConfiguration, mode="batch")
                self.signals.updateProgressBar.emit([50,'percentage'])

                #BFD method
                if self.analysisObject.modalIdentificationMethodToPerform['BFD']==True:
                    PSD.fint = np.array([self.analysisObject.fint_BFD[0]*FPP, self.analysisObject.fint_BFD[1]*FPP])
                    self.analysisObject.figureBFD, self.analysisObject.textualResult_BFD, FBFD, ZBFD_FT, VBFD, PSD_BFD = SSI.BFD(yk, PSD, plot=self.analysisObject.plotConfiguration, mode='batch', verbose=False, textualResults=True)
                    self.signals.updateProgressBar.emit([60,'percentage'])
                    #Display BFD method
                    if (self.analysisObject.plotConfiguration['typeForBFD']==True) and (self.analysisObject.modalIdentificationMethodToPerform['BFD']==True):
                        self.signals.updateGraph.emit(2)
                        self.signals.updateProgressBar.emit([70,'percentage'])

                #EFDD method
                if self.analysisObject.modalIdentificationMethodToPerform['EFDD']==True:
                    PSD.fint = np.array([self.analysisObject.fint_EFDD[0]*FPP, self.analysisObject.fint_EFDD[1]*FPP])
                    self.analysisObject.figureEFDD, self.analysisObject.figureAutc, self.analysisObject.textualResult_EFDD, FEFDD, ZEFDD, VEFDD, PSD_EFDD = SSI.EFDD(yk, PSD, plot=self.analysisObject.plotConfiguration, mode='batch', verbose=False, textualResults=True)
                    self.signals.updateProgressBar.emit([75,'percentage'])
                    #Display EFDD method
                    if (self.analysisObject.plotConfiguration['typeForEFDD']=="Autocorrelation-SVD") and (self.analysisObject.modalIdentificationMethodToPerform['EFDD']==True):
                        self.signals.updateGraph.emit(3)
                        self.signals.updateProgressBar.emit([80,'percentage'])

                #SSI method
                if self.analysisObject.modalIdentificationMethodToPerform['SSI-COV']==True:
                    yk = SSI.rearrange_data(yk,self.analysisObject.refs) 
                    FSSI_MODEL, ZSSI_MODEL, VSSI_MODEL = SSI.SSI_COV_iterator(yk,self.analysisObject.i,self.analysisObject.startingOrderNumber, self.analysisObject.endOrderNumber, self.analysisObject.incrementBetweenOrder, plot=False)
                    self.analysisObject.figureSSI, stableModes = SSI.stabilization_diagram(FSSI_MODEL,ZSSI_MODEL,VSSI_MODEL, tol=self.analysisObject.tol, plot=self.analysisObject.plotConfiguration, PSD=PSD, verbose=False)
                    self.analysisObject.textualResult_SSI, FSSI, ZSSI, VSSI, numStablePoles = SSI.stable_modes(FSSI_MODEL, ZSSI_MODEL, VSSI_MODEL, stableModes, tol=0.01, spo=10, verbose=False, textualResults=True, numStablePoles=self.analysisObject.numModesToBeConsidered)
                    FSSI_CLUSTER=np.zeros(self.analysisObject.numModesToBeConsidered)
                    ZSSI_CLUSTER=np.zeros(self.analysisObject.numModesToBeConsidered)
                    numStablePoles_CLUSTER=np.zeros(self.analysisObject.numModesToBeConsidered)
                    eigenfrequenciesIndices = np.flip(np.argsort(numStablePoles))
                    self.signals.updateProgressBar.emit([90,'percentage'])
                    
                    if FSSI.size == 0:
                        print('No sufficiently large frequency clusters could be identified') if debugActivatedForDevelopment else None 
                    else:
                        for r, l in enumerate(np.take_along_axis(FSSI, eigenfrequenciesIndices, 0)[0:self.analysisObject.numModesToBeConsidered]): FSSI_CLUSTER[r]=l

                    if ZSSI.size == 0:
                        print('No sufficiently large damping clusters could be identified') if debugActivatedForDevelopment else None 
                    else:
                        for r, l in enumerate(np.take_along_axis(ZSSI, eigenfrequenciesIndices, 0)[0:self.analysisObject.numModesToBeConsidered]): ZSSI_CLUSTER[r]=100*l
                    
                    if numStablePoles.size == 0:
                        print('No sufficiently large clusters could be identified') if debugActivatedForDevelopment else None 
                    else:
                        for r, l in enumerate(np.take_along_axis(numStablePoles, eigenfrequenciesIndices, 0)[0:self.analysisObject.numModesToBeConsidered]): numStablePoles_CLUSTER[r]=int(l)

                    #Display SSI method
                    if (self.analysisObject.plotConfiguration['typeForEFDD-AutocorrelationFitting'] is not False) and (self.    analysisObject.modalIdentificationMethodToPerform['SSI-COV']==True):
                        self.signals.updateGraph.emit(4)
                    
                self.signals.updateProgressBar.emit([100,'percentage'])

                if self.updatePlotInBatchAnalysis is True:
                    self.signals.updateProgressBar.emit([100,'Update plot complete'])
                else:
                    self.signals.updateProgressBar.emit([100,'Analysis complete'])
                    self.analysisObject.progressBar.setFormat('Analysis complete')
            elif self.analysisObject.selectedProcessing == 'batchProcessing':

                self.analysisObject.agesOfMeasurementOriginal = np.zeros(len(self.analysisObject.filesToRead)) 

                #These two will always be created because .FPP attribute is used in others modal analysis methods
                self.analysisObject.FPP = np.zeros(len(self.analysisObject.filesToRead))
                self.analysisObject.ZPP_HP = np.zeros(len(self.analysisObject.filesToRead)) #To hold half-power bandwitdh
                if self.analysisObject.modalIdentificationMethodToPerform['BFD']==True:
                    self.analysisObject.FBFD = np.zeros(len(self.analysisObject.filesToRead))
                    self.analysisObject.ZBFD_FT = np.zeros(len(self.analysisObject.filesToRead)) #To hold curve fitting
                if self.analysisObject.modalIdentificationMethodToPerform['EFDD']==True:
                    self.analysisObject.FEFDD = np.zeros(len(self.analysisObject.filesToRead))
                    self.analysisObject.ZEFDD = np.zeros(len(self.analysisObject.filesToRead))
                if self.analysisObject.modalIdentificationMethodToPerform['SSI-COV']==True:
                    self.analysisObject.numStablePoles_CLUSTER=np.zeros((len(self.analysisObject.filesToRead),self.analysisObject.numModesToBeConsidered))
                    self.analysisObject.FSSI_CLUSTER=np.zeros((len(self.analysisObject.filesToRead),self.analysisObject.numModesToBeConsidered))
                    self.analysisObject.ZSSI_CLUSTER=np.zeros((len(self.analysisObject.filesToRead),self.analysisObject.numModesToBeConsidered))

                progressStep = np.linspace(100/len(self.analysisObject.filesToRead),90,len(self.analysisObject.filesToRead))
                self.signals.updateProgressBar.emit([0,'percentage'])
                lastProgress = 0
                for iteration, files in enumerate(self.analysisObject.filesToRead):
                    currentProgressStep = int(progressStep[iteration])-lastProgress

                    accelerationDigital = auxEMMARM.readBatchFile(self.analysisObject.pathForFile, files, self.analysisObject.selectedSystem, self.analysisObject.desiredChannel)
                    self.analysisObject.agesOfMeasurementOriginal[iteration] = auxEMMARM.getAgeAtMeasurementBatchFile(self.analysisObject.pathForFile, files, self.analysisObject.filesToRead[0], self.analysisObject.selectedSystem)
                    
                    acceleration = auxEMMARM.convertToG(accelerationDigital,self.analysisObject.calibrationFactor, self.analysisObject.selectedSystem)
                    if self.analysisObject.selectedSystem == 'uEMMARM':
                        samplingFrequency = auxEMMARM.getSamplingFrequency_uEMMARM(self.analysisObject.pathForFile, files, len(acceleration))
                    else:
                        samplingFrequency = self.analysisObject.samplingFrequencyOriginal
                    
                    accelerationFiltered, samplingFrequencyFiltered  = auxEMMARM.filtering(acceleration, samplingFrequency, self.analysisObject.filterConfiguration)
                    
                    yk = MRPy(accelerationFiltered,fs=samplingFrequencyFiltered)
                    ignoreThis, PSD = SSI.SDM(yk, nperseg=self.analysisObject.nps, plot=self.analysisObject.plotConfiguration, window='hann', nfft=2*self.analysisObject.nps)
                    self.signals.updateProgressBar.emit([int(lastProgress+currentProgressStep/7),'percentage'])

                    #1) PEAK-PICKING METHOD
                    if self.analysisObject.modalIdentificationMethodToPerform['peak-picking']==True:
                        ignoreThis, self.analysisObject.FPP[iteration], self.analysisObject.ZPP_HP[iteration], PSDAveragedFrequency, PSDAveragedPeakIndex, yMaxPeakIndex=auxEMMARM.averagedPeakPickingMethod(PSD, self.analysisObject.intervalForAveragingHz, verbose=False)
                    else:
                        ignoreThis, self.analysisObject.FPP[iteration], self.analysisObject.ZPP_HP[iteration], ignoreThis, PSDAveragedPeakIndex, ignoreThis=auxEMMARM.averagedPeakPickingMethod(PSD, 0.1, verbose=False)
                    self.signals.updateProgressBar.emit([int(lastProgress+currentProgressStep/6),'percentage'])

                    #2) MAKE ALL SORTS OF COMPUTATIONS TO PREPARE FOR FURTHER COMPUTATIONS
                    PSD.pki  = np.array([PSDAveragedPeakIndex], dtype=int)
                    PSD.MGi  = np.array([0], dtype=int)
                    PSD.svi  = np.array([0], dtype=int)
                    PSD.tint = self.analysisObject.tint
                    ignoreThis, PSD = SSI.ANPSD_from_SDM(PSD,plot=self.analysisObject.plotConfiguration, mode="batch")
                    self.signals.updateProgressBar.emit([int(lastProgress+currentProgressStep/5),'percentage'])

                    #3) PERFORM BFD METHOD
                    if self.analysisObject.modalIdentificationMethodToPerform['BFD']==True:
                        PSD.fint = np.array([self.analysisObject.fint_BFD[0]*self.analysisObject.FPP[iteration], self.analysisObject.fint_BFD[1]*self.analysisObject.FPP[iteration]])
                        ignoreThis, self.analysisObject.FBFD[iteration], self.analysisObject.ZBFD_FT[iteration], VBFD, PSD_BFD = SSI.BFD(yk, PSD, plot=self.analysisObject.plotConfiguration, mode='batch', verbose=False)
                        #ZBFD[iteration,:]
                    self.signals.updateProgressBar.emit([int(lastProgress+currentProgressStep/4),'percentage'])
                    
                    #4) EFDD METHOD
                    if self.analysisObject.modalIdentificationMethodToPerform['EFDD']==True:
                        PSD.fint = np.array([self.analysisObject.fint_EFDD[0]*self.analysisObject.FPP[iteration], self.analysisObject.fint_EFDD[1]*self.analysisObject.FPP[iteration]])
                        ignoreThis, ignoreThis, self.analysisObject.FEFDD[iteration], self.analysisObject.ZEFDD[iteration], VEFDD, PSD_EFDD = SSI.EFDD(yk, PSD, plot=self.analysisObject.plotConfiguration, mode='batch', verbose=False)
                    self.signals.updateProgressBar.emit([int(lastProgress+currentProgressStep/3),'percentage'])

                    #5) PERFORM SSI-COV METHOD
                    if self.analysisObject.modalIdentificationMethodToPerform['SSI-COV']==True:
                        yk = SSI.rearrange_data(yk,self.analysisObject.refs) 
                        FSSI_MODEL, ZSSI_MODEL, VSSI_MODEL = SSI.SSI_COV_iterator(yk,self.analysisObject.i,self.analysisObject.startingOrderNumber, self.analysisObject.endOrderNumber, self.analysisObject.incrementBetweenOrder, plot=False)
                        ignoreThis, stableModes = SSI.stabilization_diagram(FSSI_MODEL,ZSSI_MODEL,VSSI_MODEL, tol=self.analysisObject.tol, PSD=PSD, verbose=False)
                        FSSI, ZSSI, VSSI, numStablePoles = SSI.stable_modes(FSSI_MODEL, ZSSI_MODEL, VSSI_MODEL, stableModes, tol=0.01, spo=10, verbose=False,numStablePoles=self.analysisObject.numModesToBeConsidered)
                        eigenfrequenciesIndices = np.flip(np.argsort(numStablePoles))
                        for r, l in enumerate(np.take_along_axis(FSSI, eigenfrequenciesIndices, 0)[0:self.analysisObject.numModesToBeConsidered]): self.analysisObject.FSSI_CLUSTER[iteration][r]=l
                        for r, l in enumerate(np.take_along_axis(ZSSI, eigenfrequenciesIndices, 0)[0:self.analysisObject.numModesToBeConsidered]): self.analysisObject.ZSSI_CLUSTER[iteration][r]=100*l
                        for r, l in enumerate(np.take_along_axis(numStablePoles, eigenfrequenciesIndices, 0)[0:self.analysisObject.numModesToBeConsidered]): self.analysisObject.numStablePoles_CLUSTER[iteration][r]=int(l)
                    self.signals.updateProgressBar.emit([int(lastProgress+currentProgressStep/2),'percentage'])

                    #6) POPULATE HEAT MAP VARIABLE
                    '''
                    if heatMap_BatchAnalysis['save'] == True:
                        if iteration == 0:
                            #The heatMap variable will be a 2D numpy array to store all PSDs in a organized way to plot a heat map of the EMM-ARM test
                            #The first row will contain the ages associated to each PSD
                            #The first column will contain the frequency bins of each PSD (they will be the same as all PSDs are processed in the same way)
                            #The number of frequency bins will depend on the specified frequencies of interest.
                            #The first element will contain a np.nan, as it will not contain anything meaningful (it is just the crossing of the row with ages and column with frequency bins)

                            #First, check if it is the first iteration, populate the first column (with the frequency bins)
                            indicesOfInterest=np.where(np.logical_and(PSD.f>=heatMap_BatchAnalysis['frequenciesOfInterest'][0], PSD.f<=heatMap_BatchAnalysis['frequenciesOfInterest'][1]))
                            heatMap=np.zeros((len(indicesOfInterest[0])+1,len(filesToRead)+1))
                            heatMap[0,0]=np.nan
                            heatMap[1:,0]=PSD.f[indicesOfInterest] #Remember that first row is dedicated for age associated to the PSDs
                        #Now, populate the next column with PSD
                        #If iteration = n, then the respective colum in heatMapa is n+1, as the first column is already populated by frequency bins)
                        heatMap[1:, iteration+1]=PSD.ANPSD[indicesOfInterest]
                        #Populate the first element of the current column being populated with the associated age of the PSD
                        heatMap[0, iteration+1]=agesOfMeasurementOriginal[iteration]
                    '''
                    #7) PROGRESS BAR
                    #Update progress bar
                    self.signals.updateProgressBar.emit([int(lastProgress+currentProgressStep),'percentage'])
                    lastProgress = lastProgress+currentProgressStep
                
                #8) PLOT RESULTS IN FREQUENCY EVOLUTION TAB
                self.signals.updateGraph.emit(5)
                '''
                self.analysisObject.graph_frequencyEvolution.update_figure(figureFrequencyEvolution)
                self.analysisObject.graph_dampingEvolution.update_figure(figureDampingEvolution)
                '''
                #Generate result files with the modal identification
                '''
                if self.analysisObject.modalIdentificationMethodToPerform['peak-picking']==True:
                    np.savetxt(resultFilePreffix+'_Frequency_PP.txt', np.vstack((agesOfMeasurementOriginal, FPP)).T, delimiter='\t', fmt='%f', header=headerResultFiles+"=======================\nAGE(DAYS)\tFREQUENCY(HZ)")
                    np.savetxt(resultFilePreffix+'_Damping_HP_PP.txt', np.vstack((agesOfMeasurementOriginal, ZPP_HP)).T, delimiter='\t', fmt='%f', header=headerResultFiles+"=======================\nAGE(DAYS)\tHALF-POWER DAMPING RATIO(%)")
                if self.analysisObject.modalIdentificationMethodToPerform['BFD']==True:
                    np.savetxt(resultFilePreffix+'_Frequency_BFD.txt', np.vstack((agesOfMeasurementOriginal, FBFD)).T, delimiter='\t', fmt='%f', header=headerResultFiles+"=======================\nAGE(DAYS)\tFREQUENCY(HZ)")
                    np.savetxt(resultFilePreffix+'_Damping_FT_BFD.txt', np.vstack((agesOfMeasurementOriginal, ZBFD_FT)).T, delimiter='\t', fmt='%f', header=headerResultFiles+"=======================\nAGE(DAYS)\tFITTING DAMPING RATIO(%)")
                if self.analysisObject.modalIdentificationMethodToPerform['EFDD']==True:
                    np.savetxt(resultFilePreffix+'_Frequency_EFFD.txt', np.vstack((agesOfMeasurementOriginal, FEFDD)).T, delimiter='\t', fmt='%f', header=headerResultFiles+"=======================\nAGE(DAYS)\tFREQUENCY(HZ)")
                    np.savetxt(resultFilePreffix+'_Damping_EFFD.txt', np.vstack((agesOfMeasurementOriginal, ZEFDD)).T, delimiter='\t', fmt='%f', header=headerResultFiles+"=======================\nAGE(DAYS)\tDAMPING RATIO(%)")
                if self.analysisObject.modalIdentificationMethodToPerform['SSI-COV']==True:
                    np.savetxt(resultFilePreffix+'_Frequency_SSI.txt', np.vstack((agesOfMeasurementOriginal, FSSI_CLUSTER.T)).T, delimiter='\t', fmt='%f', header=headerResultFiles+"=======================\nAGE(DAYS)\tFREQUENCIES(Hz)")
                    np.savetxt(resultFilePreffix+'_Damping_SSI.txt', np.vstack((agesOfMeasurementOriginal, ZSSI_CLUSTER.T)).T, delimiter='\t', fmt='%f', header=headerResultFiles+"=======================\nAGE(DAYS)\tDAMPING RATIOS(%)")
                    np.savetxt(resultFilePreffix+'_numPoles_SSI.txt', np.vstack((agesOfMeasurementOriginal, numStablePoles_CLUSTER.T)).T, delimiter='\t', fmt='%f', header=headerResultFiles+"=======================\nAGE(DAYS)\tNUMBER POLES(counts)")  
                '''
                '''
                if  heatMap_BatchAnalysis['save'] == True:
                    savez_compressed(resultFilePreffix+'_'+headerResultFiles[0:3]+'_heatMap.npz', heatMap)     
                '''
                #Set progress bar message
                self.signals.updateProgressBar.emit([100,'percentage'])
                self.signals.updateProgressBar.emit([100,'Analysis complete'])
                #Update individual file plots
                self.updatePlotInBatchAnalysis=True
                self.run()
        except Exception as e:
            # Print the exception error message
            self.signals.updateProgressBar.emit([100,'percentage'])
            self.signals.updateProgressBar.emit([0,'Error has occurred'])
            logMessage = ">{0}: An exception occurred. This may be due to bad input. Exception description: \n {1}".format(datetime.datetime.now(),traceback.format_exc())
            self.signals.updateLogMessage.emit(logMessage)

def main():

    app = QtWidgets.QApplication(sys.argv)
    main = modalAnalysis_mainWindow()
    main.show()
    sys.exit(app.exec())

if __name__ == '__main__':
    main()
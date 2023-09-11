#Attributions:
#Integration Matplotlib to QTDesigner: https://yapayzekalabs.blogspot.com/2018/11/pyqt5-gui-qt-designer-matplotlib.html
#PyQt6 and QTDesigner tutotiral: https://www.pythonguis.com/pyqt6/
#Icons: the japanese guy

from PyQt6 import QtWidgets, uic
import sys
import numpy as np
import matplotlib
matplotlib.use('QtAgg')

from PyQt6.QtCore import Qt
from PyQt6 import QtWidgets
from PyQt6.QtWidgets import QApplication, QWidget, QCheckBox, QFileDialog

from PyQt6.QtGui import QDoubleValidator, QIntValidator

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.figure import Figure


debugActivated=True

class MainWindow(QtWidgets.QMainWindow):

    def __init__(self, *args, **kwargs):
        super(MainWindow, self).__init__(*args, **kwargs)

        #Load the UI Page
        uic.loadUi(r'GUI/mainwindow.ui', self)

        #Define global properties with default values
        #These values may be used as reference for validating user inputs
        self.pathForFile=None
        self.selectedFileType = ('uEMMARM files', '*.emm')
        self.selectedFileExtension = self.selectedFileType[1][-5:]
        
        self.selectedSystem='uEMMARM'
        self.selectedProcessing='singleFile'
        self.desiredChannel = 1
        self.calibrationFactor = None
        self.samplingFrequencyOriginal = None
        self.filterConfiguration = []
        self.nps = None
        self.modalIdentificationMethodToPerform = {'peak-picking':False, 
                                            'BFD':False, 
                                            'EFDD':False, 
                                            'SSI-COV':False}
        #Peak-picking properties
        self.intervalForAveragingHz=None
        #BFD properties
        self.fint_BFD = [None,None]
        #EFDD properties
        self.fint_EFDD = [None, None]
        self.tint = np.array([None, None])
        #SSI-COV properties
        self.numModesToBeConsidered = None
        self.i = None
        self.startingOrderNumber = None
        self.endOrderNumber = None
        self.incrementBetweenOrder = None
        self.refs = [0]
        self.tol =  np.array((
                    [None,None,None],     # [allowed_variation, lower_bound, higher_bound]: limits for eigenfrequencies
                    [None,None,None],     # [allowed_variation, lower_bound, higher_bound]: limits for damping ratios
                    [None,None,None]))    # [allowed_variation, lower_bound, higher_bound]: limits for MAC
        
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
        #PP groupbox
        self.checkBox_activate_PP.stateChanged.connect(self.update_activatePP)
        self.lineEdit_averagingInterval_PP.textChanged.connect(self.update_intervalForAveragingHz)
        #BFD groupbox
        self.checkBox_activate_BFD.stateChanged.connect(self.update_activateBFD)
        self.lineEdit_curveFittingFrequencyLowerBound_BFD.textChanged.connect(self.update_lower_fint_BFD)
        self.lineEdit_curveFittingFrequencyUpperBound_BFD.textChanged.connect(self.update_upper_fint_BFD)
        #EFDD groupbox
        self.checkBox_activate_EFDD.stateChanged.connect(self.update_activateEFDD)
        self.lineEdit_curveFittingFrequencyLowerBound_EFDD.textChanged.connect(self.update_lower_fint_EFDD)
        self.lineEdit_curveFittingFrequencyUpperBound_EFDD.textChanged.connect(self.update_upper_fint_EFDD)
        self.lineEdit_curveFittingTimeLowerBound_EFDD.textChanged.connect(self.update_lower_tint)
        self.lineEdit_curveFittingTimeUpperBound_EFDD.textChanged.connect(self.update_upper_tint)
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

        ##Run analysis
        self.pushButton_runAnalysis.clicked.connect(self.runAnalysis)

    #Definition of methods to perform with the signals
    ##General info panel
    def update_selectedSystem(self, currentIndex): # i is an int
        self.selectedSystem=self.comboBox_typeOfSystem.currentText()
        #Verify which fileType is to be allowed base on the informed system
        fileTypes = [('uEMMARM files', '*.emm'), ('Text files', '*.txt'), ('CSV files', '*.csv')]
        if self.selectedSystem == "uEMMARM":
            self.selectedFileType = fileTypes[0]
            self.selectedFileExtension = self.selectedFileType[1][-5:]
        elif self.selectedSystem == "National":
            self.selectedFileType = fileTypes[1]
            self.selectedFileExtension = self.selectedFileType[1][-5:]
        elif self.selectedSystem == "old_uEMMARM":
            self.selectedFileType = fileTypes[1]
            self.selectedFileExtension = self.selectedFileType[1][-5:]
        elif self.selectedSystem == "RPi":
            self.selectedFileType = fileTypes[2]
            self.selectedFileExtension = self.selectedFileType[1][-5:]
        print(self.selectedSystem) if debugActivated else None
    def update_selectedProcessing(self):
        sender = self.sender()  # Get the sender of the signal
        if sender.isChecked():
            if sender.text() == 'Single file':
                self.selectedProcessing = 'singleFile'
            elif sender.text() == 'Multiple files':
                self.selectedProcessing = 'batchProcessing'
            else:
                #State not to be reached.
                QApplication.quit()
            print(self.selectedProcessing) if debugActivated else None
    def update_desiredChannel(self, new_value):
        self.desiredChannel=int(new_value)
        print(self.desiredChannel) if debugActivated else None
    def update_calibrationFactor(self, new_value):
        try:
            self.calibrationFactor=float(new_value)
        except Exception as e:
            None
        print(self.calibrationFactor) if debugActivated else None   
    def update_samplingFrequencyOriginal(self, new_value):
        try:
            self.samplingFrequencyOriginal=int(new_value)
        except Exception as e:
            #Expection to deal with when the field is
            self.samplingFrequencyOriginal=None
        print(self.samplingFrequencyOriginal) if debugActivated else None   
    def openFilePathDialog(self):
        #Create a QFileDialog instance
        file_dialog = QFileDialog(self)
        
        if self.selectedProcessing == 'singleFile':
            # Show the file dialog and get the selected file
            file_types = "".join((self.selectedFileType[0]," ","(",self.selectedFileExtension,");;All Files (*)"))
            selected_file_path, selected_filter = file_dialog.getOpenFileName(self, "Open File", "", file_types)
            #Verify the user has selected an appropriate file extension
            if selected_file_path and selected_file_path.endswith(self.selectedFileExtension[-3:]):
                print("Selected file:", selected_file_path)
                self.pathForFile = selected_file_path
            else:
                print("Error opening file: please select a", self.selectedFileExtension[-3:]," file.")
                self.pathForFile = None
        elif self.selectedProcessing == 'batchProcessing':
            # Show the folder dialog and get the selected folder path
            selected_folder_path = file_dialog.getExistingDirectory(self, "Select Folder", "")
            
            # Process the selected folder (e.g., print the folder path)
            if selected_folder_path:
                print("Selected folder:", selected_folder_path)
                self.pathForFile = selected_folder_path
    
    ##Filtering panel
    #Detrend groupbox
    def update_activateDetrendFilter(self, state):
        if state == Qt.CheckState.Checked.value:
            #If the checkBox is checked, user wants to include detrending filter
            #Always insert detrend in the beginning of the list
            self.filterConfiguration.insert(0,{'filter': 'detrend','type': self.comboBox_detrendType.currentText()})
            #Activate the panel's control
            self.comboBox_detrendType.setEnabled(True)
            self.label_detrendType.setEnabled(True)
            print("Detrend is checked") if debugActivated else None
        else:
            #If the checkBox is unchecked, then we need to erase the entry related to detrend from the self.filterConfiguration variable
            self.filterConfiguration = [d for d in self.filterConfiguration if not ('filter' in d and d['filter'] == 'detrend')]
            #And we need to deactivate the control
            self.comboBox_detrendType.setEnabled(False)
            self.label_detrendType.setEnabled(False)
            print("Detrend is unchecked") if debugActivated else None
        print(self.filterConfiguration) if debugActivated else None
    def update_detrendType(self, currentIndex): # i is an int
        #Change the 'type' key in the 'detrend' dictionary contained in self.filterConfiguration variable
        # Iterate through the list of dictionaries
        for filter in self.filterConfiguration:
            if filter.get('filter') == 'detrend':
                filter['type'] = self.comboBox_detrendType.currentText()
        print(self.filterConfiguration) if debugActivated else None
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
            print("Butterworth is checked") if debugActivated else None
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
            print("Butterworth is unchecked") if debugActivated else None
        print(self.filterConfiguration) if debugActivated else None
    def update_butterworthOrder(self, new_value):
        for filter in self.filterConfiguration:
            if filter.get('filter') == 'butterworth':
                filter['order'] = int(new_value)
        print(self.filterConfiguration) if debugActivated else None
    def update_butterworthLowPassFrequency(self, new_value):
        for filter in self.filterConfiguration:
            if filter.get('filter') == 'butterworth':
                try:
                    filter['frequencies'][1] = self.try_float(new_value)
                except Exception as e:
                    None
        print(self.filterConfiguration) if debugActivated else None   
    def update_butterworthHighPassFrequency(self, new_value):
        for filter in self.filterConfiguration:
            if filter.get('filter') == 'butterworth':
                try:
                    filter['frequencies'][0] = self.try_float(new_value)
                except Exception as e:
                    None
        print(self.filterConfiguration) if debugActivated else None
    #Decimation groupbox
    def update_activateDecimationFilter(self, state):
        if state == Qt.CheckState.Checked.value:
            #If the checkBox is checked, user wants to include decimation filter
            #Always insert detrend in the beginning of the list
            self.filterConfiguration.append({'filter': 'decimation','decimationFactor': self.spinBox_decimationFactor.value()})
            #Activate the panel's control
            self.label_decimationFactor.setEnabled(True)
            self.spinBox_decimationFactor.setEnabled(True)
            print("Decimation is checked") if debugActivated else None
        else:
            #If the checkBox is unchecked, then we need to erase the entry related to detrend from the self.filterConfiguration variable
            self.filterConfiguration = [d for d in self.filterConfiguration if not ('filter' in d and d['filter'] == 'decimation')]
            #And we need to deactivate the control
            self.label_decimationFactor.setEnabled(False)
            self.spinBox_decimationFactor.setEnabled(False)
            print("Decimation is unchecked") if debugActivated else None
        print(self.filterConfiguration) if debugActivated else None
    def update_decimationFactor(self,new_value):
        for filter in self.filterConfiguration:
            if filter.get('filter') == 'decimation':
                filter['decimationFactor'] = int(new_value)
        print(self.filterConfiguration) if debugActivated else None
    ##Modal analysis methods panel
    #PSD groupbox
    def update_windowLengthPSD(self,new_value):
        try:
            self.nps=int(new_value)
        except Exception as e:
            #Expection to deal with when the field is empty
            self.nps=None
        print(self.nps) if debugActivated else None
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
            #Make peak-picking active in the self.modalIdentificationMethodToPerform dict
            self.modalIdentificationMethodToPerform['peak-picking']=True
            self.intervalForAveragingHz=self.try_float(self.lineEdit_averagingInterval_PP.text())
            print("Peak-Picking is checked") if debugActivated else None
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
            print("Peak-Picking is unchecked") if debugActivated else None
        print("Modal identification dict:",self.modalIdentificationMethodToPerform,"\nPeak-picking variables: ",self.intervalForAveragingHz) if debugActivated else None
    def update_intervalForAveragingHz(self,new_value):
        try:
            self.intervalForAveragingHz=self.try_float(new_value)
        except Exception as e:
            #Expection to deal with when the field is empty
            self.intervalForAveragingHz=None
        print(self.intervalForAveragingHz) if debugActivated else None
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
            #Make BFD active in the self.modalIdentificationMethodToPerform dict
            self.modalIdentificationMethodToPerform['BFD']=True
            self.fint_BFD=[self.try_float(self.lineEdit_curveFittingFrequencyLowerBound_BFD.text()),self.try_float(self.lineEdit_curveFittingFrequencyUpperBound_BFD.text())]
            print("BFD is checked") if debugActivated else None
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
            print("BFD is unchecked") if debugActivated else None
        print("Modal identification dict:",self.modalIdentificationMethodToPerform,"\nBFD variables: ",self.fint_BFD) if debugActivated else None
    def update_lower_fint_BFD(self,new_value):
        try:
            self.fint_BFD[0]=self.try_float(new_value)
        except Exception as e:
            #Expection to deal with when the field is empty
            self.fint_BFD[0]=None
        print("fint_BFD:",self.fint_BFD) if debugActivated else None
    def update_upper_fint_BFD(self,new_value):
        try:
            self.fint_BFD[1]=self.try_float(new_value)
        except Exception as e:
            #Expection to deal with when the field is empty
            self.fint_BFD[1]=None
        print("fint_BFD:",self.fint_BFD) if debugActivated else None
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
            #Make EFDD active in the self.modalIdentificationMethodToPerform dict
            self.modalIdentificationMethodToPerform['EFDD']=True
            #Fetch parameters
            self.fint_EFDD=[self.try_float(self.lineEdit_curveFittingFrequencyLowerBound_EFDD.text()),self.try_float(self.lineEdit_curveFittingFrequencyUpperBound_EFDD.text())]
            self.tint=[self.try_float(self.lineEdit_curveFittingTimeLowerBound_EFDD.text()),self.try_float(self.lineEdit_curveFittingTimeUpperBound_EFDD.text())]
            print("EFDD is checked") if debugActivated else None
        else:
            #If the checkBox is unchecked, then we need to update the self.modalIdentificationMethodToPerform dict
            self.modalIdentificationMethodToPerform['EFDD']=False
            self.fint_EFDD=[None,None]
            self.tint=[None,None]
            #And we need to deactivate the control
            for widget in self.groupBox_EFDD.findChildren(QWidget):
                if isinstance(widget, QCheckBox):
                    if widget.text() != 'Activate':
                        widget.setEnabled(False)
                else:
                    widget.setEnabled(False)
            print("EFDD is unchecked") if debugActivated else None
        print("Modal identification dict:",self.modalIdentificationMethodToPerform,"\nEFDD variables: fint",self.fint_EFDD, ", tint", self.tint) if debugActivated else None
    def update_lower_fint_EFDD(self,new_value):
        try:
            self.fint_EFDD[0]=self.try_float(new_value)
        except Exception as e:
            #Expection to deal with when the field is empty
            self.fint_EFDD[0]=None
        print("fint_EFDD:",self.fint_EFDD) if debugActivated else None
    def update_upper_fint_EFDD(self,new_value):
        try:
            self.fint_EFDD[1]=self.try_float(new_value)
        except Exception as e:
            #Expection to deal with when the field is empty
            self.fint_EFDD[1]=None
        print("fint_EFDD:",self.fint_EFDD) if debugActivated else None
    def update_lower_tint(self,new_value):
        try:
            self.tint[0]=self.try_float(new_value)
        except Exception as e:
            #Expection to deal with when the field is empty
            self.tint[0]=None
        print("tint:",self.tint) if debugActivated else None
    def update_upper_tint(self,new_value):
        try:
            self.tint[1]=self.try_float(new_value)
        except Exception as e:
            #Expection to deal with when the field is empty
            self.tint[1]=None
        print("tint:",self.tint) if debugActivated else None
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
            #Make SSI active in the self.modalIdentificationMethodToPerform dict
            self.modalIdentificationMethodToPerform['SSI-COV']=True
            #Fetch parameters
            self.numModesToBeConsidered = self.spinBox_numModes_SSI.value()
            self.i = self.spinBox_numTimeLags_SSI.value()
            self.startingOrderNumber = self.spinBox_startingModelOrder_SSI.value()
            self.endOrderNumber = self.spinBox_endingModelOrder_SSI.value()
            self.incrementBetweenOrder = self.spinBox_orderIncrement_SSI.value()
            self.tol[0][0] = self.try_float(self.lineEdit_stabTolFrequency_allowedVariation_SSI.text()) 
            self.tol[0][1] = self.try_float(self.lineEdit_stabTolFrequency_lowerBound_SSI.text()) 
            self.tol[0][2] = self.try_float(self.lineEdit_stabTolFrequency_upperBound_SSI.text()) 
            self.tol[1][0] = self.try_float(self.lineEdit_stabTolDamping_allowedVariation_SSI.text()) 
            self.tol[1][1] = self.try_float(self.lineEdit_stabTolDamping_lowerBound_SSI.text()) 
            self.tol[1][2] = self.try_float(self.lineEdit_stabTolDamping_upperBound_SSI.text()) 
            print("SSI-COV is checked") if debugActivated else None
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
            print("SSI-COV is unchecked") if debugActivated else None
        print("Modal identification dict:",self.modalIdentificationMethodToPerform,
              "\nSSI-COV variables:",
              "\nnumModesToBeConsidered",self.numModesToBeConsidered, 
              "\nnumber of time lags (i)", self.i,
              "\nstartingOrderNumber:", self.startingOrderNumber,
              "\nendOrderNumber:",self.endOrderNumber,
              "\nincrementBetweenOrder:",self.incrementBetweenOrder,
              "\ntol:",self.tol) if debugActivated else None
    def update_numModesToBeConsidered(self, new_value):
        self.numModesToBeConsidered= int(new_value)
        print("numModesToBeConsidered:",self.numModesToBeConsidered) if debugActivated else None
    def update_i(self, new_value):
        self.i= int(new_value)
        print("i:",self.i) if debugActivated else None
    def update_startingOrderNumber(self, new_value):
        self.startingOrderNumber= int(new_value)
        print("startingOrderNumber:",self.startingOrderNumber) if debugActivated else None
    def update_endOrderNumber(self, new_value):
        self.endOrderNumber= int(new_value)
        print("update_endOrderNumber:",self.endOrderNumber) if debugActivated else None
    def update_incrementBetweenOrder(self, new_value):
        self.incrementBetweenOrder= int(new_value)
        print("update_incrementBetweenOrder:",self.incrementBetweenOrder) if debugActivated else None  
    def update_tol_frequencyVariation(self,new_value):
        try:
            self.tol[0][0]=self.try_float(new_value)/100.0
        except Exception as e:
            #Expection to deal with when the field is empty
            self.tol[0][0]=None
        print("tol:",self.tol) if debugActivated else None
    def update_tol_frequencyLowerBound(self,new_value):
        try:
            self.tol[0][1]=self.try_float(new_value)
        except Exception as e:
            #Expection to deal with when the field is empty
            self.tol[0][1]=None
        print("tol:",self.tol) if debugActivated else None
    def update_tol_frequencyUpperBound(self,new_value):
        try:
            self.tol[0][2]=self.try_float(new_value)
        except Exception as e:
            #Expection to deal with when the field is empty
            self.tol[0][2]=None
        print("tol:",self.tol) if debugActivated else None
    def update_tol_dampingVariation(self,new_value):
        try:
            self.tol[1][0]=self.try_float(new_value)/100.0
        except Exception as e:
            #Expection to deal with when the field is empty
            self.tol[1][0]=None
        print("tol:",self.tol) if debugActivated else None
    def update_tol_dampingLowerBound(self,new_value):
        try:
            self.tol[1][1]=self.try_float(new_value)/100.0
        except Exception as e:
            #Expection to deal with when the field is empty
            self.tol[1][1]=None
        print("tol:",self.tol) if debugActivated else None
    def update_tol_dampingUpperBound(self,new_value):
        try:
            self.tol[1][2]=self.try_float(new_value)/100.0
        except Exception as e:
            #Expection to deal with when the field is empty
            self.tol[1][2]=None
        print("tol:",self.tol) if debugActivated else None
    
    ##Run analysis
    def runAnalysis(self):
       # Create a Matplotlib Figure and Axes
        import matplotlib.pyplot as plt
        import random
        from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas

        fig = plt.figure()
        ax = fig.add_subplot(111)

        fs = 500
        f = random.randint(1, 100)
        ts = 1/fs
        length_of_signal = 100
        t = np.linspace(0,1,length_of_signal)
        
        cosinus_signal = np.cos(2*np.pi*f*t)
        # Customize the plot
        ax.plot(t, cosinus_signal, label="Sample Data")
        ax.set_xlabel("X-axis")
        ax.set_ylabel("Y-axis")
        self.graph_timeSeries.update_figure(fig)
        

    #Function to validate inputs
    def validateInputs (self):
        #Validate an acceleration conversion factor was inserted
        #Validate a sampling frequency was informed
        #Validate the inputs of butterworth filter
        None

    #Complementary functions
    def try_float(self,v):
        #Convert a string to float and handle the case the string is empty, converting it to None
        try:
            return float(v)
        except Exception:
            return None
def main():
    app = QtWidgets.QApplication(sys.argv)
    main = MainWindow()
    main.show()
    sys.exit(app.exec())

if __name__ == '__main__':
    main()
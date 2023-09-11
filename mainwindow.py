from PyQt6 import QtWidgets, uic
import sys
import matplotlib
matplotlib.use('QtAgg')

from PyQt6.QtCore import Qt
from PyQt6 import QtWidgets
from PyQt6.QtWidgets import QApplication, QWidget, QCheckBox

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
        #EFDD properties
        #SSI-COV properties
        
        #Define validators
        self.lineEdit_accelConvFactor.setValidator(QDoubleValidator())
        self.lineEdit_samplingFrequency.setValidator(QIntValidator())
        self.lineEdit_butterworthLowPassValue.setValidator(QDoubleValidator())
        self.lineEdit_butterworthHighPassValue.setValidator(QDoubleValidator())
        self.lineEdit_windowLength_PSD.setValidator(QIntValidator())
        self.lineEdit_averagingInterval_PP.setValidator(QDoubleValidator())
        
        #Define signals
        # Sends the current index (position) of the selected item.
        ## Filtering panel
        self.comboBox_typeOfSystem.currentIndexChanged.connect(self.update_selectedSystem)
        self.radioButton_singleFile.toggled.connect(self.update_selectedProcessing)
        self.radioButton_multipleFile.toggled.connect(self.update_selectedProcessing)
        self.spinBox_channelSelection.valueChanged.connect(self.update_desiredChannel)
        self.lineEdit_accelConvFactor.textChanged.connect(self.update_calibrationFactor)
        self.lineEdit_samplingFrequency.textChanged.connect(self.update_samplingFrequencyOriginal)
        self.checkBox_detrend.stateChanged.connect(self.update_activateDetrendFilter)
        self.comboBox_detrendType.currentIndexChanged.connect(self.update_detrendType)
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

    #Definition of methods to perform with the signals
    def update_selectedSystem(self, currentIndex): # i is an int
        self.selectedSystem=self.comboBox_typeOfSystem.currentText()
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
            None
        print(self.samplingFrequencyOriginal) if debugActivated else None
    
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
            None
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
            None
        print(self.intervalForAveragingHz) if debugActivated else None


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
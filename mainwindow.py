from PyQt6 import QtWidgets, uic
import sys
import matplotlib
matplotlib.use('QtAgg')

from PyQt6 import QtCore, QtWidgets
from PyQt6.QtWidgets import QApplication

from PyQt6.QtGui import QDoubleValidator, QIntValidator

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.figure import Figure

debugActivated=True

class MainWindow(QtWidgets.QMainWindow):

    def __init__(self, *args, **kwargs):
        super(MainWindow, self).__init__(*args, **kwargs)

        #Load the UI Page
        uic.loadUi('mainwindow.ui', self)

        #Define global variables with default values
        self.selectedSystem='uEMMARM'
        self.selectedProcessing='singleFile'
        self.desiredChannel = 1
        self.calibrationFactor = None
        self.samplingFrequencyOriginal = None
        
        #Define validators
        self.lineEdit_accelConvFactor.setValidator(QDoubleValidator())
        self.lineEdit_samplingFrequency.setValidator(QIntValidator())

        #Define signals
        # Sends the current index (position) of the selected item.
        self.comboBox_typeOfSystem.currentIndexChanged.connect(self.update_selectedSystem)
        self.radioButton_singleFile.toggled.connect(self.update_selectedProcessing)
        self.radioButton_multipleFile.toggled.connect(self.update_selectedProcessing)
        self.spinBox_channelSelection.valueChanged.connect(self.update_desiredChannel)
        self.lineEdit_accelConvFactor.textChanged.connect(self.update_calibrationFactor)
        self.lineEdit_samplingFrequency.textChanged.connect(self.update_samplingFrequencyOriginal)

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
            self.samplingFrequencyOriginal=float(new_value)
        except Exception as e:
            None
        print(self.samplingFrequencyOriginal) if debugActivated else None
    
def main():
    app = QtWidgets.QApplication(sys.argv)
    main = MainWindow()
    main.show()
    sys.exit(app.exec())

if __name__ == '__main__':
    main()
# ------------------------------------------------------
# -------------------- mplwidget.py --------------------
# ------------------------------------------------------
from PyQt6.QtWidgets import*

from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg

from matplotlib.figure import Figure

from PyQt6.QtGui import QFont, QFontInfo

import matplotlib.pyplot as plt

import matplotlib
matplotlib.use('Qt5Agg')

import copy
    
class MplWidget(QWidget):
    
    def __init__(self, parent = None):
        super(QWidget, self).__init__()
        #QWidget.__init__(self, parent)
        
        self.canvas = FigureCanvasQTAgg(Figure())
        
        vertical_layout = QVBoxLayout()
        vertical_layout.addWidget(self.canvas)
        vertical_layout.setContentsMargins(0,0,0,0)
        
        self.canvas.axes = self.canvas.figure.add_subplot(111)
        self.canvas.axes.axis('off')
        self.canvas.draw()
        
        #Get the default font style from the system and apply to the canvas
        default_font = QFont()
        default_font_info = QFontInfo(default_font)
        self.system_font_family = default_font_info.family()
        self.system_font_size = default_font_info.pointSize()
        plt.rcParams.update({'font.family': self.system_font_family})
        plt.rcParams.update({'font.size': self.system_font_size})

        self.setLayout(vertical_layout)

    
    def update_figure(self, new_figure):
        widget_width,widget_height=self.canvas.figure.get_size_inches()
        original_dpi = self.canvas.figure.get_dpi()
        # Clear the existing figure content
        self.canvas.figure.clf()
        self.canvas.axes.clear()

        # Copy the contents of the new figure onto the canvas
        self.canvas.figure = copy.deepcopy(new_figure)
        #self.canvas.figure = new_figure

        # Set the size of the new figure based on the widget's size
        self.canvas.figure.set_size_inches(widget_width, widget_height)
        self.canvas.figure.set_dpi(original_dpi)
        self.canvas.figure.tight_layout()
        
        #Set font
        all_axes = []
        for child in self.canvas.figure.get_children():
            if isinstance(child, plt.Axes):
                all_axes.append(child)
        for ax in all_axes:
            for item in [ax.xaxis.label, ax.yaxis.label, ax.title]:
                item.set_fontfamily( self.system_font_family)
                item.set_fontsize(self.system_font_size)
        #self.canvas.setMinimumSize(widget_width,widget_height)  # Adjust the scale factor as needed

        # Redraw the canvas
        self.canvas.draw()
        plt.close('all')
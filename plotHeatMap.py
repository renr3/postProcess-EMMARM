import matplotlib.cm as cmap
import matplotlib.pyplot as plt
import numpy as np 

from matplotlib.colors import LogNorm
from matplotlib.ticker import LogFormatterMathtext

def generateHeatMap(heatMapFilePath, saveFilePath, unitForAgeInPlot):
    """
    Read the .npz file containing heat map information, generates a heat map plot and save it in .svg to the save path informed
    Parameters
    -------       
    heatMapFilePath : str
        Complete path of the heat map .npz file
    saveFilePath: str
        Complete path in which the .svg figure should be saved
    unitForAgeInPlot: str
        Defines the age unit in the plot. Options are:
            - 'seconds' (default from heat map)
            - 'minutes'
            - 'hours'
            - 'days'
    selectedSystem: str
        String that defines the type of system selected, so the file can be correctly read

    Returns
    -------    
    ageAtMeasurement : int
        Age in seconds at the instant of start of the current measurement session under consideration.
    """ 
    #Convert from seconds to days
    if unitForAgeInPlot == 'seconds':
        ageConversionFactor=1
    elif unitForAgeInPlot == 'minutes':
        ageConversionFactor=60
    elif unitForAgeInPlot == 'hours':
        ageConversionFactor=60*60
    elif unitForAgeInPlot == 'days':
        ageConversionFactor=60*60*24
    else:
        print("Undefined age unit")

    heatMapNPZ = np.load(heatMapFilePath)
    heatMap=heatMapNPZ['arr_0']
    fig, ax = plt.subplots(figsize=(5,4))
    ax.set_xscale('log')
    im=ax.pcolormesh(heatMap[0,1:]/ageConversionFactor,heatMap[1:,0],heatMap[1:,1:], norm=LogNorm())
    ax.set_xlim([0.01,heatMap[0,-1]/ageConversionFactor])
    fig.colorbar(im, ax=ax)
    ax.set_xlabel("Ages (day)")
    ax.set_ylabel("Frequency (Hz)")
    plt.show()

    plt.savefig(saveFilePath, format='svg')
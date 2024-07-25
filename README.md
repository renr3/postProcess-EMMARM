# EMM-ARM post-process software
 
This is the official repository of the EMM-ARM post-process software.

This software allows to post process data obtained from EMM-ARM tests, with some specific test systems: the original implementation by José Granja [[1]](#1), a Raspberry Pi-based system [[2]](#2), a old microcontroller-based system named old_μEMMARM [[3]](#3), and a minimalist system named μEMMARM (check it here: https://github.com/renr3/microEMMARM). The software is divided in two components. 

The first component, intended for modal analysis, reads the results files storing the acceleration time series measured during testing from each of these systems, performs the required pre-processing and modal identification procedures, and outputs the following results in ".txt" files: the frequency and damping evolution along an EMM-ARM test for every modal identification method performed. 

The second component, intended for elastic modulus computation, computes the elastic modulus evolution computed from the outputs of the previous modulus. So the output ".txt" file with the frequency evolution is read, and from physical information about the test informed by the user, it estimates the elastic modulus evolution. As an output, it exports a ".txt" file with the elatic modulus evolution.

You have some options to run this software:

1. GUI_version: this is a version with a GUI. To use it, please insteall the ".exe" file found in the folder "GUI_version". At this same folder, you can check the source files of the software.
2. Jupyter_notebook: this is a Jupyter notebook version. To use it, use your favorite IDE to open the ".ipynb" file, which is already set with an example on how to use the functions of the software directly from Python
3. Direct in Python: you can simply import the Python libraries to your own code and use the functions from there. You should need all the files in the folder "Jupyter_notebook", except for the ".ipynb" file.
This code was developed on Python 3.8.10. It was not tested on other versions, so it may or may not work in different versions. To avoid any versioning errors, run it in an environment with Python 3.8.10 and the modules listed in the requirements.txt file. For installing a virtual environment using the requirements file, check the appropriate tutorial for your IDE. Example for VSCode: https://code.visualstudio.com/docs/python/environments

In the folder "**exampleFiles**" you can find exemplificative files to test the software and see how are the outputs.

This project is greatly supported by the CESSIPy module, about which you can find more information here:
- https://www.sciencedirect.com/science/article/pii/S2352711022000632
- https://github.com/MatheusCarini/CESSIPy

## References
<a id="1">[1]</a>
J. Granja and M. Azenha, “EMM-ARM User’s Guide - version 2.0.1,” Continuous characterization of stiffness of cement-based materials: experimental analysis and micro-mechanics modelling. Accessed: Dec. 26, 2023. [Online]. Available: https://repositorium.sdum.uminho.pt/handle/1822/42563

<a id="2">[2]</a>
T. Russo, R. R. Ribeiro, A. Araghi, R. de M. Lameiras, J. Granja, and M. Azenha, “Continuous Monitoring of Elastic Modulus of Mortars Using a Single-Board Computer and Cost-Effective Components,” Buildings, vol. 13, no. 5, p. 1117, Apr. 2023, doi: 10.3390/buildings13051117.

<a id="3">[3]</a>
R. Rocha Ribeiro, M. I. C. Sousa, J. H. da S. Rêgo, and R. de M. Lameiras, “Innovative low-cost system for early age E-modulus monitoring of cement pastes: validation and application to nanosilica-added and limestone-calcined clay cements,” Mater Struct, vol. 55, no. 1, p. 13, Jan. 2022, doi: 10.1617/s11527-021-01849-w.

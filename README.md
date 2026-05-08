**Radial Velocity Cross Correlator for LCO**

This repo contains code that will recover the radial velocity results from 1D echelle spectrograph data in cases where large-sclae standard pipelines fail. The code does this though cross-correlation of the observation where the pipeline failed (observed file) against an observation where the pipeline succeeded (template file). It is particularly deswigned for use on the Las Cumbres Observatory data products, but with some modifications to the reading of the data files based on their structures, it can be applied to other 1D spectra.

**Running the Code**

The following files need to be installed to the local device: ccf.py, data_prep.py, and radialVelocityCalc.py. The following common Python modules also need to be installed: scipy, numpy, astropy, and  matplotlib. 

The code can be run from the command line by running the radialVelocityCalc.py file with python. The file takes two arguments: the template file location and the observed file location:  python .\radialVelocityCalc.py [templateFileLocation], [observedFileLocation]

The files need to have headers that contain the barycentric velcocity as 'BARYCORR', and the template file must have a header that states its radial velcoity listed as 'RV'.

The program will print the radial velocity of the observed file and the associated error (both in m/s). It will also save a graph of the cross correlation values plotted against its best Gaussian fit as ccf_plot.png.

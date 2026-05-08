import sys
import ccf
import data_prep as dp
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
def getHeaderValues(observedFile, templateFile):
    hdulObserved = fits.open(observedFile)
    hdulTemplate = fits.open(templateFile)
    
    barycorrObserved = hdulObserved[0].header['BARYCORR']
    barycorrTemplate = hdulTemplate[0].header['BARYCORR']
    templateRV = hdulTemplate[0].header['RV']
    
    hdulObserved.close()
    hdulTemplate.close()
    return float(barycorrObserved), float(barycorrTemplate), float(templateRV)




def main():
    if len(sys.argv) != 3:
        print("Wrong number of inputs")
        sys.exit(1)
    templateFileName = sys.argv[1]
    observedFileName = sys.argv[2]

    barycorrObs, barycorrTemp, templateRV = getHeaderValues(observedFileName, templateFileName)
    
    obsWavelengths, obsFluxes = dp.readingAndReadyingData(observedFileName)
    tempWavelengths, tempFluxes = dp.readingAndReadyingData(templateFileName)
    obsWavelengths, obsFluxes, tempWavelengths, tempFluxes = dp.makeOverlap(obsWavelengths, obsFluxes, tempWavelengths, tempFluxes)
    newWavelengths, newObsFluxes, newTempFluxes = dp.interpolation(obsWavelengths, obsFluxes, tempWavelengths, tempFluxes)
    lags, ccfValues = ccf.performCCF(newObsFluxes, newTempFluxes, newWavelengths)
    A, k0, sigma = ccf.fitGaussian(lags, ccfValues)
    rv = ccf.calculateFinalRV(k0, newWavelengths, templateRV, barycorrObs, barycorrTemp)
    error = ccf.rvError(newTempFluxes, newWavelengths, A)
    print("Radial velocity: " + str(rv) + " and error: " + str(error))


    plt.plot(lags, ccfValues, label='CCF')
    plt.plot(lags, ccf.gaussian(lags, A, k0, sigma), label='Gaussian Fit', linestyle='--')
    plt.xlabel("Pixel Lag")
    plt.ylabel("CCF")
    plt.legend()
    plt.title("Cross-Correlation Function")
    plt.savefig("ccf_plot.png")


#Google says that the following is how you make the main function run from command line
if __name__ == "__main__":
    main()
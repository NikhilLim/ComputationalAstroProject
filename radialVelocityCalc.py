import sys
import ccf
import data_prep
def getHeaderValues(observedFile, templateFile):
    hdulObserved = fits.open(observedFile)
    hdulTemplate = fits.open(templateFile)
    
    barycorrObserved = hdulObserved[1].header['BARYCORR']
    barycorrTemplate = hdulTemplate[1].header['BARYCORR']
    templateRV = hdulTemplate[1].header['RV']
    
    hdulObserved.close()
    hdulTemplate.close()
    return float(barycorrObserved), float(barycorrTemplate), float(templateRV)




def main():
    if len(sys.argv) != 3:
        print("Wrong number of inputs")
        sys.exit(1)
    templateFileName = sys.argv[1]
    observedFileName = sys.argv[2]
    obsWavelengths, obsFluxes = readingAndReadyingData(observedFile)
    tempWavelengths, tempFluxes = readingAndReadyingData(templateFile)
    obsWavelengths, obsFluxes, tempWavelengths, tempFluxes = clipToOverlap(obsWavelengths, obsFluxes, tempWavelengths, tempFluxes)
    newWavelengths, newObsFluxes, newTempFluxes = interpolation(obsWavelengths, obsFluxes, tempWavelengths, tempFluxes)
    lags, ccfValues = performCCF(newObsFluxes, newTempFluxes, newWavelengths)
    A, k0, sigma = fitGaussian(lags, ccfValues)
    rv = calculateFinalRV(k0, newWavelengths, templateRV, barycorrObs, barycorrTemp)
    error = rvError(newTempFluxes, newWavelengths, A)
    print("Radial velocity: " + str(rv) + " and error: " + error)


    plt.plot(lags, ccfValues, label='CCF')
    plt.plot(lags, gaussian(lags, A, k0, sigma), label='Gaussian Fit', linestyle='--')
    plt.xlabel("Pixel Lag")
    plt.ylabel("CCF")
    plt.legend()
    plt.title("Cross-Correlation Function")
    plt.savefig("ccf_plot.png")


#Google says that the following is how you make the main function run from command line
if __name__ == "__main__":
    main()
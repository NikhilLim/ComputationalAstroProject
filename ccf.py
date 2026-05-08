from scipy.constants import c
from scipy.optimize import minimize
import numpy as np
def maxLagRange(velocityMax, wavelengthGrid): #this function (taking km/s) will prevent the code from taking literal hours to run by limiting lag range to physically possible values.
    #basing this off equation 7 in proposal
    minWavelength = wavelengthGrid[0]
    maxWavelength = wavelengthGrid[-1]
    rvMax = 1000*velocityMax
    numerator = np.log(rvMax/c+1)
    denominator = (np.log(maxWavelength) - np.log(minWavelength))/ (wavelengthGrid.size-1)
    return int(numerator/denominator)

def performCCF(observedFluxes, templateFluxes, wavelength):
    lagRange = maxLagRange(200,wavelength)
    lags = np.arange(-lagRange, lagRange + 1)
    ccfValues = np.zeros(len(lags))

    for j in range(len(observedFluxes)):
        for index, k in enumerate(lags):
            if 0 <= j-k < len(observedFluxes): #ensures array calls stay in bounds
                ccfValues[index] += observedFluxes[j] * templateFluxes[j-k]
    return lags, ccfValues


#from notes
def least_squares_fit(x, y, model, guess):


    def statistic(params):
        return np.sum((y - model(x, *params))**2)

    res = minimize(statistic, guess)

    return res.x
def gaussian(k, A, k0, sigma, C):
    return A * np.exp(-((k - k0)**2) / (2 * sigma**2)) + C
def fitGaussian(lags, ccfValues):
    guess = np.array([np.max(ccfValues), lags[np.argmax(ccfValues)], 10, np.median(ccfValues)]) # guess for A, k0, sigma
    A, k0, sigma, C = least_squares_fit(lags, ccfValues, gaussian, guess)
    return A, k0, sigma, C
def rvError(templateFluxes, wavelengthGrid, A):
    rMax = A / np.sum(templateFluxes**2)
    N = wavelengthGrid.size
    lnLambda = np.log(wavelengthGrid)
    sigmaTemplate  = np.std(templateFluxes)

    #this part calculates the derivative dFtemp/d(ln lambda)
    derivative = np.zeros(N)
    derivative[0] = (templateFluxes[1] - templateFluxes[0]) / (lnLambda[1] - lnLambda[0])
    derivative[-1] = (templateFluxes[-1] - templateFluxes[-2]) / (lnLambda[-1] - lnLambda[-2])
    for i in range(1, N-1):
        derivative[i] = (templateFluxes[i+1] - templateFluxes[i-1]) / (lnLambda[i+1] - lnLambda[i-1])
    squaredDerivative = np.mean(derivative**2)
    error = c*sigmaTemplate/np.sqrt(squaredDerivative)
    error = error * (np.sqrt((1/N)*(1/rMax**2 - 1)))
    return error

def calculateFinalRV(k0, wavelengthGrid, templateRV, barycorrObs, barycorrTemp):
    minWavelength = wavelengthGrid[0]
    maxWavelength = wavelengthGrid[-1]
    
    logStep = (np.log(maxWavelength) - np.log(minWavelength)) / (wavelengthGrid.size - 1)
    
    rv = c * (np.exp(k0 * logStep) - 1)
    
    # Add template RV 
    rv += templateRV
    
    # barycentric correction
    rv += (barycorrObs - barycorrTemp)
    
    return rv
    
from scipy.constants import c
from scipy.optimize import minimize
import numpy as np
def maxLagRange(velocityMax, wavelengthGrid): 

    """
    This function will prevent the code from taking literal hours to run by limiting lag range to physically possible values.

    Args:
        velocityMax (int or float): maximum velocities for the cross correlation to later consider.
        wavelengthGrid (numpy array): interpolated grid that the cross corelation will be done over.

    Returns:
        int: the max lag the cross correlation should look at.
    """
    

    #basing this off equation 7 in proposal
    minWavelength = wavelengthGrid[0]
    maxWavelength = wavelengthGrid[-1]
    rvMax = 1000*velocityMax
    numerator = np.log(rvMax/c+1)
    denominator = (np.log(maxWavelength) - np.log(minWavelength))/ (wavelengthGrid.size-1)
    return int(numerator/denominator)

def performCCF(observedFluxes, templateFluxes, wavelength):
    """
    Performs cross correlation

    Args:
        observedFluxes (numpy array): observed flux values (interpolated)
        templateFluxes (numpy array): template flux values (interpolated)
        wavelength (numpy array): interpolated grid that the cross corelation will be done over.

        Returns:
        int: the max lag the cross correlation should look at.
    """
    lagRange = maxLagRange(200,wavelength)
    lags = np.arange(-lagRange, lagRange + 1)
    ccfValues = np.zeros(len(lags))
    for j in range(len(observedFluxes)): #loops over the flux values
        for index, k in enumerate(lags): #loops through the lags and their indices
            if 0 <= j-k < len(observedFluxes): #ensures array calls stay in bounds
                ccfValues[index] += observedFluxes[j] * templateFluxes[j-k] #calculation for ccf. Adds it to the total, which for each lag will be for each pixel
    return lags, ccfValues


#stole this from class notes
def least_squares_fit(x, y, model, guess):
    """
    Perform least squares fitting using scipy.optimize.minimize.

    Parameters:
    x (numpy array): Independent variable.
    y (numpy array): Dependent variable.
    model : Model function to fit the data.
    guess (numpy array): Initial guess for the parameters of the model.

    returns the result
    """

    #defines what is to be minimized
    def statistic(params):
        return np.sum((y - model(x, *params))**2)
    #does the minimizing
    res = minimize(statistic, guess)

    return res.x
def gaussian(k, A, k0, sigma, C):
    #this function just defines the Gaussian equaion that the data will be fit to
    return A * np.exp(-((k - k0)**2) / (2 * sigma**2)) + C
def fitGaussian(lags, ccfValues):
    """
    Makes a guess and calls the least square fitting on the gaussian and the data

    Parameters:
    lags (numpy array): the difference wavelength shifts
    ccfValues (numpy array): the associated cross correlation sums


    returns the Gaussian parameters of the fit
    """
    guess = np.array([np.max(ccfValues), lags[np.argmax(ccfValues)], 10, np.median(ccfValues)]) # guess for A, k0, sigma
    A, k0, sigma, C = least_squares_fit(lags, ccfValues, gaussian, guess)
    return A, k0, sigma, C
def rvError(observedFluxes, templateFluxes, wavelengthGrid, A):
    """
    Calculates the radial velcoity error based on the formula from the proposal

    Parameters:
        observedFluxes (numpy array): observed flux values (interpolated)
        templateFluxes (numpy array): template flux values (interpolated)
        wavelengthGrid (numpy array): interpolated wavelength grid that the fluxes are currently mapped to
        A (float): parameter from the gaussian fit


    returns the error value as float
    """
    rMax = A / np.sqrt(np.sum(observedFluxes**2) * np.sum(templateFluxes**2))
    N = wavelengthGrid.size
    lnLambda = np.log(wavelengthGrid)
    sigmaTemplate  = np.std(templateFluxes)

    #this part calculates the derivative dFtemp/d(ln lambda)
    derivative = np.zeros(N)
    #these next two lines are for the edges
    derivative[0] = (templateFluxes[1] - templateFluxes[0]) / (lnLambda[1] - lnLambda[0])
    derivative[-1] = (templateFluxes[-1] - templateFluxes[-2]) / (lnLambda[-1] - lnLambda[-2])
    #central approximation of derivative
    for i in range(1, N-1):
        derivative[i] = (templateFluxes[i+1] - templateFluxes[i-1]) / (lnLambda[i+1] - lnLambda[i-1])
    #calculatios according to formula
    squaredDerivative = np.mean(derivative**2)
    error = c*sigmaTemplate/np.sqrt(squaredDerivative)
    error = error * (np.sqrt((1/N)*(1/rMax**2 - 1)))
    return error

def calculateFinalRV(k0, wavelengthGrid, templateRV, barycorrObs, barycorrTemp):
    """
    Does the final radial velocity calculation

    Parameters:
        k0 (float): result from the gaussian fit. Should be the shift at which the cross correlation value was highest/ the fluxes overlapped the most
        wavelengthGrid (numpy array): interpolated wavelength grid that the fluxes are currently mapped to
        templateRV (float): known radial velocity of template from header
        barycorrObs (float): known barycentric velocity of observed spectrum from header
        barycorrTemp(float): known barycentric velocity of template from header


    returns the rv value as float
    """
    #based on the RV formula from proposal
    minWavelength = wavelengthGrid[0]
    maxWavelength = wavelengthGrid[-1]
    
    logStep = (np.log(maxWavelength) - np.log(minWavelength)) / (wavelengthGrid.size - 1)
    
    rv = c * (np.exp(k0 * logStep) - 1)
    
    # Add template RV 
    rv += templateRV
    
    # barycentric correction
    rv += (barycorrObs - barycorrTemp)
    
    return rv
    
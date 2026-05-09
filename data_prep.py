from astropy.io import fits
import numpy as np
def sigmaclipping(wavelengths,fluxes , alpha, tolerance, maxIterations):
    """
    removes peaks from wavelengths

    Parameters:
        wavelengths (numpy array): wavelength values
        fluxes (numpy array): the associated flux at eacdh wavelength
        alpha (float): user defined threshold
        tolerance (float): if the standard deviation change between iteations falls below this value, stop the clipping
         maxIterations (int): max loops


    returns the wavelengths and fluxes but without their peaks
    """
    #based on the RV 
    clippedFluxes = fluxes.copy()
    clippedWavelengths = wavelengths.copy()
    iteration = 0
    #loop
    while(iteration < maxIterations):
        median = np.median(clippedFluxes)
        std = np.std(clippedFluxes)
        #clipping condition
        clippingMask = np.abs(clippedFluxes - median) < alpha * std
        newClippedFluxes = clippedFluxes[clippingMask]
        newClippedWavelengths = clippedWavelengths[clippingMask]
        #sees if standard deviation isn't really changing
        if(std -  np.std(newClippedFluxes) < tolerance):
            clippedFluxes = newClippedFluxes
            clippedWavelengths = newClippedWavelengths
            break
        #makes the newclipped fluxes the main ones and startsthe loop again
        clippedFluxes = newClippedFluxes
        clippedWavelengths = newClippedWavelengths
        
        iteration = iteration +1
        
    return clippedWavelengths, clippedFluxes
def readingAndReadyingData(fileAddress):
    #this just takes a file address and reads that data into two numpy arrays of wavelengtha nd fluxes
    hdul = fits.open(fileAddress)
    hdul.info()
    spec = hdul[1]

    data = hdul[1].data
    scifiber = hdul[0].header['SCIFIBER']

    
    allWavelengths = []
    allNormFlux = []
    #goes through all rows because this is echelle spectro
    for row in data:
        #makes sure its science data ad good data
        if row['fiber'] != scifiber:
            continue
        m = row['mask'] == 0
        good = (
        m &
        np.isfinite(row['normflux']) &
        np.isfinite(row['wavelength']) &
        (row['blaze'] > 0)
        )
        wavelengthsOfRow, fluxesOfRow = sigmaclipping(row['wavelength'][good],row['normflux'][good], 3, 0.000001, 15)
        allWavelengths.append(wavelengthsOfRow)
        allNormFlux.append(fluxesOfRow)
    
    
        #print(row['blaze'])
        #ax.plot(row['wavelength'][good], (row['normflux'][good]))
    #combining and ordering everything
    allWavelengths = np.concatenate(allWavelengths)
    allNormFlux = np.concatenate(allNormFlux)
    ordering = np.argsort(allWavelengths)
    allWavelengths = allWavelengths[ordering]
    allNormFlux = allNormFlux[ordering]

    return allWavelengths, allNormFlux

def makeOverlap(observedWavelengths, observedFluxes, templateWavelengths, templateFluxes):
    #simple function that gets rid of any values for which the wavelengths of the two data sets don't overlap
    minWavelength = max(observedWavelengths[0], templateWavelengths[0])
    maxWavelength = min(observedWavelengths[-1], templateWavelengths[-1])
    observedMask = (observedWavelengths >= minWavelength) & (observedWavelengths <= maxWavelength)
    templateMask = (templateWavelengths >= minWavelength) & (templateWavelengths <= maxWavelength)
    return observedWavelengths[observedMask], observedFluxes[observedMask], templateWavelengths[templateMask], templateFluxes[templateMask]
    

def interpolation(observedWavelengths, observedFluxes, templateWavelengths, templateFluxes):
    """
    moves both flux arrays ontoa a common, logarithimic wavelength grid

    Parameters:
        observedWavelengths (numpy array): observed spectrum wavelengths
        observedFluxes (numpy array): associated flux values
        templateWavelengths (numpy array): template spectrum wavelengths
        templateFluxes (numpy array): associated flux values


    returns both the fluxes with only one wavelength grid that they share
    """
    N = observedWavelengths.size

    minWavelength = observedWavelengths[0]
    maxWavelength = observedWavelengths[-1]
    j = np.arange(N)
    #calculating the new wavelength grid that is logarithimic
    exponentMultiplier = (np.log(maxWavelength)-np.log(minWavelength))/(N-1)
    newWavelengths = minWavelength * np.exp(j*exponentMultiplier)
    #loops through all obsered fluxes
    newObservedFluxes = np.zeros(N)
    i=0
    for j in range(N):
        if j == N - 1:  #edge case of last point
            newObservedFluxes[j] = observedFluxes[-1]
            break
        while i < N - 2 and observedWavelengths[i+1] < newWavelengths[j]: #what i is doing, is its basically finding what original wavelength the new wavelength value is right before. The orginal wavelength values that the new wvaelngth value is in between is used in the new flux calculation
            i += 1 
        newObservedFluxes[j] = observedFluxes[i] + ((observedFluxes[i+1]-observedFluxes[i])/(observedWavelengths[i+1]-observedWavelengths[i])) * (newWavelengths[j]-observedWavelengths[i]) #does the calculation for the flux that is associated with the new grid and puts it in that spot


    #literally just the code above but for template values
    newTemplateFluxes = np.zeros(N)
    i=0
    for j in range(N):
        if j == N - 1:  #edge case of last point
            newTemplateFluxes[j] = templateFluxes[-1]
            break
        while i < N - 2 and templateWavelengths[i+1] < newWavelengths[j]:
            i += 1 
        newTemplateFluxes[j] = templateFluxes[i] + ((templateFluxes[i+1]-templateFluxes[i])/(templateWavelengths[i+1]-templateWavelengths[i])) * (newWavelengths[j]-templateWavelengths[i])
    return newWavelengths, newObservedFluxes, newTemplateFluxes
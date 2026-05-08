def sigmaclipping(wavelengths,fluxes , alpha, tolerance, maxIterations):
    clippedFluxes = fluxes.copy()
    clippedWavelengths = wavelengths.copy()
    iteration = 0
    while(iteration < maxIterations):
        median = np.median(clippedFluxes)
        std = np.std(clippedFluxes)
        clippingMask = np.abs(clippedFluxes - median) < alpha * std
        newClippedFluxes = clippedFluxes[clippingMask]
        newClippedWavelengths = clippedWavelengths[clippingMask]
        
        if(std -  np.std(newClippedFluxes) < tolerance):
            clippedFluxes = newClippedFluxes
            clippedWavelengths = newClippedWavelengths
            break
            
        clippedFluxes = newClippedFluxes
        clippedWavelengths = newClippedWavelengths
        
        iteration = iteration +1
        
    return clippedWavelengths, clippedFluxes
def readingAndReadyingData(fileAddress):
    hdul = fits.open(fileAddress)
    hdul.info()
    spec = hdul[1]
    print(spec.columns)
    data = hdul[1].data
    scifiber = hdul[0].header['SCIFIBER']
    print(scifiber)
    fig, ax = plt.subplots(figsize=(12, 5))
    
    allWavelengths = []
    allNormFlux = []
    for row in data:
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
    
    
        np.set_printoptions(threshold=sys.maxsize)
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
    minWavelength = max(observedWavelengths[0], templateWavelength[0])
    maxWavelength = min(observedWavelengths[-1], templateWavelength[-1])
    observedMask = (observedWavelengths >= minWavelength) & (observedWavelengths <= maxWavelength)
    templateMask = (templateWavelengths >= minWavelength) & (templateWavelengths <= maxWavelength)
    return observedWavelengths[observedMask], observedFluxes[observedMask], templateWavelengths[templateMask], templateFluxes[templateMask]
    

def interpolation(observedWavelengths, observedFluxes, templateWavelengths, templateFluxes):
    N = observedWavelengths.size
    print(N)
    minWavelength = observedWavelengths[0]
    maxWavelength = observedWavelengths[-1]
    j = np.arange(N)
    exponentMultiplier = (np.log(maxWavelength)-np.log(minWavelength))/(N-1)
    newWavelengths = minWavelength * np.exp(j*exponentMultiplier)
    newObservedFluxes = np.zeros(N)
    i=0
    for j in range(N):
        if j == N - 1:  #edge case of last point
            newObservedFluxes[j] = observedFluxes[-1]
            break
        while i < N - 2 and observedWavelengths[i+1] < newWavelengths[j]: #what i is doing, is its basically finding what original wavelength the new wavelength value is right before. The orginal wavelength values that the new wvaelngth value is in between is used in the new flux calculation
            i += 1 
        newObservedFluxes[j] = observedFluxes[i] + ((observedFluxes[i+1]-observedFluxes[i])/(observedWavelengths[i+1]-observedWavelengths[i])) * (newWavelengths[j]-observedWavelengths[i])

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
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
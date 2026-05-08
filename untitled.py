from scipy.constants import c
def maxLagRange(velocityMax, wavelengthGrid): #this function (taking km/s) will prevent the code from taking literal hours to run by limiting lag range to physically possible values.
    #basing this off equation 7 in proposal
    rvMax = 1000*velocityMax
    numerator = np.log(rvMax/c+1)
    denominator = (N-1)/(

def performCCF(observedFlux, templateFlux, wavelength):
    
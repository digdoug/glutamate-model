from scipy.special import erf

# Connection parameters
Jp_eeA         = 1.64         # maximum connection weight from pyramidal to pyramidal
sigma_eeA    = 14.4         # width of the connectivity footprint from pyramidal to pyramidal

Jp_eeB         = 1.73          # maximum connection weight from pyramidal to pyramidal
sigma_eeB    = 12.76         # width of the connectivity footprint from pyramidal to pyramidal

def circ_distance(fi):
    y = fi
    if (0 <= fi <= 180.):
        y = fi
    elif (180. < fi <= 360.):
        y = -360. + fi
    elif (-180. < fi < 0):
        y = fi
    elif (-360. <= fi < -180.):
        y = 360 + fi
    elif (fi > 360.):
        y = circ_distance(fi - 360.)
    elif (fi < -360.):
        y = circ_distance(fi + 360.)
    return y


# connectivity footprint
tmpA   = sqrt(2*pi)*sigma_eeA*erf(360.*0.5/sqrt(2.)/sigma_eeA)/360.
Jm_eeA = (1.-Jp_eeA*tmpA)/(1.-tmpA)


weight_fnA=lambda delta: (Jm_eeA+(Jp_eeA-Jm_eeA)*exp(-0.5*((circ_distance(delta))**2)/sigma_eeA**2))

weight_cross= lambda theta_i,theta_j,sigma: 1./(sigma*(2*pi)**0.5)*exp(-0.5*(circ_distance(theta_i-theta_j)**2./sigma**2.))

weight_eA=zeros(NE)
weight_eB=zeros(NE)

for i in xrange(NE): 
    weight_eA[i]=weight_fnA(360./NE*i)
    weight_eB[i]=1.


fweightA = rfft(weight_eA) # Fourier transform


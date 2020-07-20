import numpy as np
from scipy.optimize import minimize

def f_cmf(M,R,si=0,fe=0.1):
    """
    Calculate core mass fraction (cmf) for a planet given
    its observed mass, radius and assuming the silica and iron 
    amount in the core and mantle respectivly (by mol). 
    Note:
        For nominal case one should assume 0.1 iron and 0 silica.
    Parameters
    -----
        M: float/int
            Mass of the planet in Earth masses.
        R: float/int
            Radius of the planet in Earth radius.
        si: float
            Amount of silica in core by mol.
        fe: float
            Amount of iron in mantle by mol.
    Returns: float
        Core mass fraction for the planet, as a decimal (Mcore/Mtotal)
    """
    cnm = np.array([ [0.360172,3.071785],
        [-11.164546,-2.225710],
        [6.686932,-0.806478] ])
    alpha = 0.55076
    beta = np.array([1.361339,-3.039501])
    a = -0.567941
    b = -2.976432
    m = 0.263157
    c = 0.031716
    args = np.array([alpha,beta,cnm,a,b,c,m])
    
    M = np.log10(M)
    R = np.log10(R)
    fe = fe-0.1
    x0 = (M/m+R-c)/(m+1/m)
    y0 = x0*m+c
    lam = np.sqrt((x0-M)**2+(y0-R)**2)
        
    res = 0
    for i in range(3):
        for j in range(2):
            res += R**i*M**j*cnm[i,j]
    return res*np.exp(lam*b)*(si*beta[1]+si**2*beta[0]+1)*(fe*alpha+1)+fe*a

def f_rho0(cmf,rhoc=8278,rhom=4000):
    '''
    Parameters
    -----
        cmf: float
            Core mass fraction.
        rhoc: float/int
            Reference density of core.
        rhom: float/int
            Reference density of mantle.
    Returns: float
        Uncompressed density in [kg/m^3].
    '''
    return rhoc*rhom/(rhoc+cmf*(rhom-rhoc))

def f_FeSi(cmf,si=0,fe=0.1,py=0.6):
    '''
    Calculate refactory chemical ratio of Fe/Si for the planet,
    given its cmf.
    Parameters
    -----
        cmf: float
            Core mass fraction.
        si: float
            Amount of silica in core by mol.
        fe: float
            Amount of iron in mantle by mol.
        py: float
            Amount of olivine and pyroxene in the mantle.
    Returns: float
        Fe/Si ratio by weight.
    '''
    Fe,Ni,Si,Mg,O = [55.85e-3,58.69e-3,28.09e-3,24.31e-3,16e-3]
    km = cmf*(2*(1-fe)*Mg+2*fe*Fe+Si*(1+py)+2*O*(2+py))
    kc = (1-cmf)*((0.88-si)*Fe+0.1*Ni+si*Si)  
    return Fe/Si*( (0.88-si)*km+2*fe*kc )/( si*km+(1+py)*kc )

def star_to_planet(FeSi,si=0,fe=0.1,py=0.6):
    """
    Convert Fe/Si values of the stars to cmf of planet depending
    on the assumed structure.
    Parameters
    -----
        FeSi: array_like
            Fe/Si ratio by weight of the stars.
        si: float
            Amount of silica in core by mol.
        fe: float
            Amount of iron in mantle by mol.
        py: float
            Amount of olivine and pyroxene in the mantle.
    Returns: float
        Fe/Si ratio by weight.
    """
    func = lambda cmf,FeSi: (f_FeSi(cmf,si,fe,py)-FeSi)**2
    cmf = np.zeros(len(FeSi))
    for i,item in enumerate(FeSi):
        res = minimize(func,0.5,args=item)
        cmf[i] = res.x
    return cmf

def plot_cont(cmf):
    pdf,bins = np.histogram(cmf,int(np.sqrt(len(cmf))),
                        range=(0,max(cmf)),density=True)
    cmf_x = (bins[:-1]+bins[1:])/2
    N = len(cmf_x)
    X1 = np.zeros([N,N])
    Y1 = np.zeros([N,N])
    Z1 = np.zeros([N,N])
    x = np.linspace(-0.5,1.6,N) #mass in log space
    m = 0.2544
    args = np.array([-0.08996526, -0.18532778,  1.0564978 ])
    for i in range(N):
        c = np.log10(args[2]+args[1]*cmf_x[i]+args[0]*cmf_x[i]**2)
        y = m*x+c+0.005
        X1[i] = x
        Y1[i] = y
        Z1[i] += pdf[i]
    return 10**X1,10**Y1,Z1
    
    
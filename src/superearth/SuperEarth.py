import numpy as np
from scipy.optimize import minimize

def f_cmf(M,R,si=0,fe=0.1):
    '''
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
            Radius of the planet in Earth radii.
        si: float
            Amount of silica in core by mol.
        fe: float
            Amount of iron in mantle by mol.

    Returns: float
        Core mass fraction for the planet, as a decimal (Mcore/Mtotal)

    '''
    cnm = np.array([ [0.360172,3.071785],
        [-11.164546,-2.225710],
        [6.686932,-0.806478] ])
    alpha = 0.55076
    beta = np.array([1.361339,-3.039501])
    a = -0.567941
    b = -2.976432
    m = 0.263157
    c = 0.031716

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
    Calculate uncompressed density for the planet,
    given its cmf.

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
    return Fe/Si*((0.88-si)*km+2*fe*kc)/(si*km+(1+py)*kc)

def f_FeMg(cmf,si=0,fe=0.1,py=0.6):
    '''
    Calculate refactory chemical ratio of Fe/Mg for the planet,
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
        Fe/Mg ratio by weight.

    '''
    Fe,Ni,Si,Mg,O = [55.85e-3,58.69e-3,28.09e-3,24.31e-3,16e-3]
    km = cmf*(2*(1-fe)*Mg+2*fe*Fe+Si*(1+py)+2*O*(2+py))
    kc = (1-cmf)*((0.88-si)*Fe+0.1*Ni+si*Si)  
    return Fe/Mg*((0.88-si)*km+2*fe*kc)/( 2*(1-fe)*kc )

def st_pl(FeX,si=0,fe=0.1,py=0.6,method='Fe2Mg',f_FeX=None):
    '''
    Convert Fe/X values of the stars to cmf of planet depending
    on the assumed structure.
    
    Parameters
    -----
        FeX: float/array_like
            Fe/X ratio by weight of the stars.
        si: float
            Amount of silica in core by mol.
        fe: float
            Amount of iron in mantle by mol.
        py: float
            Amount of olivine and pyroxene in the mantle.
        method: string
            Specify the ratio, Fe2Mg or Fe2Si.
        f_FeX: function
            Pass a function to convert f_FeX.
    
    Return: array_like
        cmf of converted stellar Fe/X ratio.

    '''
    if not isinstance(FeX, (list, tuple, np.ndarray)): FeX = [FeX]
    FeX = np.asarray(FeX) 
    if method=='Fe2Mg':
        f_FeX = f_FeMg
    elif method=='Fe2Si':
        f_FeX = f_FeSi
        
    func = lambda cmf,FeX: (f_FeX(cmf,si,fe,py)-FeX)**2
    cmf = np.zeros(len(FeX))
    for i,item in enumerate(FeX):
        res = minimize(func,0.5,args=item)
        cmf[i] = res.x
    return cmf

def guess_pl(x,cmf,si,fe,py,method):
    if method=='R':
        def func(R):
            if 0.8<R<5.:
                return (f_cmf(x,R,si,fe)-cmf)**2
            else:
                return 1e6
    elif method=='M':
        def func(M):
            if 0.8<M<25:
                return (f_cmf(M,x,si,fe)-cmf)**2
            else:
                return 1e6
    res = minimize(func,2)
    return res.x[0]
def guess_R(M,FeSi,si=0,fe=0.1,py=0.6,cmf=None):
    '''
    Guess planets radius [Re] given its mass [Me], based on estimate of 
    Fe/Si ratio of the star or cmf.
    '''
    if not cmf:
        cmf = float(st_pl(FeSi,si,fe,py))
    R = guess_pl(M,cmf,si,fe,py,'R')
    return R
def guess_M(R,FeSi,si=0,fe=0.1,py=0.6,cmf=None):
    '''
    Guess planets mass [Me] given its radius [Re], based on estimate of 
    Fe/Si ratio of the star or cmf.
    '''
    if not cmf:
        cmf = float(st_pl(FeSi,si,fe,py))
    M = guess_pl(R,cmf,si,fe,py,'M')
    return M






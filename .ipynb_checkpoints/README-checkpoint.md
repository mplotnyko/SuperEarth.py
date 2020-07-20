# SuperEarth 

SuperEarth package is an analytical model for apparoximating core mass fraction (cmf M_core/M_total) and other interior parameters for a Super-Earth sized planet (1-20M_e), assuming Earth like composition and interior structure.
The analytical function is a fit to a set of generated data using interior structure code which follows [Valencia et al. (2007)](https://iopscience.iop.org/article/10.1086/509800), for more infromation visit [Plotnykov & Valencia (2020)] work.

## Installation

1. Clone the repository

```git clone https://github.com/mplotnyko/SuperEarth.py```

## Example 

In order to approximate cmf one will need Radius and Mass of the planet, as well as assume the iron (x_{Fe}) amount in the mantle by mol and silica (x_{Si}) amount in the core by mol. 
For the nominal case one should assume 0.1 x_{Fe} (10\%) and 0 x_{Si}.

Using Earth as a test case, we can compare the CMF to the approximate solution given by f_cmf and assuming nominal case.
For this test case the function performs the best, where the worst case is \sim 0.02 
    
    >>> import SuperEarth
    >>> cmf_E = 0.325
    >>> cmf = SuperEarth.f_cmf(1,1,si=0,fe=0.1)
    >>> cmf-cmF_E
    0.0037472894465270246

Once the cmf is know, other useful quentites can be calculated such as the iron to silica refractory ratio (Fe/Si) or the uncompressed density (\rho_0) of the planet.
For example, using simulated data from structural code ([Plotnykov & Valencia (2020)]), one can calculate the uncompressed denisty.

    import numpy as np
    M,R,CMF = np.loadtxt('data.csv',skiprows=1,delimiter=',').T
    cmf = SuperEarth.f_cmf(M,R,si=0,fe=0.1))
    rho_0 = SuperEarth.f_rho0(cmf,rhoc=8278,rhom=4000)

Where rhoc and rhom are the nominal values for core and mantle density under reference conditions. Here is how the data looks like.

![](images/MR_rho0.jpg)

Same procedure can be done to estimate Fe/Si ratio or to use Stellar Fe/Si ratio and convert it to cmf.
For example, assuming some Fe/Si distribtuion for stars and finding the corresponding cmf, we look at distribution of stars given that they are planets. 
Assuming nominal case for the planet composition.

    FeSi = np.random.normal(2,1,500)
    cmf = SuperEarth.star_to_planet(FeSi,si=0,fe=0.1,py=0.6)

Here is the result of such experimnet, where the density plot are converted stars and the points are some simulated planets.
![](images/FeSi_star.jpg)

These examples can be found in the example.py file.



    from LC_GenMock import generate_mock_data, plot_mock_data
    import LC_model
    import numpy as np

The first step is to generate transit data with a Starspot.
Then the duration of transit and the strength of the perturbation due to the Starspot needs to be specified.

    transit_window = 10 #hours
    spot_stength = 2000 #ppm
    spot_width = 0.03 #days
    spot_location = transit_window/24*0.5 #days
    starspot = [transit_window,spot_location,spot_stength,spot_width]

Using quadratic limb darkening model the transit parameters (theta_0 & ld_coef) can be specified to generate the data.

    theta_0 = [140.5,0.0576,1.571,224.22] 
    ld_coef = [0.47,0.19]
    results = generate_mock_data(theta_0,ld_coef,transit_window,*starspot)

These results can be easily ploted using the following command:

    fig = plot_mock_data( *results)

If one wishes then to obtain a best fit results for the generated mock data, here is an example of using [emcee](https://github.com/dfm/emcee).
First get an estimate of the parameters using least-squares method:

    theta = np.array([140,0.1,np.pi/2,224,0.5,0.5])
    time = x_mock
    flux = y_mock*1e-6+1

    fit = LC_model.transit_fit(theta[:-2],theta[-2:],method="quad_ld")

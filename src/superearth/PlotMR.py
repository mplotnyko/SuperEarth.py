import os
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import numpy as np
import pandas as pd
from superearth import st_pl,guess_R
from superearth.utils import check_cache_archive,plot_cont,plot_cmf,MR_H2O,MR_HHe

package_dir = os.path.dirname(__file__)
check_cache_archive()
def exoplanets(Merr, Rerr, default_pl=True, best_pl=False, Mrange=[0, 30], Rrange=[0, 4.5], Teffrange=None,
               onlymultiplanet=False, onlyMplanets=False, rocky=False, query=None):
    """
    Retrieves exoplanet data from the NASA Exoplanet Archive.

    Parameters
    ----------
    Merr: float
        Maximum relative error allowed for planet mass.
    Rerr: float
        Maximum relative error allowed for planet radius.
    default_pl: bool, default: True
        Flag whether to use default NASA planetary data.
    best_pl: bool, default: False
        Falg whether to use the most up to date reference.
        When True ignores default flag, can have slighly different planets.
    Mrange: list, default: [0, 30]
        Range of planet mass [min_mass, max_mass].
    Rrange: list, default: [0, 4.5]
        Range of planet radius [min_radius, max_radius].
    Teffrange: list, default: None
        Range of stellar effective temperatures [min_Teff, max_Teff].
    onlymultiplanet: bool, default: False
        Flag whether to use only multi-planet systems.
    onlyMplanets: bool, default: False
        Flag whether to use only M dwarf systems.
    rocky: bool, default: True
        Flag whether to use only rocky planets.
    query: str, default: None
        Additional custom query string to filter the data.

    Returns
    -------
    class: pandas.DataFrame
        pandas DataFrame object containing the filtered exoplanet data.

    """
    df_exo = pd.read_csv(package_dir+"/Data/ExoplanetArchiveData.csv")    
    all_query = f'{Mrange[0]}<pl_masse<{Mrange[1]} and {Rrange[0]}<pl_rade<{Rrange[1]} '+\
                f'and pl_masseerr1/pl_masse<{Merr} and pl_radeerr1/pl_rade<{Rerr} and  '+\
                f'pl_masselim>-1 and pl_radelim>-1'
    if best_pl:
        df_exo.query('pl_masse>0 and pl_rade>0 and pl_masselim>-1 and pl_radelim>-1',inplace=True)
        df_exo.drop_duplicates(subset=['pl_name'],inplace=True)
    elif default_pl:
        all_query += ' and default_flag==1'    

    if Teffrange:
        all_query += f' and {Teffrange[0]}<st_teff<{Teffrange[1]}'
    if onlymultiplanet:
        all_query += ' and sy_pnum>1'
    if onlyMplanets:
        all_query += ' and st_teff<3700'
    df_exo.query(all_query,inplace=True)


    if rocky:
        rock = f'(pl_rade-pl_radeerr1)<(pl_masse)**0.27*10**0.03 | '+\
               f'(pl_rade)<(pl_masse+pl_masseerr1)**0.27*10**0.03'
            #    f'(pl_rade-pl_radeerr1*0.5)>(pl_masse+pl_masseerr1*0.5)**0.27*10**0.03'
        df_exo.query(rock,inplace=True)
    if query:
        df_exo.query(query)
    df_exo.reset_index(inplace=True,drop=True)
    return df_exo

def plot_pl(data, color='b', marker='o', Teq=True, show_Teq=True, Teq_param=[4, 0.3, 0], Teq_kwargs=None,
            Mrange=[1, 20], Rrange=[1, 3], axes_yscale="log", axes_xscale="log", show_cmf=True,
            cmf_param=None, show_H2O=True, show_H=False, show_stars=False, stars_param=[2,1,'Fe2Mg'],
            stars_quantiles=None,star_kwargs=None, fig=None, **data_kwargs):
    """
    Plot exoplanet data from NASA archive and add isochemical lines.

    Parameters
    ----------
    data: pandas.DataFrame
        pandas DataFrame object containing the exoplanet data from NASA.
        Note, the table needs to have same column names as NASA.
    color: str, default='b'
        Color of the scatter plot markers.
    marker: str, default='o'
        Marker style of the scatter plot.
    Teq: bool, default=True
        Flag whether to plot Teq.
    show_Teq: bool, default=True
        Flag whether to show the colorbar for Teq.
    Teq_param: list, default=[4, 0.3, 0]
        Parameters for calculating planet Teq.
        The list object should contain [f,A,C]:
            - Scaling factor (f).
            - Fractional albedo (A).
            - Constant/greenhose offset (C).
    Teq_kwargs: dict, default=None
        Additional keyword arguments for customizing Teq scatter plot.
    Mrange: list, default=[1, 20]
        Range of allowed planet masses [min_mass, max_mass].
    Rrange: list, default=[1, 3]
        Range of allowed planet radii [min_radius, max_radius].
    axes_yscale: str, default='log'
        Scale of the y-axis ('linear' or 'log').
    axes_xscale: str, default='log'
        Scale of the x-axis ('linear' or 'log').
    show_cmf: bool, default=True
        Flag whether to plot cmf lines.
    cmf_param: dict, default=None
        Parameters for customizing the cmf plot.
        The dictionary should have the following keys:
            - 'cmf': list, with cmf values to use.
            - 'label': list, with labels matching cmf.
            - 'color': list, with colors matching cmf.
        Example:
        cmf_param = {'cmf': [0., 0.33, 0.67, 1.],
                     'label': ['RTR', 'Earth', 'Mercury', 'Fe'],
                     'color': ['k', 'green', 'red', 'k']}
    show_H2O: bool, default=True
        Flag whether to plot H2O line.
    show_H: bool, default=False
        Flag whether to plot H-He line.
    show_stars: bool, default=False
        Flag whether to plot stellar composition.
    stars_param: list, default=[2, 0.7, 'Fe2Mg']
        Parameters for calculating stellar composition.
        The list object should contain [mu,sigma,'Fe2X']:
            - The mean value of star (mu).
            - The standard deviation of star (sigma).
            - Specify the chemical ratio of the star ('Fe2Mg' or 'Fe2Si').
    stars_quantiles: list, default=None
        List of quantiles to use when ploting stellar composition. 
        If None than use a continuous spectrum.
        Example:
        stars_quantiles = [50,16,2.3] 
    star_kwargs: dict, default=None
        Additional keyword arguments for customizing the stellar composition plot.
    fig: matplotlib.figure.Figure, default=None
        Figure object to be used for the plot.
    **data_kwargs: dict
        Additional keyword arguments for customizing the planet scatter plot.

    Returns
    -------
    matplotlib.figure.Figure
        The created figure object.
    matplotlib.axes._subplots.AxesSubplot
        The created axes object.

    """
    if fig:
        ax = fig.get_axes()[0]
    else:
        fig,ax = plt.subplots(1,1,dpi=160)
    if not Teq_kwargs:
        Teq_kwargs = {}
    if not star_kwargs:
        star_kwargs = {'color':'purple'}
    Mmin,Mmax = Mrange
    Rmin,Rmax = Rrange
    Mass = np.linspace(Mmin,Mmax,30)
    
    if show_stars:
        FeX = np.random.normal(stars_param[0],stars_param[1],1000)
        cmf_st = st_pl(FeX,method=stars_param[2])
        if not stars_quantiles:
            X,Y,Z = plot_cont(cmf_st,Mrange)
            ax.contourf(X,Y,Z,len(X),zorder=0,**star_kwargs)    
        else:
            for tmp,percentile in enumerate(stars_quantiles):
                cmf_perc = np.percentile(cmf_st,[percentile,100-percentile])
                R1 = [guess_R(mass,0,cmf=cmf_perc[0]) for mass in Mass]
                R2 = [guess_R(mass,0,cmf=cmf_perc[1]) for mass in Mass]
                if R1==R2:
                    ax.plot(Mass,R1,alpha=1,zorder=0,**star_kwargs)
                else:
                    ax.fill_between(Mass,R1,R2,alpha=0.9/(tmp+1),zorder=0,**star_kwargs)
    if show_cmf:
        if not cmf_param:
            cmf_param = {'cmf':[0.,0.33,0.67,1.],
                         'label':['RTR','Earth','Merc','Fe'],
                         'color':['k','green','red','k']}
        for cmf,label,c in zip(cmf_param['cmf'], cmf_param['label'], cmf_param['color']):
            plot_cmf(ax,Mass,cmf,label,c)
    if show_H2O:
        #plot H2O envelope planet
        M,R = MR_H2O()
        ax.plot(M,R,color='C0')
    if show_H:
        #plot H-He envelope planet
        M,R = MR_HHe()
        ax.plot(M,R,color='blue')

    #plot the data
    ax.errorbar(data.pl_masse,data.pl_rade,xerr=[-data.pl_masseerr2,data.pl_masseerr1],
                yerr=[-data.pl_radeerr2,data.pl_radeerr1],fmt=' ',color='k',alpha=0.6,zorder=1)
    
    if Teq:
        #some magic numbers
        Rsun = 6.95508e8 #in m
        Au = 1.496e11 #in m
        a = ((data.pl_orbper/365.25)**2*data.st_mass)**(1/3) #in Au
        Flux = (data.st_teff)**4*(data.st_rad*Rsun)**2/(Au*a)**2 #no constants
        f,A,C = Teq_param
        Teq = C+(Flux/f)**0.25*(1-A)**0.25 
        
        cm = ax.scatter(data.pl_masse,data.pl_rade,c=Teq,zorder=1,**Teq_kwargs)
        if show_Teq:
            cb = fig.colorbar(cm,label=r'$T_{eq}$ (K)')
    else:
        ax.scatter(data.pl_masse,data.pl_rade,color=color,marker=marker,zorder=1,**data_kwargs)
    ax.set_yscale(axes_yscale)
    ax.set_xscale(axes_xscale)
    ax.set_xlim(Mmin,Mmax)
    ax.set_ylim(Rmin,Rmax)
    ax.set_xlabel(r'M/M$_\oplus$')
    ax.set_ylabel(r'R/R$_\oplus$')
    formatter = FuncFormatter(lambda y, _: '{:.16g}'.format(y))
    if axes_yscale=='log':
        yticks = np.concatenate([(np.linspace(Rmin,2,4)).round(1),
                                 (np.linspace(2,Rmax,4)).round(1)])
        ax.set_yticks(yticks)
        ax.yaxis.set_major_formatter(formatter)
    if axes_xscale=='log':
        xticks = (10**np.linspace(np.log10(Mmin),np.log10(Mmax),8)).round(0)
        ax.set_xticks(xticks)
        ax.xaxis.set_major_formatter(formatter)

    return fig,ax


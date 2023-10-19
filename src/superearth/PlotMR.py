import os
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter,NullFormatter
import numpy as np
import pandas as pd
from superearth import st_pl,guess_R
from superearth.utils import check_cache_archive,plot_cont,plot_cmf,MR_H2O,MR_HHe

package_dir = os.path.dirname(__file__)
check_cache_archive()
def exoplanets(Merr, Rerr, default_pl=True, best_pl=False, Mrange=[0, 30], Rrange=[0, 4.5], Teffrange=None,
               onlymultiplanet=False, onlyMplanets=False, rocky=False, limit=None, exclude=None,query=None):
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
        Custom query string to filter the data, replaces error requests.
    limit: list, default: None
        Query only for the planets in the list.
    exclude: list, default: None
        Exclude planets in the list from query.
    Returns
    -------
    class: pandas.DataFrame
        pandas DataFrame object containing the filtered exoplanet data.

    """
    df_exo = pd.read_csv(package_dir+"/Data/ExoplanetArchiveData.csv")    
    if query:
        all_query = query
    else:
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
    if exclude:
        df_exo.query('pl_name!=@exclude',inplace=True)
    if limit:
        df_exo.query('pl_name==@limit',inplace=True)
    df_exo.reset_index(inplace=True,drop=True)
    return df_exo

def plot_pl(data, color='b', marker='o', Teq=True, show_Teq=True, Teq_param=[4, 0.3, 0], Teq_kwargs=None,
            Mrange=None, Rrange=None, axes_yscale="log", axes_xscale="log", show_cmf=True,
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
    Mmin,Mmax = min(data.pl_masse)*0.9,max(data.pl_masse)*1.1
    Rmin,Rmax = min(data.pl_rade)*0.9,max(data.pl_rade)*1.1
    if Mrange:
        Mmin,Mmax = Mrange
    if Rrange:
        Rmin,Rmax = Rrange
    Mass = 10**np.linspace(np.log10(Mmin),np.log10(min([Mmax,50])),128)
    
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
        if data.st_rad[0]>0:
            #some magic numbers
            Rsun = 6.95508e8 #in m
            Au = 1.496e11 #in m
            a = ((data.pl_orbper/365.25)**2*data.st_mass)**(1/3) #in Au
            Flux = (data.st_teff)**4*(data.st_rad*Rsun)**2/(Au*a)**2 #no constants
            f,A,C = Teq_param
            Teq_pl = C+(Flux/f)**0.25*(1-A)**0.25 
        else:
            Teq_pl = data.pl_eqt
        cm = ax.scatter(data.pl_masse,data.pl_rade,c=Teq_pl,zorder=1,**Teq_kwargs)
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
    formatter = FuncFormatter(lambda y, _: '{:.15g}'.format(y))
    if axes_yscale=='log':
        # yticks = np.concatenate([(np.linspace(Rmin,2,4)).round(1),
        #                          (np.linspace(2,Rmax,4)).round(1)])
        yticks = np.concatenate([(np.linspace(Rmin,Rmax,8)).round(1)])
        ax.set_yticks(yticks)
        ax.yaxis.set_major_formatter(formatter)
        ax.yaxis.set_minor_formatter(NullFormatter())
    if axes_xscale=='log':
        xticks = (10**np.linspace(np.log10(Mmin),np.log10(Mmax),8)).round(0)
        xticks[0] = np.max([xticks[0],Mmin]).round(1)
        ax.set_xticks(xticks)
        ax.xaxis.set_major_formatter(formatter)

    return fig,ax


def plotly_pl(data, color='black', marker='circle',size=10, Teq=True, Teq_param=[4, 0.3, 0], Teq_kwargs=None,
            Mrange=None, Rrange=None, axes_yscale="log", axes_xscale="log", show_cmf=True,
            cmf_param=None, show_H2O=True, show_H=False, show_stars=False, stars_param=[2, 1, 'Fe2Mg'],
            stars_quantiles=[2.3,16,50], star_kwargs=None, fig=None):
    """
    Plotly version of the plotting exoplanet data.
    For plotly to work in Jupyterlab, need to have ipywidgets and jupyter-dash modules.
    Cehck plotly getting-started section for more information.
    """
    try:
        import plotly.graph_objects as go
    except:
        print("Importing plotly failed, make sure you have all the modules installed.")
        return None
    
    if fig:
        fig = go.FigureWidget(fig)
        ax = fig.data[0]
    else:
        fig = go.FigureWidget()
        ax = go.Scatter()

    if not Teq_kwargs:
        Teq_kwargs = {}
    if not star_kwargs:
        star_kwargs = {'color': 'indigo'}
    Mmin,Mmax = min(data.pl_masse)*0.9,max(data.pl_masse)*1.1
    Rmin,Rmax = min(data.pl_rade)*0.9,max(data.pl_rade)*1.1
    if Mrange:
        Mmin,Mmax = Mrange
    if Rrange:
        Rmin,Rmax = Rrange
    Mass = 10**np.linspace(np.log10(Mmin),np.log10(min([Mmax,50])),128)

    if show_stars:
        FeX = np.random.normal(stars_param[0], stars_param[1], 1000)
        cmf_st = st_pl(FeX, method=stars_param[2])
        for tmp, percentile in enumerate(stars_quantiles):
            cmf_perc = np.percentile(cmf_st, [percentile, 100 - percentile])
            R1 = [guess_R(mass, 0, cmf=cmf_perc[0]) for mass in Mass]
            R2 = [guess_R(mass, 0, cmf=cmf_perc[1]) for mass in Mass]

            fig.add_trace(go.Scatter(x=Mass, y=R1, fill=None, mode='lines', line_color='indigo', 
                                     showlegend=False, hoverinfo='skip'))
            if R1 != R2:
                # Add R2 curve with fill
                fig.add_trace(go.Scatter(x=Mass, y=R2, fill='tonexty', mode='lines', line_color=star_kwargs['color'],
                                         showlegend=False, hoverinfo='skip'))

    if show_cmf:
        if not cmf_param:
            cmf_param = {'cmf': [0., 0.33, 0.67, 1.],
                         'label': ['RTR', 'Earth', 'Merc', 'Fe'],
                         'color': ['black', 'green', 'red', 'black']}
        for cmf, label, color in zip(cmf_param['cmf'], cmf_param['label'], cmf_param['color']):
            Radius = np.empty(len(Mass))
            pos = (Mass>1)  # get the Masses that can be evaluted by the function
            Radius[pos] = np.array([guess_R(mass,0,cmf=cmf) for mass in Mass[pos]]) 
            ii = np.argmax(pos)
            m,c = np.polyfit(np.log10(Mass[ii:ii+5]),np.log10(Radius[ii:ii+5]),1) # get the best-fit parameters for straight line
            Radius[np.invert(pos)] = 10**(m*np.log10(Mass[np.invert(pos)])+c) # get the Radius that is outside the function 
            fig.add_trace(go.Scatter(x=Mass, y=Radius, mode='lines', line=dict(color=color), 
                                     showlegend=False, hoverinfo='skip'))
            # fig.add_annotation(x=np.log10(Mass[-1]),y=np.log10(Radius[-1]),
            #     text=label,font=dict(size=15, color=color),showarrow=False,
            #     xanchor='center',yanchor='bottom')
    if show_stars:
        fig.add_trace(go.Scatter(x=[Mass[-1]], y=[R1[-1]], showlegend=False, 
                    mode='markers', marker=dict(symbol='star', size=15, color='white'),       
                    hovertemplate='%{customdata[2]}:   %{customdata[0]} +/- %{customdata[1]}',
                    customdata=[(stars_param[0], stars_param[1], stars_param[2])]))

    if show_H2O:
        # plot H2O envelope planet
        M, R = MR_H2O()
        fig.add_trace(go.Scatter(x=M, y=R, mode='lines', line=dict(color='blue'),showlegend=False,hoverinfo='skip'))

    if show_H:
        # plot H-He envelope planet
        M, R = MR_HHe()
        fig.add_trace(go.Scatter(x=M, y=R, mode='lines', line=dict(color='blue'),showlegend=False,hoverinfo='skip'))

    if Teq:
        # some magic numbers
        Rsun = 6.95508e8  # in m
        Au = 1.496e11  # in m
        a = ((data.pl_orbper / 365.25) ** 2 * data.st_mass) ** (1 / 3)  # in Au
        Flux = (data.st_teff) ** 4 * (data.st_rad * Rsun) ** 2 / (Au * a) ** 2  # no constants
        f, A, C = Teq_param
        Teq_pl = C + (Flux / f) ** 0.25 * (1 - A) ** 0.25
        color = Teq_pl
        Teq_pl = np.round(Teq_pl, decimals=0)
        marker_dict = dict(color=color, colorscale="Magma",symbol=marker, opacity=1, size=size, 
                                colorbar=dict(thickness=10))
    else:
        Teq_pl = data.pl_eqt
        marker_dict = dict(color=color, symbol=marker, opacity=1, size=size)


    # plot the data
    error_x = dict(type='data', array=data.pl_masseerr1, arrayminus=-data.pl_masseerr2, visible=True)
    error_y = dict(type='data', array=data.pl_radeerr1, arrayminus=-data.pl_radeerr2, visible=True)
    scatter_trace = go.Scatter(
        x=data.pl_masse,
        y=data.pl_rade,
        mode='markers',
        showlegend=False,
        error_x=error_x,
        error_y=error_y,
        marker=marker_dict,
        hovertemplate='Name: %{text}<br>Mass: %{x}<br>Radius:%{y}<br>Teq: %{customdata[0]}<br>Reference: %{customdata[1]}',
        text=data.pl_name,
        customdata=list(zip(Teq_pl, data.pl_refname))
    )

    # Add the scatter plot trace to the figure
    fig.add_trace(scatter_trace)

    fig.update_layout(width=1000, height=600,
        xaxis=dict(type=axes_xscale, range=[np.log10(Mmin), np.log10(Mmax)]),
        yaxis=dict(type=axes_yscale, range=[np.log10(Rmin), np.log10(Rmax)]),
        xaxis_title=r'M',
        yaxis_title=r'R'
    )
    return fig
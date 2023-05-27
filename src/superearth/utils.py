import os
from datetime import datetime, timedelta
import pandas as pd
import numpy as np
from superearth import guess_R
import warnings


package_dir = os.path.dirname(__file__)
def check_cache_archive():
    #Check if the cached table file exists
    cache_file = package_dir+"/Data/ExoplanetArchiveData.csv"
    if not os.path.isfile(cache_file):
        update()
    else:
        #Check when the file was last modified
        modified_time = datetime.fromtimestamp(os.path.getmtime(cache_file))
        if datetime.now() - modified_time > timedelta(days=30): 
            message = "\nThe NASA archive data was last updated more than a month ago.\n" \
                      "If you wish to update, run se.utils.update()"
            warnings.warn(message)
def update():
    #quering ipac nasa database
    listdb = """default_flag,pl_name,hostname,sy_pnum,discoverymethod,disc_year,disc_facility,
                pl_refname,pl_orbper,pl_orbpererr1,pl_orbpererr2,pl_orbperlim,
                pl_rade,pl_radeerr1,pl_radeerr2,pl_radelim,pl_masse,pl_masseerr1,pl_masseerr2,
                pl_masselim,pl_orbeccen,pl_orbeccenerr1,pl_orbeccenerr2,pl_orbeccenlim,
                ttv_flag,st_refname,st_spectype,st_teff,st_tefferr1,st_tefferr2,
                st_tefflim,st_rad,st_raderr1,st_raderr2,st_radlim,st_mass,st_masserr1,
                st_masserr2,st_masslim,rowupdate"""
    listdb = listdb.replace("\n", "").replace(" ", "") 
    nasa_url = "https://exoplanetarchive.ipac.caltech.edu/TAP/sync?query=select" 
    query = "where+default_flag+between+0+and+1" #get all published data for planets
    url = f"{nasa_url}+{listdb}+from+ps+{query}&format=csv"
    df_orig = pd.read_csv(url)
    df_orig.sort_values(by=['pl_name', 'rowupdate'],inplace=True)
    df_orig.reset_index(inplace=True,drop=True)                        
    df_orig.to_csv(package_dir+"/Data/ExoplanetArchiveData.csv") # save updated table
    return 
def plot_cont(cmf,Mrange=[1,20]):
    pdf,bins = np.histogram(cmf,int(np.sqrt(len(cmf))),
                        range=(0,max(cmf)),density=True)
    cmf_x = (bins[:-1]+bins[1:])/2
    N = len(cmf_x)
    X,Y,Z = np.zeros([3,N,N])
    M = np.linspace(Mrange[0],Mrange[1],N) 
    for i in range(N):
        R = [guess_R(mass,0,cmf=cmf_x[i]) for mass in M]
        X[i] = M
        Y[i] = R
        Z[i] += pdf[i] #pdf is same along the M-R line
    return X,Y,Z    
def plot_cmf(ax,Mass,cmf,label,color):
    Radius = np.array([guess_R(mass,0,cmf=cmf) for mass in Mass])
    ax.plot(Mass,Radius,color=color,zorder=0)
    dx,dy = Mass[1]-Mass[0],Radius[2]-Radius[0]
    deg = np.arctan(dy/dx)/2/np.pi*360
    ax.text(Mass[-1]-len(label)*0.5-3,Radius[-1]*0.97,label,rotation=deg,fontsize=8,color=color,
        horizontalalignment='left',verticalalignment='bottom')
    return ax
def MR_H2O():
    #plot H2O envelope planet
    M,R = np.loadtxt(package_dir+'/Data/MR_H2O.txt')
    return M,R
def MR_HHe():
    #plot H-He envelope planet
    M,R = np.loadtxt(package_dir+'/Data/MR_H1.txt')
    return M,R

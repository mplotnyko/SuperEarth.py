import matplotlib.pyplot as plt
import superearth as se

#get only rocky planets
df = se.exoplanets(0.25,0.08,rocky=True)
fig,ax = se.plot_pl(df,color='grey',Teq=False,alpha=0.9,marker='.') #plot M-R data

#get onyl planets around Mdwarfs
df = se.exoplanets(0.25,0.08,onlyMplanets=True)
fig,ax = se.plot_pl(df,fig=fig) #plot M-R data, on the same fig
plt.savefig('MR_rockyMpl.jpg',bbox_inches='tight')

#get only rocky planets sorted by the latest upload data to NASA
df = se.exoplanets(0.25,0.08,best_pl=True,rocky=True)
fig,ax = se.plot_pl(df,color='grey',Teq=False,alpha=0.9,marker='.',
                    show_stars=True,stars_quantiles=[50,16,2.3],show_cmf=True,show_H2O=False)
plt.savefig('MR_rockybest.jpg',bbox_inches='tight')

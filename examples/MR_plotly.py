import superearth as se

df = se.exoplanets(0.25,0.08,rocky=False)
fig = se.plotly_pl(df,color='black', marker='circle') #plot M-R data using plotly

# Display the figure
fig.show()
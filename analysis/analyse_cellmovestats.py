# conda install -c conda-forge trackpy
# conda install -c conda-forge pims
import pandas as pd
import trackpy as tp
import numpy as np
from matplotlib import pyplot as plt

df=pd.read_csv('/home/anton/chirality/program/prod_v9_ar1.3_L20.0_layer13.0_desync0.4_seed2.0_b4e3_twist1.0_adv.dat',delimiter='\s+',header=None,
               names=['x','y','th(i)','d(i)','d(i)*alpha(i)','depth(i)','c(i)','particle','frame','fx','fy','vx','vy'])
dfgr=df.groupby('particle')

for cell,data in dfgr:
    fig, ax = plt.subplots(1,2)
    #print(cell)
    #if cell==10:
    #    data.plot.scatter(x='frame',y='th(i)')
    #    plt.figure()
    #    data.plot.scatter(x='x',y='y')
    #data.plot.scatter(x='frame',y='th(i)')
    ax[0].scatter(data['x'],data['y'],c=data['th(i)'])
    ax[1].scatter(data['frame'],data['th(i)'],c=data['th(i)'])
    plt.show()

#tp.plot_traj(df)
#dfuf=df[(df['y']>0) & (df['depth(i)']<1)] # upper front
#dfgr=dfuf.groupby('particle')
#for cell,data in dfgr:
#        data.plot.scatter(x='frame',y='th(i)')
#tp.plot_traj(dfuf)

#d=tp.compute_drift(dfuf)

#d.plot()

#d.plot(x='x',y='y')

# http://soft-matter.github.io/trackpy/v0.4.2/tutorial/walkthrough.html
#for part in [4,5,6,7]:
#  xs=df[df['particle']==part]['x'].values
#  ys=df[df['particle']==part]['y'].values
#  fx=df[df['particle']==part]['fx'].values
#  fy=df[df['particle']==part]['fy'].values
#  ths=df[df['particle']==part]['th(i)'].values
#  thsmov=np.arctan(np.gradient(ys)/np.gradient(xs))
#  #plt.plot(ths-thsmov,'.')
#  vxg=np.gradient(xs)
#  vyg=np.gradient(ys)
#  vx=df[df['particle']==part]['vx'].values
#  vy=df[df['particle']==part]['vy'].values
  #plt.plot(np.arccos((vx*fx+vy*fy)/np.sqrt(vx*vx+vy*vy)/np.sqrt(fx*fx+fy*fy)))
#  b=4e3
#  bani=-4e3
  #plt.plot(vx-b*fx,bani*(fx*np.cos(2*ths)+fy*np.sin(2*ths)),'.')
  #plt.plot(vy-b*fy,bani*(fx*np.sin(2*ths)-fy*np.cos(2*ths)),'.')
#  plt.plot(vx,vxg)
#  plt.plot(vy,vyg)


#dfuf.plot(x='x',y='th(i)')


# conda install -c conda-forge trackpy
# conda install -c conda-forge pims
import pandas as pd
import trackpy as tp
import numpy as np
from matplotlib import pyplot as plt

df=pd.read_csv('myfile.dat',delimiter='\s+',header=None,
               names=['x','y','th(i)','d(i)','d(i)*alpha(i)','depth(i)','c(i)','particle','frame','vx(i)','vy(i)','vth(i)'])

df=df[df['y']>0]

dfgr=df.groupby('particle')

#plt.plot(df['x'].values[-10000:],df['y'].values[-10000:],'.')
#plt.show()

cellrot_front=np.array([])
cellrotvel_front=np.array([])
cellvx_front=np.array([])
cellvy_front=np.array([])
cellrot_bulk=np.array([])
cellrotvel_bulk=np.array([])
cellvx_bulk=np.array([])
cellvy_bulk=np.array([])

fig, ax = plt.subplots(2,4)
for cell,data in dfgr:
    #fig, ax = plt.subplots(1,2)
    #print(cell)
    #if cell==10:
    #    data.plot.scatter(x='frame',y='th(i)')
    #    plt.figure()
    #    data.plot.scatter(x='x',y='y')
    #data.plot.scatter(x='frame',y='th(i)')
    #ax[0].scatter(data['x'],data['y'],c=data['th(i)'])
    #ax[1].scatter(data['frame'],data['th(i)'],c=data['th(i)'])
    #plt.show()
    cellrot_front=np.append(cellrot_front,data['th(i)'].values[data['depth(i)'].values<0.1])
    cellrotvel_front=np.append(cellrotvel_front,data['vth(i)'].values[data['depth(i)'].values<0.1])
    cellvx_front=np.append(cellvx_front,data['vx(i)'].values[data['depth(i)'].values<0.1])
    cellvy_front=np.append(cellvy_front,data['vy(i)'].values[data['depth(i)'].values<0.1])
    cellrot_bulk=np.append(cellrot_bulk,data['th(i)'].values[data['depth(i)'].values>5.1])
    cellrotvel_bulk=np.append(cellrotvel_bulk,data['vth(i)'].values[data['depth(i)'].values>5.1])
    cellvx_bulk=np.append(cellvx_bulk,data['vx(i)'].values[data['depth(i)'].values>5.1])
    cellvy_bulk=np.append(cellvy_bulk,data['vy(i)'].values[data['depth(i)'].values>5.1])

ax[0,0].hist(np.fmod(cellrot_front+2*np.pi,2*np.pi),bins=1000)
ax[0,1].hist(cellrotvel_front,bins=10000)
ax[0,1].set_yscale('log')
ax[0,1].set_xlim([-2.5,2.5])
ax[0,2].hist(cellvx_front,bins=10000)
ax[0,2].set_yscale('log')
ax[0,2].set_xlim([-2.5,2.5])
ax[0,3].hist(cellvy_front,bins=10000)
ax[0,3].set_yscale('log')
ax[0,3].set_xlim([-2.5,2.5])
ax[1,0].hist(np.fmod(cellrot_bulk+2*np.pi,2*np.pi),bins=1000)
ax[1,1].hist(cellrotvel_bulk,bins=10000)
ax[1,1].set_yscale('log')
ax[1,1].set_xlim([-2.5,2.5])
ax[1,2].hist(cellvx_bulk,bins=10000)
ax[1,2].set_yscale('log')
ax[1,2].set_xlim([-2.5,2.5])
ax[1,3].hist(cellvy_bulk,bins=10000)
ax[1,3].set_yscale('log')
ax[1,3].set_xlim([-2.5,2.5])
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

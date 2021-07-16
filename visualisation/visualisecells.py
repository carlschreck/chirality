import pandas as pd
from matplotlib.patches import Ellipse
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
import matplotlib.cm as cm
from ipywidgets import interactive, IntSlider, fixed, interact
import sys

df=pd.read_csv(sys.argv[1],delimiter='\s+',header=None,
               names=['x','y','th(i)','d(i)','d(i)*alpha(i)','depth(i)','c(i)','particle','frame','vx(i)','vy(i)','vth(i)'])


df=df[df['y']>=0]

def showcells(df, frame = 10):
    #norm = matplotlib.colors.Normalize(vmin=-np.max(np.abs(df['vth(i)'])), vmax=np.max(np.abs(df['vth(i)'])), clip=True)
    df=df[df['frame']==frame]
    #mapper = cm.ScalarMappable(norm=norm, cmap=cm.seismic)
    fig, ax = plt.subplots(subplot_kw={'aspect': 'equal'}, figsize = (10, 4))
    for index, row in df.iterrows():
        e=Ellipse(xy=[row['x'],row['y']],
                  height=row['d(i)'],
                  width=row['d(i)*alpha(i)'],
                  angle=row['th(i)']/2/np.pi*360)
        ax.add_artist(e)
        #e.set_facecolor(mapper.to_rgba(row['vx(i)']))
    #ax.set_xlim(np.min(df['x'])-10,np.max(df['x'])+10)
    #ax.set_ylim(np.min(df['y'])-10,np.max(df['y'])+10)
    ax.set_xlim(np.mean(df['x'])-100,np.mean(df['x'])+100)
    ax.set_ylim(np.mean(df['y'])-50,np.mean(df['y'])+50)
    ax.set_aspect('equal', 'box')
    fig.tight_layout()
    plt.savefig('myfile{:04d}.png'.format(frame))
    #plt.show()
    plt.close()

for i in range(1,np.max(df['frame'])+1):
    print(i)
    showcells(df, frame=i)

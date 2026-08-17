#!/usr/bin/env python

import sys
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

observed_bc = sys.argv[1]
read_threshold = int(sys.argv[2])

df = pd.read_csv(observed_bc, sep='\s+', names=['bc_count','bc_name'])
print("length DF",str(len(df))) ### gremlin

print(df['bc_name'].value_counts()) ### gremlin

df1 = df[~df['bc_name'].str.contains('=')]
print("length DF1",str(len(df1))) ### gremlin


cter = df1.set_index('bc_name')['bc_count'].sort_values(ascending=False) #ascending=False to get largest at the top

print("length cter:",str(len(cter))) ### gremlin
xmin = min(cter)
xmax = max(cter)*10.


def startfig(w=4,h=2,rows=1,columns=1,wrs=None,hrs=None,frameon=True,return_first_ax=True):

    '''
    for initiating figures, w and h in centimeters
    example of use:
    a,fig,gs = startfig(w=10,h=2.2,rows=1,columns=3,wr=[4,50,1],hrs=None,frameon=True)
    hrs - height ratios
    wrs - width ratios
    frameon - whether first axes with frame
    
    returns:
    if return_first_ax=True
    a,fig,gs
    else
    fig,gs
    '''
    
    ratio = 0.393701 #1 cm in inch
    myfigsize = (w*ratio,h*ratio)
    fig = plt.figure(figsize = (myfigsize))
    gs = mpl.gridspec.GridSpec(rows, columns ,width_ratios=wrs,height_ratios=hrs)
    if return_first_ax==True:
        a = fig.add_subplot(gs[0,0],frameon=frameon)
        return a,fig,gs
    else:
        return fig,gs


a,fig,gs = startfig(7,5)

# histogram
bins=np.logspace(np.log10(xmin),np.log10(xmax),51)
hs, bins,patches = plt.hist(cter,bins=bins)

#plot barchart
lefts = bins[:-1]
rights = bins[1:]
a.bar(x = lefts,width = rights-lefts,height = hs*rights,
      align='edge',
      lw=0.,color = 'c')
a.set_xscale('log')
a.set_xlim(xmin,xmax)

# threshold:
a.plot((read_threshold,read_threshold),(a.get_ylim()[0],a.set_ylim()[1]),lw=1,color='r')
a.set_xlabel('Reads per barcode')

# print pct reads
tot = cter.sum()
hi = (cter[cter>=read_threshold]).sum()/tot*100.
lo = (cter[cter<read_threshold]).sum()/tot*100.

pad = 0.1
xpad = (a.get_xlim()[1]-a.get_xlim()[0])*pad
ypad = (a.get_ylim()[1]-a.get_ylim()[0])*pad

a.text(0.1,0.9,"%.1f%%\nof reads"%lo,va='center',ha='center',transform=a.transAxes)
a.text(0.9,0.9,"%.1f%%\nof reads"%hi,va='center',ha='center',transform=a.transAxes)

a.set_ylabel('# reads from bin')
plt.savefig("4_observed_bc_dist.png")
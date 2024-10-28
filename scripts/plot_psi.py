#!/opt/python/bin/python3
from sys import argv
import pandas as pd
import numpy as np
from scipy.stats import gaussian_kde,pearsonr
import matplotlib.pyplot as plt

simfile = argv[1]
nmdjdir = argv[2]
outfile = argv[3]

catfile = f'{nmdjdir}/classification.tsv'
psifile = f'{nmdjdir}/PSI.tsv'


sim = pd.read_csv(simfile,sep='\t')
cat = pd.read_csv(catfile,sep='\t')
psi = pd.read_csv(psifile,sep='\t')

cat = cat[cat.nmd_id_list.notna()&cat.coding_id_list.notna()]

events = []
for _,row in cat.iterrows():
    for tr in row.nmd_id_list.split(', '):
        events.append([row.event_id,'nmd',tr])
    for tr in row.coding_id_list.split(', '):
        events.append([row.event_id,'coding',tr])
cat = pd.DataFrame(events,columns=['event_id','biotype','transcript_id'])

sim = cat.merge(sim).groupby(['event_id',
                              'biotype']).TPM.agg(np.nansum).reset_index().pivot(index='event_id',
                                                                                 columns='biotype',
                                                                                 values='TPM').reset_index()
sim['ground truth $\Psi$'] = sim.nmd/(sim.nmd+sim.coding)
sim.drop(['nmd','coding'],axis=1,inplace=True)
m = psi[psi.denom>30].rename(columns={'psi': 'NMDj $\Psi$'}).merge(sim)

def MSE(x,y):
    return np.nanmean((x-y)**2)
def corcoef(m,col1,col2):
    c = m[(m[col1].notna())&(m[col2].notna())]
    return pearsonr(c[col1],c[col2])[0]

def plot_density(temp,ax,col1='ground truth $\Psi$',col2='NMDj $\Psi$',log=True,**kwargs):
    x = temp[(temp[col1].notna())&(temp[col2].notna())][col1].to_numpy()
    y = temp[(temp[col1].notna())&(temp[col2].notna())][col2].to_numpy()
    xy = np.vstack([x,y])
    z = gaussian_kde(xy)(xy)
    if log:
        z = np.log2(z)
    idx = z.argsort()
    x, y, z = x[idx], y[idx], z[idx]
    ax.scatter(x,y,c=z,**kwargs)
    ax.set_xlabel(col1)
    ax.set_ylabel(col2)
    ax.plot([0,1],[0,1],lw=1,ls=':',c='black')
    mse = MSE(temp[col1],temp[col2])
    cor = corcoef(m,col1,col2)
    ax.set_title(f"MSE = {round(mse,3)}, Pearson corcoef = {round(cor,2)}")

plt.rcParams.update({'font.size': 14})
f,ax = plt.subplots()
plot_density(m,ax)
f.savefig(outfile,transparent=False,facecolor='white',dpi=300)
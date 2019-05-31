'''
to run : 
  python -i basecode.py
simulation from McDougal RA, Hines ML, Lytton WW. Reaction-diffusion in the NEURON simulator. 
Front Neuroinformatics 2013 7:28 
http://www.frontiersin.org/neuroinformatics/10.3389/fninf.2013.00028/abstract
http://neurosimlab.org/pdfs/fninf-07-00028.pdf
'''

from neuron import h, rxd, gui
import numpy as np 
import pylab as plt
import gzip
import pickle as pk

# set parameters
dend = h.Section()
dend.L = 1000
dend.nseg = 200 
alpha = 0.3

# Define where?: whole dend;  what?: c stuff; how?: cprime reaction
all_secs = rxd.Region(h.allsec(), nrn_region='i')
c = rxd.Species(all_secs, d=10, name='c')
reaction = rxd.Rate(c, (0 - c) * (alpha - c) * (1 - c)) # a positive feedback system: c'>0 if c>alpha

def initial ():
  'set up initial conditions'
  h.finitialize()
  for nd in c.nodes:
    nd.concentration = 1 if nd.x < 0.2 else 0 # start with all the 

def cprime (c):
  "for setting up c v c' plot"
  return (0-c)*(alpha-c)*(1-c)

def show (): 
  'output sections of nodes'
  from  pprint import pprint as pp
  pp([(x,format(c.nodes.concentration[x-10:x+10])) for x in range(100,1000,100)])

def Vtvec(vlen, Dt):
  return np.linspace(0, vlen*Dt, vlen, endpoint=True)

# creating vectors to record across multiple nodes
vecdict = {}
def recnodes (steps=5):
  'set up vectors for recording from multiple nodes'
  vecdict.clear()
  noderange = [x-1 for x in range(0,dend.nseg+1,dend.nseg/steps)]
  noderange[0]=0
  for loc in noderange:
    vec = h.Vector(1000/h.dt)
    vecdict.update({'node# %d'%(loc):vec})
    vec.record(c.nodes[loc]._ref_concentration)

lendict = {}
def runtimes ():
  'create a list of concentration values at different times throughout the simulation'
  initial()
  lendict.clear()
  lendict.update({'time = 0':c.nodes.concentration})
  for t in [200, 400, 600, 800, 1000]:
    h.continuerun(t)
    lendict.update({'time = %d'%(t):c.nodes.concentration})

figlist=[]
def plotgraph (dct, show = None, xvals = None, xlabel = 'xaxis', newfig=False):
  'Create a graph of concentrations for a data set'
  global figlist
  if show is None: show = range(len(dct))
  if xvals is None: xvals = range(len(dct.values()[0]))
  if newfig or len(figlist)==0:
    figlist.append(plt.figure(figsize=(6,4), tight_layout=True))
  else:
    plt.figure(0)
    plt.clf()
  plt.xlabel(xlabel); plt.ylabel('[c]')
  for i,x in enumerate(sorted(dct.keys(),key=lambda(x): int(x.split()[-1]))): # sort on labels numerically
    if i in show: 
      plt.plot(xvals, dct[x], label=x)  
  plt.legend()
  plt.xlim(0, xvals.max())
  plt.ylim(0, 1.1)
  
def savedata(myfile = 'simdata'):
  datalist = [vecdict, lendict]
  savelist = pk.dumps(datalist)
  datazip = gzip.open(myfile, 'w')
  datazip.write(savelist)
  datazip.close() 

def loaddata(fname = 'simdata'):
  global simlists
  rdzipdata = gzip.open(fname, 'r')
  rdsimdata = rdzipdata.read()
  simlists = pk.loads(rdsimdata)
  return simlists 

def runseq(plotflag=True, saveflag=False):
  'Running simulation and plot'
  recnodes()
  runtimes()
  if plotflag:
    plotgraph(vecdict, xvals = Vtvec(len(vecdict.values()[0]), h.dt), xlabel='time (ms)', newfig=True)
    plotgraph(lendict, xvals = np.linspace(0, dend.L, dend.nseg), xlabel='distance (microns)', newfig=True)
  if saveflag:
    savedata()

if __name__ == '__main__':
  runseq()
  plt.show()

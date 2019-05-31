from neuron import h, rxd

alpha = 0.3

# define the geometry
dend = h.Section()
dend.L = 1000
dend.nseg = 1000

# define the reaction-diffusion
all_secs = rxd.Region(h.allsec(), nrn_region='i')
c = rxd.Species(all_secs, d=1, name='c')
reaction = rxd.Rate(c, (0 - c) * (alpha - c) * (1 - c))

# initialize
h.finitialize()
for seg in c.nodes:
    seg.concentration = 1 if seg.x < 0.2 else 0

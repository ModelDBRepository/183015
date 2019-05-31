from neuron import h, rxd
import numpy
import os
import sys

h.load_file('stdrun.hoc')

rxd.options.use_reaction_contribution_to_jacobian = False
sqrt2 = numpy.sqrt(2)

# record times (must be integers)
t_lo = 200
t_hi = 600

length = 1000

nsegs = [250, 500, 1000, 2000]
alphas = list(numpy.arange(0.01, 0.5, 0.05))

all_errors = []
for nseg in nsegs:
    errors = []
    for alpha in alphas:
        print 'with alpha = %g' % alpha
        if not os.fork():
            # define the geometry
            dend = h.Section()
            dend.L = length
            dend.nseg = nseg

            # define the reaction-diffusion
            all_secs = rxd.Region(h.allsec(), nrn_region='i')
            c = rxd.Species(all_secs, d=1, name='c')
            reaction = rxd.Rate(c, (0 - c) * (alpha - c) * (1 - c))

            # initialize
            h.finitialize()
            for seg in c.nodes:
                seg.concentration = 1 if seg.x < 0.2 else 0

            h.CVode().active(1)
            h.CVode().atol(1e-13)

            positions = {}

            def record_position():
                positions[int(h.t)] = waveposition()

            def waveposition():
                x = numpy.linspace(0, dend.L, dend.nseg + 1)
                x = x[:-1] + dend.L / (2. * dend.nseg)
                y = c.states
                for i, value in enumerate(y):
                    if value < 0.3:
                        break
                else:
                    raise Exception('entire domain above threshold')
                
                print 'y(%g) = %g' % (x[i - 1], y[i - 1])
                print 'y(%g) = %g' % (x[i], y[i])
                # y - y1 = m * (x - x1), so x = (y - y1) / m + x1
                x1, x2 = x[i - 1 : i + 1]
                y1, y2 = y[i - 1 : i + 1]
                m = (y2 - y1) / (x2 - x1)
                x = (alpha - y1) / m + x1
                print 'y(%g) = %g (estimated)' % (x, alpha)
                return x
                

            for time in [t_lo, t_hi]:
                h.CVode().event(time, record_position)

            h.continuerun(t_hi + 1)

            print 'distance traveled: %g' % (positions[t_hi] - positions[t_lo])
            speed = (positions[t_hi] - positions[t_lo]) / (t_hi - t_lo)
            analytic_speed = sqrt2 * (0.5 - alpha)
            print 'speed = %g' % speed
            print 'analytic speed: %g' %  analytic_speed
            print 'speed error: %g' % (speed - analytic_speed)
            with open('data.txt', 'w') as f:
                f.write('%g' % (speed - analytic_speed))
            sys.exit()
        # wait for the forked process to complete
        os.wait()
        
        # get the value
        with open('data.txt') as f:
            speed_error = float(f.read())
        errors.append(speed_error)
        print 'back in parent: speed error: %g' % speed_error

    print alphas
    print errors
    all_errors.append(errors)


from matplotlib import pyplot
for errors, nseg in zip(all_errors, nsegs):
    pyplot.semilogy(alphas, numpy.abs(errors), label='dx=%g' % (float(length) / nseg))
pyplot.xlabel('alpha', fontsize=18)
pyplot.ylabel('error in wavespeed', fontsize=18)
pyplot.legend()
pyplot.setp(pyplot.gca().get_xticklabels(), fontsize=18)
pyplot.setp(pyplot.gca().get_yticklabels(), fontsize=18)
pyplot.show()

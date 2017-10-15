"""
A model of spatial working memory from Compte et al. (2000) implemented in Brian
"""

import brian_no_units
from brian import *
from numpy.fft import rfft,irfft
from numpy.random import *
import time
import numpy as np


# File names for monitors

spikes_eA_file = 'spikes_eA.dat'
spikes_iA_file = 'spikes_iA.dat'
spikes_eB_file = 'spikes_eB.dat'
spikes_iB_file = 'spikes_iB.dat'
lfpA_file= 'lfpA.dat'
lfpB_file= 'lfpB.dat'


execfile('parameters.py')
execfile('connectivity.py')
execfile('stimulus.py')

# Run main code

t0=time.clock()

# Define equations for neuron groups

eqs_eA = '''
dv/dt = (i_l + i_ext_eA +i_ee_BA +i_eeA + i_ieA + i_e) / Cm_e: volt

i_l = -gl_e*(v-El_e): amp
i_ext_eA = - g_ext_eA*s_ext*(v-E_ampa) : amp
i_ee_BA = - Gee_BA*s_sumB*(v-E_nmda)/(1+b*exp(-a*v)) : amp
i_eeA = - GeeA*s_tot*(v-E_nmda)/(1+b*exp(-a*v)) : amp
i_ieA = - GieA*s_gaba*(v-E_gaba) : amp

ds_ext/dt = -s_ext / t_ampa : 1
ds_rec_ampa/dt = -s_rec_ampa / t_ampa : 1
ds_cross/dt = -s_cross / t_ampa : 1
ds_gaba/dt = -s_gaba / t_gaba : 1
ds_nmda/dt = -s_nmda / t_nmda + alfa * x * (1 - s_nmda) : 1
dx/dt = -x / t_x : 1

i_lfp = np.abs(i_ext_eA) + np.abs(i_ee_BA) + np.abs(i_eeA) + np.abs(i_ieA) + np.abs(i_e) : amp

s_tot : 1
s_sumB: 1
i_e : amp
'''

eqs_iA = '''
dv/dt = (-gl_i*(v-El_i) - g_ext_iA*s_ext*(v-E_ampa) - Gei_BA*s_sumB*(v-E_nmda)/(1+b*exp(-a*v)) - GeiA*(v-E_nmda)*s_tot/(1+b*exp(-a*v)) - GiiA*s_gaba*(v-E_gaba) + i_e) / Cm_i: volt

ds_ext/dt = -s_ext / t_ampa : 1
ds_rec_ampa/dt = -s_rec_ampa / t_ampa : 1
ds_cross/dt = -s_cross / t_ampa : 1
ds_gaba/dt = -s_gaba / t_gaba : 1


s_tot : 1
s_sumB : 1
i_e : amp

'''

eqs_eB = '''
dv/dt = (i_l + i_ext_eB +i_ee_AB +i_eeB + i_ieB + i_e) / Cm_e: volt

i_l = -gl_e*(v-El_e) : amp
i_ext_eB = - g_ext_eB*s_ext*(v-E_ampa) : amp
i_ee_AB = - Gee_AB*s_sumA*(v-E_nmda)/(1+b*exp(-a*v)) : amp
i_eeB = - GeeB*s_tot*(v-E_nmda)/(1+b*exp(-a*v)) : amp
i_ieB = - GieB*s_gaba*(v-E_gaba) : amp

ds_ext/dt = -s_ext / t_ampa : 1
ds_rec_ampa/dt = -s_rec_ampa / t_ampa : 1
ds_cross/dt = -s_cross / t_ampa : 1
ds_gaba/dt = -s_gaba / t_gaba : 1
ds_nmda/dt = -s_nmda / t_nmda + alfa * x * (1 - s_nmda) : 1
dx/dt = -x / t_x : 1

i_lfp = np.abs(i_ext_eB) + np.abs(i_ee_AB) + np.abs(i_eeB) + np.abs(i_ieB) + np.abs(i_e) : amp

s_tot : 1
s_sumA : 1
i_e : amp
'''

eqs_iB = '''
dv/dt = (-gl_i*(v-El_i) - g_ext_iB*s_ext*(v-E_ampa) - Gei_AB*s_sumA*(v-E_nmda)/(1+b*exp(-a*v))- GeiB*(v-E_nmda)*s_tot/(1+b*exp(-a*v)) - GiiB*s_gaba*(v-E_gaba) + i_e) / Cm_i: volt

ds_ext/dt = -s_ext / t_ampa : 1
ds_rec_ampa/dt = -s_rec_ampa / t_ampa : 1
ds_gaba/dt = -s_gaba / t_gaba : 1

s_tot : 1
s_sumA :1
i_e : amp
'''

eqs_lfp = '''
lfp : amp
'''

### BEGIN Groups ###
print "Setting up the populations ..."

PeA = NeuronGroup(NE, eqs_eA, threshold=Vth_e, reset=Vr_e, refractory=tr_e, clock=simulation_clock, order=2, freeze=True)
PeA.v = El_e
PeA.s_ext = 0
PeA.s_gaba = 0
PeA.s_nmda = 0
PeA.s_rec_ampa = 0
PeA.x = 0
PeA.i_e = 0*nA

PiA = NeuronGroup(NI, eqs_iA, threshold=Vth_i, reset=Vr_i, refractory=tr_i, clock=simulation_clock, order=2, freeze=True)
PiA.v = El_i
PiA.s_ext = 0
PiA.s_rec_ampa = 0
PiA.s_gaba = 0
PiA.i_e = 0*nA

PeB = NeuronGroup(NE, eqs_eB, threshold=Vth_e, reset=Vr_e, refractory=tr_e, clock=simulation_clock, order=2, freeze=True)
PeB.v = El_e
PeB.s_ext = 0
PeB.s_gaba = 0
PeB.s_nmda = 0.5
PeB.s_rec_ampa = 0
PeB.x = 0
PeB.i_e = 0*nA

PiB = NeuronGroup(NI, eqs_iB, threshold=Vth_i, reset=Vr_i, refractory=tr_i, clock=simulation_clock, order=2, freeze=True)
PiB.v = El_i
PiB.s_ext = 0
PiB.s_rec_ampa = 0
PiB.s_gaba = 0
PiB.i_e = 0*nA

MA_lfp = NeuronGroup(1,eqs_lfp)
MA_lfp.lfp = 0

MB_lfp = NeuronGroup(1,eqs_lfp)
MB_lfp.lfp = 0


# external background Poisson input
PGeA  = PoissonGroup(NE, fext_eA, clock=simulation_clock)
PGiA  = PoissonGroup(NI, fext_iA, clock=simulation_clock)
PGeB  = PoissonGroup(NE, fext_eB, clock=simulation_clock)
PGiB  = PoissonGroup(NI, fext_iB, clock=simulation_clock)

### END Groups ###

### BEGIN Connections ###
print "Creating connections ..."

# NMDA update (presynaptic) via self-connection
selfnmdaA = IdentityConnection(PeA,PeA,'x',weight=1.0)
selfnmdaB = IdentityConnection(PeB,PeB,'x',weight=1.0)

# Ring A

# background poisson
CpeA = IdentityConnection(PGeA, PeA, 's_ext', weight=1.0)
CpiA = IdentityConnection(PGiA, PiA, 's_ext', weight=1.0)

# GABA
CieA = Connection(PiA, PeA, 's_gaba', weight=1.0)
CiiA = Connection(PiA, PiA, 's_gaba', weight=1.0)

# Ring B

# background poisson
CpeB = IdentityConnection(PGeB, PeB, 's_ext', weight=1.0)
CpiB = IdentityConnection(PGiB, PiB, 's_ext', weight=1.0)

# recurrent GABA
CieB = Connection(PiB, PeB, 's_gaba', weight=1.0)
CiiB = Connection(PiB, PiB, 's_gaba', weight=1.0)

### END Connections ###


# Update NMDA via FFT
@network_operation(simulation_clock, when='start')
def update_nmda(simulation_clock):
    # from A
    fsnmdaA = rfft(PeA.s_nmda)
    fstotA = fsnmdaA*fweightA
    # from A to A
    PeA.s_tot = irfft(fstotA, NE)
    PiA.s_tot = fsnmdaA[0]

    PeB.s_sumA = PeA.s_nmda.sum()
    PiB.s_sumA = PeA.s_nmda.sum()
    PeA.s_sumB = PeB.s_nmda.sum()
    PiA.s_sumB = PeB.s_nmda.sum()
    PeB.s_tot = PeB.s_nmda.sum()
    PiB.s_tot = PeB.s_nmda.sum()

# Update injected current stimulus input
@network_operation(current_clock, when='start')
def update_currents(current_clock):
    c_time = current_clock.t
    if (t_cue_start < c_time < t_cue_start + t_cue_dur):
        PeA.i_e = current_e
        PiA.i_e = 0*nA
        PeB.i_e = 0*nA
        PiB.i_e = 0*nA
    else:
        PeA.i_e = 0*nA
        PiA.i_e = 0*nA
        PeB.i_e = 0*nA
        PiB.i_e = 0*nA

dt_lfp = 1.*ms
lfp_clock = Clock(dt=dt_lfp)


@network_operation(lfp_clock,when='start')
def measure_lfp(lfp_clock):
    MA_lfp.lfp = mean(PeA.i_lfp)
    MB_lfp.lfp = mean(PeB.i_lfp)

# initiating monitors

# Record to file
traceA_lfp = StateMonitor(MA_lfp,'lfp',record=True)
traceB_lfp = StateMonitor(MB_lfp,'lfp',record=True)

spikes_eA = FileSpikeMonitor(PeA,spikes_eA_file,record=True)
spikes_iA = FileSpikeMonitor(PiA,spikes_iA_file,record=True)
spikes_eB = FileSpikeMonitor(PeB,spikes_eB_file,record=True)
spikes_iB = FileSpikeMonitor(PiB,spikes_iB_file,record=True)

M_eA=PopulationSpikeCounter(PeA)
M_iA=PopulationSpikeCounter(PiA)
M_eB=PopulationSpikeCounter(PeB)
M_iB=PopulationSpikeCounter(PiB)

print "Running ..."
run(simtime)
print "Done!"


np.savetxt(lfpA_file,np.vstack((traceA_lfp.times/second,traceA_lfp[0]/nA)).T,delimiter='\t')
np.savetxt(lfpB_file,np.vstack((traceB_lfp.times/second,traceB_lfp[0]/nA)).T,delimiter='\t')


rpop_eA = 1.*M_eA.nspikes/NE/simtime
rpop_iA = 1.*M_iA.nspikes/NI/simtime
rpop_eB = 1.*M_eB.nspikes/NE/simtime
rpop_iB = 1.*M_iB.nspikes/NI/simtime

print 'rate pop e A = ', rpop_eA
print 'rate pop e B = ', rpop_eB
print 'rate pop i A = ', rpop_iA
print 'rate pop i B = ', rpop_iB

figure()
subplot(211)
raster_plot(spikes_eA)
subplot(212)
raster_plot(spikes_eB)
show()

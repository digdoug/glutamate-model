
# size of the network
NE          = 2048               # number of E cells
NI          = 512                  # number of I cells

# simulation-related parameters
simtime     = 6000.*ms                # total simulation time [ms]
t_end = simtime
dt          = 0.02*ms         # simulation step length [ms]
simulation_clock=Clock(dt=dt)   # clock for the time steps

# pyramidal cells
Cm_e        = 0.5*nF  		# [nF] total capacitance
gl_e        = 25.0*nS 		# [ns] total leak conductance
El_e        = -70.0*mV		# [mV] leak reversal potential
Vth_e       = -50.0*mV		# [mV] threshold potential
Vr_e	    = -60.0*mV		# [mV] resting potential
tr_e	    = 2.0*ms		# [ms] refractory time

# interneuron cells
Cm_i        = 0.2*nF		# [nF] total capacitance
gl_i        = 20.0*nS		# [ns] total leak conductance
El_i        = -70.0*mV 		# [mV] leak reversal potential
Vth_i       = -50.0*mV 		# [mV] threshold potential
Vr_i	    = -60.0*mV 		# [mV] resting potential
tr_i	    = 1.0*ms  		# [ms] refractory time

# external background input
fext_eA      = 1800.0*Hz        # [Hz] external input frequency (poisson train)
fext_iA      = 1800.0*Hz
fext_eB      = 1800.0*Hz        # [Hz] external input frequency (poisson train)
fext_iB      = 1800.0*Hz

# AMPA receptor (APMAR)
E_ampa       = 0.0*mV          # [mV] synaptic reversial potential
t_ampa       = 2.0*ms          # [ms] exponential decay time constant
g_ext_eA     = 3.1*nS        # [nS] maximum conductance from external to pyramidal cells
g_ext_iA     = 2.38*nS         # [nS] maximum conductance from external to interneuron cells
g_ext_eB     = 3.1*nS           # [nS] maximum conductance from external to pyramidal cells
g_ext_iB     = 2.38*nS        # [nS] maximum conductance from external to interneuron cells

Gee_ampaA    = 0*nS#0.801*1024./NE*nS      # [nS] maximum conductance from pyramidal to pyramidal cells
Gei_ampaA    = 0*nS#0.684*1024./NE*nS      # [nS] maximum conductance from pyramidal to inhbitory cells
Gee_ampaB    = 0*nS#0.2486/NE*uS#0.391*nS      # [nS] maximum conductance from pyramidal to pyramidal cells
Gei_ampaB    = 0*nS#0.1958/NE*uS#0.293*nS      # [nS] maximum conductance from pyramidal to inhbitory cells

# GABA receptor (GABAR)
E_gaba          = -70.0*mV          # [mV] synaptic reversial potential
t_gaba          = 10.0*ms        # [ms] exponential decay time constant
GieA            = 1.4*1.336*2048/NE*nS         # [ns] synaptic conductance interneuron to pyramidal cells
GiiA            = 1.4*1.024*2048/NE*nS          # [ns] synaptic conductance interneuron to interneuron cells
GieB            = 1.4*1.336*2048/NE*nS          # [ns] synaptic conductance interneuron to pyramidal cells #0.99
GiiB            = 1.4*1.024*2048/NE*nS          # [ns] synaptic conductance interneuron to interneuron cells



# NMDA receptor (NMDAR)
E_nmda          = 0.0*mV         # [mV] synaptic reversial potential
t_nmda          = 100.0*ms        # [ms] decay time of NMDA currents
t_x             = 2.0*ms         # [ms] controls the rise time of NMDAR channels
alfa            = 0.5*kHz          # [kHz] controls the saturation properties of NMDAR channels
a               = 0.062/mV          # [1/mV] control the voltage dependance of NMDAR channel
b               = 1.0/3.57          # [1] control the voltage dependance of NMDAR channel ([Mg2+]=1mM )
GeeA            = 1.2*0.381*2048/NE*nS          # [ns] synaptic conductance from pyramidal to pyramidal cells
# this also got reduced by 1.25%
GeiA        = 0.8*0.292*2048/NE*nS          # [ns] synaptic conductance from pyramidal to interneuron cells
GeeB        = 1.02*1.2*0.381*2048/NE*nS         # [ns] synaptic conductance from pyramidal to pyramidal cells 1.8
# this is apparently getting reduced by 1.25%
GeiB        = 1.2*0.292*2048/NE*nS          # [ns] synaptic conductance from pyramidal to interneuron cells 1.8


Gei_AB = 200./NE*nS
Gei_BA = 60./NE*nS

Gee_AB = Gei_AB
Gee_BA = Gei_BA

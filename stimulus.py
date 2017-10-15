# cue parameters
 
t_cue_start      = 2000*ms          # start of the cue signal
t_cue_dur        = 1000*ms          # length period of the cue signal
i_cue_amp     = 200*pA       # cue current amplitude 180
i_cue_ang     = 180               # mean angle of the cue current
i_cue_width   = 14.4              # sigma/width of the current signal in degrees
dt_curr       = 50*ms             # current_clock step length (should be the highest common denominator between tc_start and tr_stop)
t_cue_stop       = t_cue_start + t_cue_dur
current_clock = Clock(dt=dt_curr)


# return the normalised "distance" of the neurons in a circle
def circ_distance(fi):
    if (fi > 0):
        return min(fi,360-fi)
    else:
        return max(fi,fi-360)
    

# calculating the cue currents
currents = lambda i,j: i_cue_amp*exp(-0.5*circ_distance((i-j)*360./NE)**2/i_cue_width**2)

current_e=zeros(NE)
j = i_cue_ang*NE/360.
for i in xrange(NE): 
    current_e[i]=currents(i,j)

    
   
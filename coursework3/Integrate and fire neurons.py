import numpy as np
import matplotlib.pyplot as plt
import random
import math

def intAndFire(leakPot, Vrest, Vthresh, Rm, Taum, Ie, duration, timestep):
    V0=Vrest
    data=[]
    x=[]
    for t in range(0,int(duration/timestep)):
        V=V0+((leakPot-V0+(Rm*Ie))/Taum)*timestep
        data.append(V/mV)
        x.append(t*timestep)
        if V>Vthresh:
            V0=Vrest
        else:
            V0=V

    plt.plot(x,data)
    plt.ylabel('V / mV')
    plt.xlabel('time / seconds')
    plt.show()
    return

def synapse2Neruons(Taum, E_L, Vrest, Vth, RmIe, RmGs, deltaS, Taus, Es, duration, timestep):
    V0n1 = random.randint(int(Vrest/mV), int(Vth/mV))*mV
    V0n2 = random.randint(int(Vrest/mV), int(Vth/mV))*mV
    S1to2_prev = 0
    S2to1_prev = 0
    RmIs1to2 = 0
    RmIs2to1 = 0
    data_N1 = []
    data_N1.append(V0n1)
    data_N2 = []
    data_N2.append(V0n2)
    x = []

    for t in range(1,1+int(duration/timestep)):
        #synapse n2 to n1
        if data_N2[t-1] == Vrest:
            S2to1 = 0.5 + S2to1_prev - (timestep*(S2to1_prev))/Taus
        else:
            S2to1 = S2to1_prev - (timestep*(S2to1_prev))/Taus

        #N1
        RmIs2to1 = RmGs*(Es - data_N2[t-1])*S2to1
        Vn1 = data_N1[t-1] + ((E_L- data_N1[t-1] + RmIe + RmIs2to1)/Taum)*timestep

        #synapse n1 to n2
        if data_N1[t-1] == Vrest:
            S1to2 = 0.5 + S1to2_prev - (timestep*(S1to2_prev))/Taus
        else:
            S1to2 = S1to2_prev - (timestep*(S1to2_prev))/Taus

        #N2
        RmIs1to2 = RmGs*(Es - data_N1[t-1])*S1to2
        Vn2 = data_N2[t-1] + ((E_L- data_N2[t-1] + RmIe + RmIs1to2)/Taum)*timestep

        if Vn1>Vth:
            Vn1=Vrest
        data_N1.append(Vn1)

        if Vn2>Vth:
            Vn2=Vrest
        data_N2.append(Vn2)

        S1to2_prev = S1to2
        S2to1_prev = S2to1
        x.append(t*timestep)

    data_N1[:] = [x / mV for x in data_N1]
    data_N2[:] = [x / mV for x in data_N2]
    plt.plot(x,data_N1[1:], label="Neuron 1")
    plt.plot(x,data_N2[1:], label="Neuron 2")
    plt.ylabel('V / mV')
    plt.xlabel('time / seconds')
    plt.legend(loc='upper right')
    plt.show()

    return

def STDP(E_L, Vrest, Vth, Rm, Taum, Ie, N_synapses, Taus, g_bar, Es, deltaS, timestep, firingRate, duration, STDP_flag):
    s_array = [0] * N_synapses
    g_array = [g_bar] * N_synapses
    data = []
    data.append(Vrest)   
    x = []

    t_pre = [0] * N_synapses
    t_post = -1000
    t_post_bool = False

    for t in range(1,1+int(duration/timestep)):

        sumG_bar_s = 0
        for i in range(0, N_synapses):
            sumG_bar_s += g_array[i]*s_array[i]

        RmIs = sumG_bar_s * Rm * (Es - data[t-1])
        v = data[t-1] + ((E_L- data[t-1] + (Rm*Ie) + RmIs)/Taum)*timestep

        for i in range (0, N_synapses):
            if random.random() < (firingRate * timestep):
                s_array[i] = 0.5 + s_array[i] - (timestep*s_array[i])/Taus
                t_pre[i] = t*timestep
            else:
                s_array[i] = s_array[i] - (timestep*s_array[i])/Taus

        if v>Vth:
            t_post = (t)*timestep
            t_post_bool = True
            v=Vrest
        data.append(v)
        x.append(t*timestep)

        if(STDP_flag == 'on'):
            for i in range (0, N_synapses):
                if t_post_bool == True :
                    g_array[i] = f_delta_t(g_array[i], (t_post - t_pre[i]), 0.2*nSeim, 0.25*nSeim, 20*ms, 20*ms)
                elif t_pre[i] == t*timestep :
                    g_array[i] = f_delta_t(g_array[i], (t_post - t_pre[i]), 0.2*nSeim, 0.25*nSeim, 20*ms, 20*ms)
        
        t_post_bool = False

    #data[:] = [x / mV for x in data]
    #plt.plot(x,data[1:])
    #plt.ylabel('V / mV')
    #plt.xlabel('time / seconds')
    plt.hist(g_array)
    plt.xlabel('Synaptic Strengths / nano-Seimens')
    plt.ylabel('Quantity')
    plt.show()
    
    return

def f_delta_t(g_bar,delta_t, A_plus, A_minus, Tau_plus, Tau_minus):
    new_g_bar = 0
    if delta_t > 0:
        new_g_bar = g_bar + A_plus*math.exp(-delta_t/Tau_plus)
    else:
        new_g_bar = g_bar -A_minus*math.exp(delta_t/Tau_minus)

    if(new_g_bar > 4*nSeim):
        new_g_bar = 4*nSeim
    if(new_g_bar < 0):
        new_g_bar = 0

    return new_g_bar

sec=1.0
ms=0.001
Volt=1.0
mV=0.001
Ohm=1.0
MOhm=1000000
Amp=1.0
nA=0.000000001
nSeim=0.000000001
hz=1.0

#PART A
#intAndFire(-70*mV, -70*mV, -40*mV, 10*MOhm, 10*ms, 3.1*nA, 1*sec, 0.25*ms)
#synapse2Neruons(20*ms, -70*mV, -80*mV, -54*mV, 18*mV, 0.15, 0.5, 10*ms, -80*mV, 1*sec, 0.25*ms)

#PART B
STDP(-65*mV, -65*mV, -50*mV, 100*MOhm, 10*ms, 0, 40, 2*ms, 4*nSeim, 0, 0.5, 0.25*ms, 15*hz, 300*sec, 'on')
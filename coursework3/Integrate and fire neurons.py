import numpy as np
import matplotlib.pyplot as plt
import random

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

sec=1.0
ms=0.001
Volt=1.0
mV=0.001
Ohm=1.0
MOhm=1000000
Amp=1.0
nA=0.000000001

#intAndFire(-70*mV, -70*mV, -40*mV, 10*MOhm, 10*ms, 3.1*nA, 1*sec, 0.25*ms)
synapse2Neruons(20*ms, -70*mV, -80*mV, -54*mV, 18*mV, 0.15, 0.5, 10*ms, -80*mV, 1*sec, 0.25*ms)

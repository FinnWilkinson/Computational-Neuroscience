import numpy as np
import matplotlib.pyplot as plt

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

sec=1.0
ms=0.001
Volt=1.0
mV=0.001
Ohm=1.0
MOhm=1000000
Amp=1.0
nA=0.000000001

intAndFire(-70*mV, -70*mV, -40*mV, 10*MOhm, 10*ms, 3.1*nA, 1*sec, 0.25*ms)

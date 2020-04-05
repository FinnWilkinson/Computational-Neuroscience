import random as rnd
import numpy as np

def get_spike_train(rate,big_t,tau_ref):

    if 1<=rate*tau_ref:
        print("firing rate not possible given refractory period f/p")
        return []


    exp_rate=rate/(1-tau_ref*rate)

    spike_train=[]

    t=rnd.expovariate(exp_rate)

    while t< big_t:
        spike_train.append(t)
        t+=tau_ref+rnd.expovariate(exp_rate)

    return spike_train

def get_fano_factor(spike_train, window_Width):
    current_count=0
    interval_counts=[]
    spike_train_index=0
    current_window = window_Width

    #get interval counts
    while spike_train_index < len(spike_train):
        if spike_train[spike_train_index] < (current_window):
            current_count+=1
            spike_train_index+=1
        else :
            current_window+=window_Width
            interval_counts.append(current_count)
            current_count = 0

    #get mean of interval counts
    mu=sum(interval_counts)/len(interval_counts)
    #get varience of interval counts
    sigma_squared=np.var(interval_counts)
    #calc fano_factor
    fano_factor = sigma_squared/mu
    return fano_factor

def get_coef_of_var(spike_train):
    intervals=[]
    for i in range(0, len(spike_train)-1):
        intervals.append(spike_train[i+1]-spike_train[i])

    #get mean
    mu=sum(intervals)/len(intervals)
    #get standard deviation
    sigma=np.std(intervals)
    #calculate Coefficient of variation
    coef_of_var = sigma/mu

    return coef_of_var


Hz=1.0
sec=1.0
ms=0.001

rate=35.0 *Hz #firing rate
tau_ref=5*ms #refactory period
window_Width=100*ms

big_t=1000*sec #seconds of spike train

spike_train=get_spike_train(rate,big_t,tau_ref)
fano_factor = get_fano_factor(spike_train, window_Width)
coef_of_var = get_coef_of_var(spike_train)


print(len(spike_train)/big_t)
print(spike_train)

print("Spike train length: {}"  .format(len(spike_train)))
print("Firing rate: {}" .format(rate))
print("Refractory Period: {}" .format(tau_ref))
print("Window Width: {}" .format(window_Width))
print("Fano Factor: {}"  .format(fano_factor))
print("coefficient of varience: {}" .format(coef_of_var))

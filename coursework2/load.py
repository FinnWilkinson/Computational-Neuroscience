import numpy as np
import matplotlib.pyplot as plt

def load_data(filename,T):

    data_array = [T(line.strip()) for line in open(filename, 'r')]

    return data_array

def get_fano_factor(spike_train, window_Width, time_between_samples):
    interval_counts=[]
    current_count=0

    for i in range(0, len(spike_train)):
        if spike_train[i] == 1:
            current_count+=1
        else:
            pass

        if (i%(window_Width/time_between_samples) == 4):
            interval_counts.append(current_count)
            current_count=0


    #get mean of interval counts
    mu=sum(interval_counts)/len(interval_counts)
    #get varience of interval counts
    sigma_squared=np.var(interval_counts)
    #calc fano_factor
    fano_factor = sigma_squared/mu
    return fano_factor

def get_coef_of_var(spike_train, time_between_samples):
    intervals=[]
    last_spike_pos=0
    for i in range(0, len(spike_train)):
        if(spike_train[i]==1):
            intervals.append((i-last_spike_pos)*time_between_samples)
            last_spike_pos=i

    #get mean
    mu=sum(intervals)/len(intervals)
    #get standard deviation
    sigma=np.std(intervals)
    #calculate Coefficient of variation
    coef_of_var = sigma/mu

    return coef_of_var


#spikes=[int(x) for x in load_data("rho.dat")]
spikes=load_data("rho.dat",int)

print(len(spikes))
print(spikes[0:5])

#stimulus=[float(x) for x in load_data("stim.dat")]
stimulus=load_data("stim.dat",float)

print(len(stimulus))
print(stimulus[0:5])

ms=0.001
window_Width=100*ms
time_between_samples = 2*ms

fano_factor=get_fano_factor(spikes, window_Width, time_between_samples)
coef_of_var=get_coef_of_var(spikes, time_between_samples)
print("Window Width: {}" .format(window_Width))
print("Fano Factor: {}" .format(fano_factor))
print("Coefficient of Varience: {}" .format(coef_of_var))

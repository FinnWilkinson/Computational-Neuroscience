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

def get_autocorrelogram(spike_train, look_range):
    positiveSpikes=[0]*int(look_range/2)
    negativeSpikes=[0]*int(look_range/2)
    spikeCount=0
    for i in range(0, len(spike_train)):
        if(spike_train[i]==1):
            spikeCount+=1
            for j in range(1,int(look_range/2)+1):
                if(i-j>=0):
                    if(spike_train[i-j]==1):
                        negativeSpikes[j-1]+=1
                if(i+j<len(spike_train)):
                    if(spike_train[i+j]==1):
                        positiveSpikes[j-1]+=1

    spikeCount+=1
    negativeSpikes.reverse()
    negativeSpikes.append(spikeCount)
    for x in positiveSpikes:
        negativeSpikes.append(x)
    finalData=[x/spikeCount for x in negativeSpikes]
    fig, ax = plt.subplots()
    ax.bar([i*2 for i in range(-int(look_range/2), int(look_range/2)+1)], finalData )
    ax.set_ylabel('Normalised Number of Spikes at Each Interval')
    ax.set_xlabel('Interval (ms)')
    ax.set_title('Autocorrelogram of Fly H1 Neuron Responding to Approximate White-Noise Visual Motion Stimulus - Over Range -100ms to +100ms')
    plt.show()
    return

def get_spike_trigger_average(spike_train, stimulus, look_range):
    spikeCount = 0
    data=[0]*(int(look_range/2)+1)
    for i in range(0, len(spike_train)):
        if spike_train[i] == 1:
            spikeCount+=1
            for j in range(0, len(data)):
                if i-j >=0 :
                    data[j]+=stimulus[i-j]
    print(data)
    data.reverse()
    finalData=[x/spikeCount for x in data]
    print(finalData)
    fig, ax = plt.subplots()
    ax.bar([i*2 for i in range(-len(data)+1, 1)], finalData)
    ax.set_ylabel('Normalised Spike Triggered Average')
    ax.set_xlabel('Interval (ms)')
    ax.set_title('Spike Triggered Average Over a 100ms Window for Fly H1 Neuron Responding to Approximate White-Noise Visual Motion Stimulus')
    plt.show()
    return


#spikes=[int(x) for x in load_data("rho.dat")]
spikes=load_data("rho.dat",int)

print(len(spikes))
print(spikes[10:20])

#stimulus=[float(x) for x in load_data("stim.dat")]
stimulus=load_data("stim.dat",float)

print(len(stimulus))
print(stimulus[10:20])

ms=0.001
window_Width=100*ms
time_between_samples = 2*ms

fano_factor=get_fano_factor(spikes, window_Width, time_between_samples)
coef_of_var=get_coef_of_var(spikes, time_between_samples)
print("Window Width: {}" .format(window_Width))
print("Fano Factor: {}" .format(fano_factor))
print("Coefficient of Varience: {}" .format(coef_of_var))
get_autocorrelogram(spikes, 100)
get_spike_trigger_average(spikes, stimulus, 100)

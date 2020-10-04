import numpy as np
import pyaudio
import os
import wave, struct
import matplotlib.pyplot as plt
from scipy.io import wavfile

def binary(A, T,time,signal_rate, folder):
    
    signal = []
    amount = int(signal_rate/T)
    times =int(time*T)

    for elem in range(0, times):
        shot = np.random.rand()
        if shot>0.5:
            for elem in range(0, amount):
                signal.append(1)
        else:
            for elem in range(0, amount):
                signal.append(0)
    times = np.linspace(0, time, len(signal))
    signal = A*np.array(signal)
    plt.figure(figsize=(30, 4))
    plt.plot(times, signal) 
    plt.xlim(times[0], times[-1])
    plt.xlabel('time (s)')
    plt.ylabel('amplitude')
    # You can set the format by changing the extension
    # like .pdf, .svg, .eps
    pat = os.getcwd() + "\\Sound" + "\\" + folder + "\\Recording.png"
    plt.savefig(pat, dpi=100)
    plt.show()
    signal = np.array([times,signal])
    return signal
def record(time,rate, pather):
    set_chunk()
    # the file name output you want to record into
    filename = pather
    # set the chunk size of 1024 samples
    # sample format
    FORMAT = pyaudio.paInt16
    # mono, change to 2 if you want stereo
    channels = 1
    # 44100 samples per second
    sample_rate = rate
    record_seconds = time
    # initialize PyAudio object
    p = pyaudio.PyAudio()
    # open stream object as input & output
    stream = p.open(format=FORMAT,
                    channels=channels,
                    rate=sample_rate,
                    input=True,
                    output=True,
                    frames_per_buffer=chunk)
    frames = []
    print("Recording...")
    for i in range(int(sample_rate / chunk * record_seconds)):
        data = stream.read(chunk)
        # if you want to hear your voice while recording
        stream.write(data)
        frames.append(data)
    print("Finished recording.")
    # stop and close stream
    stream.stop_stream()
    stream.close()
    # terminate pyaudio object
    p.terminate()
    # save audio file
    # open the file in 'write bytes' mode
    wf = wave.open(filename, "wb")
    # set the channels
    wf.setnchannels(channels)
    # set the sample format
    wf.setsampwidth(p.get_sample_size(FORMAT))
    # set the sample rate
    wf.setframerate(sample_rate)
    # write the frames as bytes
    wf.writeframes(b"".join(frames))
    # close the file
    wf.close()
def playback(path):
    print("MOODY BLUES :REWIND")
    CHUNK = chunk
    wf = wave.open(path, 'rb')

    # instantiate PyAudio (1)
    p = pyaudio.PyAudio()

    # open stream (2)
    stream = p.open(format=p.get_format_from_width(wf.getsampwidth()),
                    channels=wf.getnchannels(),
                    rate=wf.getframerate(),
                    output=True)

    # read data
    data = wf.readframes(CHUNK)

    # play stream (3)
    while len(data) > 0:
        stream.write(data)
        data = wf.readframes(CHUNK)

    # stop stream (4)
    stream.stop_stream()
    stream.close()

    # close PyAudio (5)
    p.terminate()
def waveform(pather, folder):
    samplerate, data = wavfile.read(pather)
    times = np.arange(len(data))/float(samplerate)

    # Make the plot
    # You can tweak the figsize (width, height) in inches
    plt.figure(figsize=(30, 4))
    plt.fill_between(times, data) 
    plt.xlim(times[0], times[-1])
    plt.xlabel('time (s)')
    plt.ylabel('amplitude')
    # You can set the format by changing the extension
    # like .pdf, .svg, .eps
    pat = os.getcwd() + "\\Sound" + "\\" + folder + "\\recording.png"
    plt.savefig(pat, dpi=100)
    plt.show()
    sound = np.array([times,data])
    return sound
def play(sound, foldername, name):
    chunk = 1000
    sampleRate = 50000.0 # hertz
    duration = 2.0 # seconds

    obj = wave.open(os.getcwd() +"\\Sound\\" + foldername+ '\\' + name +  '.wav','w')
    obj.setnchannels(1) # mono
    obj.setsampwidth(2)
    obj.setframerate(sampleRate)
    for i in range(len(sound)):
        value = int(sound[i])
        data = struct.pack('<h', value)
        obj.writeframesraw( data )
    obj.close()
    playback(os.getcwd() +"\\Sound\\" + foldername+ '\\' + name +  '.wav')



def set_chunk():
    global chunk
    
    chunk = 1000


def Setter(folder, mode, **kwargs): #Provide Folder Name(inside sound Directory) #Provide a name for the wav file #Mode 1- binary 2 - record
    set_chunk()
    name = "Recording"
    for key in kwargs:
        if (key=='A'): #Amplitude
            A = float(kwargs[key])
        if (key=='T'): #Period
            T = float(kwargs[key])
        if (key=='rate'):
            rate = int(kwargs[key])
        if (key=='time'):
            time = int(kwargs[key])
     

    path = os.getcwd() + "\\Sound\\" + str(folder)
    try:
        os.mkdir(path)
    except:
        pass
    pather = path + "\\" + name
    if mode == 1:
        signal = binary(A,T,time,rate, folder)
    else:
        record(T,rate, pather)   
        playback(pather)
        signal = waveform(pather,folder)
    return signal


#Setter('Tester', 2, A=5, T=0.1)

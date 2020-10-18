import wave, struct, math, random
import os
import numpy as np
import pyaudio
import os
import wave, struct
import matplotlib.pyplot as plt
from scipy.io import wavfile

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



chunk = 1000
sampleRate = 44100.0 # hertz
duration = 2.0 # seconds
#frequency = 440.0 # hertz
obj = wave.open('sound.wav','w')
obj.setnchannels(1) # mono
obj.setsampwidth(2)
obj.setframerate(sampleRate)
for i in range(99999):
   value = random.randint(-32767, 32767)
   data = struct.pack('<h', value)
   obj.writeframesraw( data )
obj.close()
playback(os.getcwd() + "\sound.wav")

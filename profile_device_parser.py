import pickle
import numpy as np
import matplotlib.pyplot as plt
import os
from math import log10

def AV_fft(fs_p, signal_p, is_subtract_mean_p=False):
    """
    The function returns only the positive-frequency term.
    See https://docs.scipy.org/doc/scipy/reference/generated/scipy.fftpack.fft.html

    This function fixs the bug in AV_fft(). That is when len(signal_p) is even, in AV_fft() we forget to include the result at the nyquist frequency fs_p/2.
    But the result is correct when len(signal_p) is odd.
    """
    try:
        print("AV_fft: Wrong signal_p dimension %d" % (int(np.shape(signal_p)[1])))
    except IndexError:    

        if is_subtract_mean_p:
            signal_p = copy.deepcopy(signal_p - np.mean(signal_p))

        N_signal_l = len(signal_p)
        fft_l = np.fft.fft(signal_p)
        file_l = np.fft.fftfreq(N_signal_l, d = (1/fs_p))

        if ((N_signal_l % 2) == 0):
            freq_positive_l = np.abs(file_l[0:( (N_signal_l//2) + 1 )])
            fft_l           = fft_l[0:( (N_signal_l//2) + 1)]
        else:
            freq_positive_l = file_l[0:(((N_signal_l - 1)//2) + 1)]
            fft_l           = fft_l[0:(((N_signal_l - 1)//2) + 1)]


        fft_amp_l = 2.0*np.abs(fft_l)/N_signal_l
        fft_amp_l[0] /= 2.0 # We divide by 2.0 because there is no negative frequency of 0 Hz
        if ((N_signal_l % 2) == 0):
            fft_amp_l[-1] /= 2.0 # We divide by 2.0 because there is no negative frequency of the Nyquist frequency.

        fft_pw_l = fft_amp_l**2
        fft_gain_l = np.abs(fft_l)

        # fft_amp_l has the same unit as signal_p. [fft_amp_l] = [signal_p]
        # [fft_pw_l] = [fft_amp_l]*[fft_amp_l], i.e. if the unit of the signal is microvolt, the unit of the power (fft_pw_l) is microvolt^2
        # fft_gain_l is calculated based on the amplitude. Therefore, to convert to dB, we do 20*np.log10(fft_gain_l) because dB = 10*np.log10(power_output/power_input)
        # -3dB refers to the a frequency with which the power of output is reduced by half or the amplitude is reduced from 1 to sqrt(2).
        return freq_positive_l, fft_l, fft_amp_l, fft_pw_l, fft_gain_l

def areaUnderCurve(x,y) :
  areaUnderCurve = sum(y)
  return areaUnderCurve

def getAmplitudeAllFreq(directory,parsing_range) :
  freq_quantity = len(os.listdir(directory))
  rFreq_l = list(range(1,freq_quantity+1))
  file_l = []
  rAmp_l = [-1]*freq_quantity

  print(directory)
  
  for filename in [filename for filename in os.listdir(directory) if filename != 'gain' ]:
    file_l.append(filename)

  #for loop in each freq
  for filename in file_l : 
    data = loadData(directory,filename)
    profile = EEGProfile(data=data)
    this_freq = profile.func_gen_file_l-1

    [x_pw,y_pw] = profile.getPeakPowerGraph(scope=21)

    area = areaUnderCurve(x_pw,y_pw)

    currentAmp = area**0.5
    rAmp_l[this_freq] = currentAmp

  return rFreq_l[0:parsing_range],rAmp_l[0:parsing_range]

def loadData(directory,freq_inspect) :
  fileName = '%s/%s'%(directory,freq_inspect)
  f = open(fileName,'rb')
  data = pickle.load(f)
  f.close()
  return data

def openGainFile(path) :
  data = {}
  try :
    f = open(path,'rb')
    data = pickle.load(f)
    # print("open",data)
    f.close()
  except :
    pass
  return data

def saveGainFile(path,gain) :
  if not os.path.exists(os.path.dirname(path)):
    try:
        os.makedirs(os.path.dirname(path))
    except OSError as exc: # Guard against race condition
        if exc.errno != errno.EEXIST:
            raise
  f = open(path,'wb')
  pickle.dump(gain,f)
  f.close()
  print('save %s complete'%(path))


def getDevicePath(root_path,deviceName,deviceCH) :
  devices_path = getAllDevicePath(root_path)
  for device_path in devices_path :
    ext = device_path.strip('.').strip('./').split('/')
    device_name = ext[0]
    device_ch = ext[-1] if len(ext) > 1 else 'AF7'
    if (deviceName == device_name and deviceCH == device_ch) :
      return device_path

def getAllDevicePath(this_directory) :
  file_l = os.listdir(this_directory)
  device_l = [device for device in file_l if os.path.isdir('%s/%s'%(this_directory,device)) and device[0] != '_' ]
  if (len(device_l) != 0) :
    for device in device_l :
      next_directory = '%s/%s'%(this_directory,device)
      yield from getAllDevicePath(next_directory)
  else :
    yield this_directory

def ampToVoltage(data_l,gain_l,parsing_range) :
  data_inverse = []
  freq_l = []
  for i in range(parsing_range) :
    gainFactor = 10**(gain_l[i]/20)
    current_data_inverse = data_l[i]/gainFactor
    data_inverse.append(current_data_inverse)
    freq_l.append(i+1)
  return data_inverse,freq_l

def voltageToAmp(data_l,gain_l,parsing_range) :
  data_new = []
  freq_l = []
  for i in range(parsing_range) :
    gainFactor = 10**(gain_l[i]/20)
    current_data_new = data_l[i]*gainFactor
    data_new.append(current_data_new)
    freq_l.append(i+1)
  return data_new,freq_l


def parseGainToList(gain_l) :
  return [gain[1] for gain in gain_l]

class EEGProfile : 
  def __init__(
    self,
    data,
  ) :
    self.data = data
    self.Fs_l = data[0]
    self.samples_duration_l = data[1]
    self.sorted_samples_l = data[2]
    self.func_gen_amp_l = data[3]
    self.func_gen_file_l = data[4]
    self.func_gen_offs_l = data[5]
    self.func_gen_duty_l = data[6]
    self.func_gen_phase_l = data[7]
    [freq_positive_l, fft_l, fft_amp_l, fft_pw_l, fft_gain_l] = AV_fft(fs_p=self.Fs_l,signal_p=self.sorted_samples_l)
    self.freq_positive_l = freq_positive_l
    self.fft_l = fft_l
    self.fft_amp_l = fft_amp_l
    self.fft_pw_l = fft_pw_l
    self.fft_gain_l = fft_gain_l

  def __str__(self) :
    return self.data
  
  def toDiagram(self) :
        fig, axs = plt.subplots(4)
        fig.suptitle('%s - rawValue ,fft_amp_l , fft_pw_l - freq %s Hz'%("mindwave",self.func_gen_file_l))
        axs[0].plot(self.sorted_samples_l)
        axs[1].set_title("amp")
        axs[1].plot(self.freq_positive_l,self.fft_amp_l,'go-')
        axs[2].set_title("gain")
        axs[2].plot(self.freq_positive_l,self.fft_gain_l,'go-')
        axs[3].set_title("pw")
        axs[3].plot(self.freq_positive_l,self.fft_pw_l,'go-')
        plt.show()

  def getFFTIndexByFreq(self) :
    return self.func_gen_file_l*60 

  def getFFTGain(self) :
    return self.fft_gain_l

  def getPeakPositiveFreq(self,scope=7) :
    index = self.getFFTIndexByFreq()
    delta = scope//2
    return self.freq_positive_l[index-delta:index+delta+1] 

  def getPeakPower(self,scope=7) :
    index = self.getFFTIndexByFreq()
    delta = scope//2
    return self.fft_pw_l[index-delta:index+delta+1] 
  
  def getPeakAmp(self,scope=7) :
    index = self.getFFTIndexByFreq()
    delta = scope//2
    return self.fft_amp_l[index-delta:index+delta+1] 
  
  def getPeakPowerGraph(self,scope=7) :
    x = self.getPeakPositiveFreq(scope)
    y_pw = self.getPeakPower(scope)
    return x,y_pw

  def getPeakAmplitudeGraph(self,scope=7) :
    x = self.getPeakPositiveFreq(scope)
    y_amp = self.getPeakAmp(scope)
    return x,y_amp


INPUT_VOLTAGE = 1 #mV
GAIN_ROOT_PATH = './_gains'
SEP = ', '
DEVICE_ROOT_PATH = '.'

gain_d = {}
devices_support = os.listdir(GAIN_ROOT_PATH)

for device in devices_support :
  device_path = '%s/%s'%(GAIN_ROOT_PATH,device)
  tmp = openGainFile(device_path)
  gain_d[device] = tmp


#print device info
print()
print("--------------------------------------")
for device in devices_support :
  print('%s ( %s ) support in freq range 1 - %d'%(
    device,
    ', '.join(gain_d[device]["ch"].keys()),
    len(gain_d[device]["ch"]["AF7"])
    ))
print("--------------------------------------")
print()


root_file_l = os.listdir(DEVICE_ROOT_PATH)
device_l = [device for device in devices_support if os.path.isdir('./%s'%device) and device[0].isupper()]


#-------------------------- select master device section --------------------------
print("pls select master device data")
device_default = device_l
device_master = input('%s (%s) : '%(SEP.join(device_l),device_default[0])).strip().capitalize() or device_default[0]
device_master_channel = [ *gain_d[device_master]["ch"] ]

if ( len(device_master_channel)==1 ) :
  print("use device channel %s "%(device_master_channel[0]))
  device_master_channel = device_master_channel[0]
else :
  print("pls select master device channel")
  device_master_channel = input('%s (%s) : '%(SEP.join(device_master_channel),device_master_channel[0])).strip() or device_master_channel[0]

#-------------------------- select master device section --------------------------




#-------------------------- select slave device section  --------------------------
print()
print("which device format you want to parse to ?")
device_default = device_l
device_parse_to = input('%s (%s) : '%(SEP.join(device_l),device_default[0])).strip().capitalize() or device_default[0]
device_parse_to_channel = [ *gain_d[device_parse_to]["ch"] ]

if ( len(device_parse_to_channel)==1 ) :
  print("use device channel %s "%(device_parse_to_channel[0]))
  device_parse_to_channel = device_parse_to_channel[0]
else :
  print("pls select parse to device channel")
  device_parse_to_channel = input('%s (%s) : '%(SEP.join(device_parse_to_channel),device_parse_to_channel[0])).strip() or device_parse_to_channel[0]

#-------------------------- select slave device section  --------------------------


parsing_range = min(
  len( gain_d[device_master]["ch"][device_master_channel] ), 
  len( gain_d[device_parse_to]["ch"][device_parse_to_channel] )
)


#show summary before simulate
print()
print("--------------------------------------")
print('(original) %s ( %s ) support in freq range 1 - %d'%(device_master,device_master_channel,len(gain_d[device_master]["ch"][device_master_channel])))
print('(parse to) %s ( %s ) support in freq range 1 - %d'%(device_parse_to,device_parse_to_channel,len(gain_d[device_parse_to]["ch"][device_parse_to_channel])))
print('simulate only 1 - %d Hz'%( parsing_range ))
print("--------------------------------------")
print()
input("comfirm (ENTER)")


#get device path
master_device_path = getDevicePath(DEVICE_ROOT_PATH, device_master, device_master_channel)
parse_to_device_path = getDevicePath(DEVICE_ROOT_PATH, device_parse_to, device_parse_to_channel)

#get device gain
master_gain = parseGainToList( gain_d[device_master]["ch"][device_master_channel] )
parse_to_gain = parseGainToList( gain_d[device_parse_to]["ch"][device_parse_to_channel] )

#get amplitude from raw data
master_freq_l,master_amp_l = getAmplitudeAllFreq(master_device_path,parsing_range)
parse_to_freq_l,parse_to_amp_l = getAmplitudeAllFreq(parse_to_device_path,parsing_range)

#parsing
master_voltage,master_voltage_freq = ampToVoltage(master_amp_l,master_gain,parsing_range)
parsed_freq,parsed_amp = voltageToAmp(master_voltage,parse_to_gain,parsing_range)



#test convert simulate back to master
test_inverse,test_inverse_freq = ampToVoltage(parsed_freq,parse_to_gain,parsing_range)
test_output,test_output_freq = voltageToAmp(test_inverse,master_gain,parsing_range)



#plot
fig, axs = plt.subplots(3,2)
fig.suptitle('parse amp from %s (%s) to %s (%s)'%(
  device_master, device_master_channel,
  device_parse_to, device_parse_to_channel
))
axs[0,0].set_title('amp master (%s)'%(device_master))
axs[0,0].plot(master_freq_l,master_amp_l,'go-')
axs[1,0].set_title("voltage master")
axs[1,0].plot(master_voltage_freq,master_voltage,'go-')
axs[2,0].set_title('result -simulate %s from %s'%(device_parse_to,device_master))
axs[2,0].plot(parsed_amp,parsed_freq,'go-')
axs[2,1].set_title('simulate ideal')
axs[2,1].plot(parse_to_freq_l,parse_to_amp_l,'go-')
axs[0,1].set_title('test -simulate %s from simulated %s'%(device_master,device_parse_to))
axs[0,1].plot(test_output_freq,test_output,'go-')
plt.show()




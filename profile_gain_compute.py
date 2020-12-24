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

def merge_two_dicts(x, y):
  z = x.copy()   # start with x's keys and values
  z.update(y)    # modifies z with y's keys and values & returns None
  return z

def loadData(directory,freq_inspect) :
  fileName = '%s/%s'%(directory,freq_inspect)
  f = open(fileName,'rb')
  data = pickle.load(f)
  f.close()
  return data

def openGainFile(root_path,device_name) :
  data = {
    "ch" : {

    }
  }
  path = '%s/%s'%(root_path,device_name)
  try :
    f = open(path,'rb')
    data = pickle.load(f)
    f.close()
  except :
    pass

  return data

def saveGainFile(device_name,device_ch,new_gain_l,device_profile) :
  global gain_d
  global GAIN_ROOT_PATH
  gain_d = openGainFile(GAIN_ROOT_PATH,device_name)
  gain_d["ch"][device_ch] = [(i+1,new_gain_l[i]) for i in range(len(new_gain_l))]
  gain_data = merge_two_dicts(gain_d,device_profile)

  save_path = '%s/%s'%(GAIN_ROOT_PATH,device_name)
  if not os.path.exists(os.path.dirname(save_path)):
    try:
        os.makedirs(os.path.dirname(save_path))
    except OSError as exc: # Guard against race condition
        if exc.errno != errno.EEXIST:
            raise
  f = open(save_path,'wb')
  pickle.dump(gain_data,f)
  f.close()
  print('save %s complete'%(save_path))

def getGainAllFreq(directory) :
  freq_quantity = len(os.listdir(directory))
  ext = directory.strip('.').strip('./').split('/')
  device_name = ext[0]
  device_ch = ext[-1] if len(ext) > 1 else 'AF7'
  device_profile = {}
  freq_l = list(range(1,freq_quantity+1))
  file_l = []
  gain_l = []
  # freq_inspect = input("freq inspect (all) : ").strip()
  
  for filename in [filename for filename in os.listdir(directory) if filename != 'gain' ]:
    file_l.append(filename)
  gain_l = [-1]*freq_quantity

  for filename in file_l : 
    data = loadData(directory,filename)
    # print(data)
    profile = EEGProfile(data=data)
    device_profile["params"] = {
      "amplitude" : {
        "value" : profile.func_gen_amp_l,
        "unit" : "V"
      },
      "offset" : {
        "value" :  profile.func_gen_offs_l,
        "unit" : "V"
      },
      "duty" : {
        "value" : profile.func_gen_duty_l,
        "unit" : None
      },
      "phase" : {
        "value" : profile.func_gen_phase_l,
        "unit" : "degree"
      }, 
      "desc" : "Gain of amplitude "
    }
    this_freq = profile.func_gen_freq_l-1

    [x_pw,y_pw] = profile.getPeakPowerGraph(scope=21)
    area = areaUnderCurve(x_pw,y_pw)

    currentGain = 20*log10(area**0.5/INPUT_VOLTAGE)
    gain_l[this_freq] = currentGain
  
  return freq_l,gain_l,device_name,device_ch,device_profile


def getDevicePath(DEVICE_ROOT_PATH,device_name) :
  filePathTmp = '%s/%s'%(DEVICE_ROOT_PATH,device_name)

  while True :
    this_directory = os.listdir(filePathTmp)
    if len(this_directory) == 1 :
      this_file = this_directory[0]
    elif not os.path.isdir(filePathTmp+'/%s'%this_directory[0]) :
      this_file = this_directory[0] 
    else :
      print("pls select file in",filePathTmp)
      this_file = input( '%s (%s) : '%(
        SEP.join(this_directory),
        this_directory[0]) 
      ) or this_directory[0]

    next_directory = filePathTmp + '/%s'%this_file
    isDir = os.path.isdir(next_directory)
    if not isDir:
      return filePathTmp
    else :
      filePathTmp = next_directory

def getAllDevicePath(this_directory) :
  file_l = os.listdir(this_directory)
  device_l = [device for device in file_l if os.path.isdir('%s/%s'%(this_directory,device)) and device[0] != '_' ]
  if (len(device_l) != 0) :
    for device in device_l :
      next_directory = '%s/%s'%(this_directory,device)
      yield from getAllDevicePath(next_directory)
  else :
    yield this_directory

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
    self.func_gen_freq_l = data[4]
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
        fig.suptitle('%s - rawValue ,fft_amp_l , fft_pw_l - freq %s Hz'%("mindwave",self.func_gen_freq_l))
        axs[0].plot(self.sorted_samples_l)
        axs[1].set_title("amp")
        axs[1].plot(self.freq_positive_l,self.fft_amp_l,'go-')
        axs[2].set_title("gain")
        axs[2].plot(self.freq_positive_l,self.fft_gain_l,'go-')
        axs[3].set_title("pw")
        axs[3].plot(self.freq_positive_l,self.fft_pw_l,'go-')
        plt.show()

  def getFFTIndexByFreq(self) :
    return self.func_gen_freq_l*60 

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
SEP = ', '
GAIN_ROOT_PATH = './_gains'
DEVICE_ROOT_PATH = '.'

all_devices_path = getAllDevicePath(DEVICE_ROOT_PATH)
root_file_l = os.listdir(DEVICE_ROOT_PATH)
device_l = [device for device in root_file_l if os.path.isdir('%s/%s'%(DEVICE_ROOT_PATH,device)) and device[0].isupper()]


print("pls select device data to compute gain")
device_default = device_l
device_name = input('%s (%s) : '%(SEP.join(device_l),"all")).strip().capitalize() or device_default
if (device_name == 'All') :
  device_name = device_default


if (isinstance(device_name,str)) :
  directory = getDevicePath(DEVICE_ROOT_PATH,device_name)

  freq_l, new_gain_l, device_name, device_ch,device_profile = getGainAllFreq(directory)
  print(device_profile)

  saveGainFile(device_name,device_ch,new_gain_l,device_profile)

  plt.plot(freq_l,new_gain_l)
  plt.ylabel('gain')
  plt.xlabel('freq')
  plt.show()
else :
  for device_path in all_devices_path :
    print('data found %s'%(device_path))
    freq_l, new_gain_l, device_name, device_ch,device_profile = getGainAllFreq(device_path)

    saveGainFile(device_name,device_ch,new_gain_l,device_profile)
    print()
  
  input("enter to close")
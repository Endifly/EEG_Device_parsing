import pickle
import numpy as np
import matplotlib.pyplot as plt
import os

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

def loadData(device_path,file_name) :
  fileName = '%s/%s'%(device_path,file_name)
  f = open(fileName,'rb')
  data = pickle.load(f)
  f.close()
  return data

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
    fft_gain = self.getFFTGain()
    y_pw = self.getPeakPower(scope)
    return x,y_pw,fft_gain

SEP = ', '
DEVICE_ROOT_PATH = '.'

root_file_l = os.listdir(DEVICE_ROOT_PATH)
device_l = [device for device in root_file_l if os.path.isdir('%s/%s'%(DEVICE_ROOT_PATH, device)) and device[0].isupper()]


print("pls select device data")
device_default = device_l[0]
deviceName = input('%s (%s) : '%(
  SEP.join(device_l),
  device_default
  )).strip().capitalize() or device_default
device_path = getDevicePath(DEVICE_ROOT_PATH,deviceName)


freq_quantity = len(os.listdir(device_path))
freq_l = list(range(1,freq_quantity+1))
file_l = []
freq_inspect = input("freq inspect (all) : ").strip()

if (freq_inspect == "" or freq_inspect == 'all' ) :
  for filename in [filename for filename in os.listdir(device_path) if filename != 'gain' ]:
    file_l.append(filename)
  plt_pw_l = [-1]*freq_quantity
  plt_amp_l = [-1]*freq_quantity
  plt_gain_l = [-1]*freq_quantity

  for filename in file_l : 
    data = loadData(device_path,filename)
    profile = EEGProfile(data=data)
    this_freq = profile.func_gen_freq_l-1

    [x,y_pw,gain] = profile.getPeakPowerGraph(scope=21)
    area = areaUnderCurve(x,y_pw)

    plt_pw_l[this_freq]=area
    plt_amp_l[this_freq]=area**0.5
    plt_gain_l[this_freq]=gain

  print("amp_l",plt_amp_l)
  fig, axs = plt.subplots(2)
  fig.suptitle('%s'%(device_path.strip('.').strip('/').replace('/',' - ')))
  axs[0].set_title("fft_amp by freq")
  axs[0].plot(freq_l,plt_amp_l,'go-')
  axs[1].set_title("fft_pw by freq")
  axs[1].plot(freq_l,plt_pw_l,'go-')
  plt.show()
  
else :
  filename_l = os.listdir(device_path)
  for filename in filename_l :
    data = loadData(device_path,filename)
    profile = EEGProfile(data=data)
    if (profile.func_gen_freq_l == int(freq_inspect)) :
      profile.toDiagram()
      break
from t14re import T14re, Config, ConvertDataFromBin2Mat
import datetime


t14re = T14re()
# T14RE 2D
cfgname = 'T14RE_2D_100fps.cfg'
# T14RE 3D
cfgname = 'T14RE_3D_100fps.cfg'

t14re.load_and_send_config(cfgname)


CaptureDeviceIndex = 0 # Disable any cameras except for a radar. Check Windows Device Manager.
t14re.connect_and_setup_radar(CaptureDeviceIndex)


t14re.capture(0)

t14re.disconnect()



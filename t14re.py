import serial
from serial.tools import list_ports
import cv2
import time
import numpy as np
import datetime
import queue
import scipy.io
import threading
import os
import glob
import sys

class T14re():
    def __init__(self):
        self.config = 0
        self.q = queue.Queue(-1)
        self.q_name = queue.Queue(-1)
        self.isStopped = False
        self.isConnected = False
        self.cap = 0
        self.DefaultFolderName = 'data_t14re'
        self.framecnt = 0
        self.getfilecnt = 0
        self.savefilecnt = 0
        self.MaxFiles = 0
    
    def get_com_port(self):
        ports = list_ports.comports()
        if not ports:
            print('Error: No COM ports. Check connection.',file=sys.stderr)
        else:
            for p in ports:
                if 'USB' in p.description and p.location[-1:] == '2': # for MMIC (MI_02)
                    print('Selected port: ' + p.device)
                    return p.device
        return False
        
    def load_and_send_config(self,cfgname):
        self.config = Config()
        COM_port_for_MMIC = self.get_com_port()
        self.config.send_config(cfgname,COM_port_for_MMIC)
    
    def connect_and_setup_radar(self,CaptureDeviceIndex):
        # get VideoCapture object
        # Window: settings for using Microsoft Media Foundation
        width = self.config.profileCfg['numAdcSamples'] * self.config.NRx * 2
        height = self.config.NumOfChirpsPerChirpset * self.config.frameCfg['ChirpsetsPerFrame'] * self.config.FramesPerData
        print('width: ' + str(width))
        print('height: ' + str(height))
        
        camera_num = self.count_camera_num()
        if camera_num > 1:
            print('\nCAUTION: %d cameras are connected.' % camera_num)
            s = 0
            for _ in range(32):
                s = input('\nPlease disable any cameras except for a radar from the Device Manager and set CaptureDeviceIndex=0.\n'
                'Or, have you set an appropriate CaptureDeviceIndex? [y/n]: ')
                if s=='y':
                    break
                elif s=='n':
                    return False
                else:
                    print('\nInvalid input!')

        self.cap = cv2.VideoCapture(CaptureDeviceIndex, cv2.CAP_MSMF)
        TF1 = self.cap.set(cv2.CAP_PROP_CONVERT_RGB, 0)
        TF2 = self.cap.set(cv2.CAP_PROP_FORMAT, -1)
        TF3 = self.cap.set(cv2.CAP_PROP_FRAME_WIDTH, width)
        TF4 = self.cap.set(cv2.CAP_PROP_FRAME_HEIGHT, height)
        if TF1 & TF2 & TF3 & TF4 & self.cap.isOpened():
            self.isConnected = True
            print('\nConnection has succeeded!')
        else:
            self.isConnected = False
            print('\nError: Connection failed. Please unplug cables and reconnect.',file=sys.stderr)
        return True
        
    def get_and_put_data(self):
        while not self.isStopped:
            self.getfilecnt += 1
            ret, data = self.cap.read()
            self.q.put(data)
            filename = datetime.datetime.now().strftime('%Y-%m-%d-%H-%M-%S-%f_')
            filename = filename + str(int(self.config.Fs)) + 'fps_' + str(self.getfilecnt) + '.bin'
            self.q_name.put(filename)
            print('\rfilecnt:',self.getfilecnt, end='')
        print('\nFinish')
    
    def save_data(self):
        '''
        Save data with binary.
        '''
        foldername = datetime.datetime.now().strftime('%Y-%m-%d-%H-%M-%S-%f_t14re')
        os.makedirs(self.DefaultFolderName + '/' + foldername)
        print('foldername:',foldername)
        
        while not self.isStopped:
            if self.q.qsize() > 1:
                for _ in range(self.q.qsize()):
                    filename = self.q_name.get(block=False)
                    f = open(self.DefaultFolderName + '/' + foldername + '/' + filename,mode='ab')
                    frame = self.q.get(block=False)
                    f.write(frame.tobytes())
                    f.close()
                    self.savefilecnt += 1
                    if self.getfilecnt >= self.MaxFiles and self.MaxFiles > 0:
                        self.isStopped = True
                        break
                        
        while self.savefilecnt < self.MaxFiles and self.MaxFiles > 0:
            filename = self.q_name.get(block=False)
            f = open(self.DefaultFolderName + '/' + foldername + '/' + filename,mode='ab')
            frame = self.q.get(block=False)
            f.write(frame.tobytes())
            f.close()
            self.savefilecnt += 1
        
    def capture(self,MaxFiles=0,startDatetime=0,endDatetime=0):
        '''
        Start capturing
        if MaxFiles >= 1, radar captures MaxFiles files.
        if MaxFiles <  1, PC waits for q key to stop (be careful with the data size).
        '''
        self.isStopped = False
        self.getfilecnt = 0
        self.savefilecnt = 0
        while not self.q.empty():
            self.q.get(block=False)
        while not self.q_name.empty():
            self.q_name.get(block=False)
        
        if startDatetime==0:
            startDatetime = datetime.datetime.now()
        
        if self.waitUntilGivenDatetime(startDatetime):
            now = datetime.datetime.now()
            d = now.strftime('%Y/%m/%d %H:%M:%S')
            print('Start datetime:',d)
            
            if MaxFiles<1:
                self.MaxFiles = -1
                if endDatetime==0:
                    print('Capture until press q')
                    thread1 = threading.Thread(target=self.get_and_put_data)
                    thread2 = threading.Thread(target=self.save_data)
                    thread3 = threading.Thread(target=self.PressQuitKeyEvent)
                    thread1.start()
                    thread2.start()
                    thread3.start()
                    thread1.join()
                    thread2.join()
                    thread3.join()
                else:
                    print('Capture until the specified time')
                    thread1 = threading.Thread(target=self.get_and_put_data)
                    thread2 = threading.Thread(target=self.save_data)
                    thread3 = threading.Thread(target=self.stopAtSpecifiedTime,args=(endDatetime,))
                    thread1.start()
                    thread2.start()
                    thread3.start()
                    thread1.join()
                    thread2.join()
                    thread3.join()
            else:
                print('Capture for decided number of files')
                self.MaxFiles = MaxFiles
                thread1 = threading.Thread(target=self.get_and_put_data)
                thread2 = threading.Thread(target=self.save_data)
                thread1.start()
                thread2.start()
                thread1.join()
                thread2.join()

    def disconnect(self):
        '''
        Disconnect radar.
        '''
        if self.isConnected:
            self.cap.release()
            self.isConnected = False
            print('Disconnect.')
        else:
            print('Error: Already disconnected.', file=sys.stderr)
    
    def count_camera_num(self):
        cameraNum = 0
        for cameraId in range(5):
            cap = cv2.VideoCapture(cameraId,cv2.CAP_MSMF)
            if cap.isOpened():
                cameraNum += 1
            cap.release()
        return cameraNum
    
    def PressQuitKeyEvent(self):
        s = 0
        while s!='q':
            s = input('\nPress q to stop capturing: ')
        self.isStopped = True
        
    def waitUntilGivenDatetime(self,Datetime):
        while True:
            diff_start = datetime.datetime.now() - Datetime
            if diff_start.days < 0:
                time.sleep(1)
            else:
                break
        return True
    
    def stopAtSpecifiedTime(self,endDatetime):
        if self.waitUntilGivenDatetime(endDatetime):
            now = datetime.datetime.now()
            d = now.strftime('%Y/%m/%d %H:%M:%S')
            print('\nStop datetime:',d)
            self.isStopped = True
        return True
            
class Config():
    def __init__(self):
        self.profileCfg = {}
        self.frameCfg = {}
        self.FramesPerData = 0
        self.NumOfChirpsPerChirpset = 0
        self.NRx = 0
        self.NTx = 0
        self.Fs = 0
        self.TxOrder = []
    
    def load_config(self,cfgname):
        # get contents of config file
        file = open(cfgname)
        commands = file.readlines()
        file.close()
        
        for i in range(len(commands)):
            tmp_cfg = commands[i].replace('\r','')
            tmp_cfg = commands[i].replace('\n','')
            tmp_cfg = tmp_cfg.split()
            if tmp_cfg[0]=='profileCfg':
                self.profileCfg['profileId'] = int(tmp_cfg[1])
                self.profileCfg['startFreq'] = float(tmp_cfg[2]) # GHz
                self.profileCfg['idleTime'] = float(tmp_cfg[3]) # us
                self.profileCfg['adcStartTime'] = float(tmp_cfg[4]) # us
                self.profileCfg['rampEndTime'] = float(tmp_cfg[5]) # us
                self.profileCfg['txOutPower'] = float(tmp_cfg[6])
                self.profileCfg['txPhaseShifter'] = float(tmp_cfg[7])
                self.profileCfg['freqSlopeConst'] = float(tmp_cfg[8]) # MHz/us
                self.profileCfg['txStartTime'] = float(tmp_cfg[9]) # us
                self.profileCfg['numAdcSamples'] = int(tmp_cfg[10]) 
                self.profileCfg['digOutSampleRate'] = int(tmp_cfg[11]) # kilo samples per sec

                hpfCornerFreq1_idx = int(tmp_cfg[12])
                if hpfCornerFreq1_idx==0:
                    self.profileCfg['hpfCornerFreq1'] = 175 # kHz
                elif hpfCornerFreq1_idx==1:
                    self.profileCfg['hpfCornerFreq1'] = 235
                elif hpfCornerFreq1_idx==2:
                    self.profileCfg['hpfCornerFreq1'] = 350
                elif hpfCornerFreq1_idx==3:
                    self.profileCfg['hpfCornerFreq1'] = 700

                hpfCornerFreq2_idx = int(tmp_cfg[13])
                if hpfCornerFreq2_idx==0:
                    self.profileCfg['hpfCornerFreq2'] = 350 # kHz
                elif hpfCornerFreq2_idx==1:
                    self.profileCfg['hpfCornerFreq2'] = 700
                elif hpfCornerFreq2_idx==2:
                    self.profileCfg['hpfCornerFreq2'] = 1400
                elif hpfCornerFreq2_idx==3:
                    self.profileCfg['hpfCornerFreq2'] = 2800
                
                self.profileCfg['rxGain'] = float(tmp_cfg[14])
                
            elif tmp_cfg[0]=='frameCfg':
                self.frameCfg['ChirpStartIndex'] = int(tmp_cfg[1])
                self.frameCfg['ChirpEndIndex'] = int(tmp_cfg[2])
                self.frameCfg['ChirpsetsPerFrame'] = int(tmp_cfg[3]) # same as 'number of loops' in TI document
                self.frameCfg['NumberOfFrames'] = int(tmp_cfg[4]) # 0 means infinite
                self.frameCfg['Periodicity'] = float(tmp_cfg[5]) # ms
                self.Fs = 1/(self.frameCfg['Periodicity']*1e-3) # framerate 
                self.frameCfg['TriggerSelect'] = int(tmp_cfg[6]) # 1: software 2: hardware
                self.frameCfg['FrameTriggerDelay'] = float(tmp_cfg[7]) # ms
            
            elif tmp_cfg[0]=='channelCfg':
                self.NRx = self.count_ones_by_bin(int(tmp_cfg[1]))
                self.NTx = self.count_ones_by_bin(int(tmp_cfg[2]))
                
            elif tmp_cfg[0]=='chirpCfg':
                self.NumOfChirpsPerChirpset += 1
                self.TxOrder.append(tmp_cfg[-1])
                
            elif 'framesPerData' in tmp_cfg[0]:
                self.FramesPerData = int(tmp_cfg[-1])
        return commands
    
    def send_config(self,cfgname,COM_port):
        commands = self.load_config(cfgname)
        
        baudrate = 115200
        serialCommand = serial.Serial(COM_port,baudrate)
        serialCommand.reset_output_buffer()
        serialCommand.reset_input_buffer()

        # send command
        for command in commands:
            if self.SendCommand(serialCommand, command) == False:
                print('Error: Failed to send command.', file=sys.stderr)
                serialCommand.close()
                return False

        serialCommand.close()
        print("Configuration data was sent successfully")
        return True
    
    def SendCommand(self,serialCommand,command):
        # send command replacing linebreak character appropriately
        command = command.replace('\r', '')
        command = command.replace('\n', '')
        command += '\n'
        serialCommand.write(command.encode())

        # wait for response
        ret = serialCommand.read_until(b'\nmmwDemo:/>')

        # confirm response
        # blank, comment, sensorReset do not output Done
        if (
            (command.strip() != '') and
            (command.strip()[0] != '%') and
            (command != 'sensorReset\n')
        ):
            rets = ret.decode().split('\n')
            if (rets[len(rets) - 2] != 'Done'):
                return False
        return True
    
    def count_ones_by_bin(self,num):
        bin_num = bin(num)[2:]
        count = 0
        for i in bin_num:
            count += int(i)
        return count
        
def ConvertDataFromBin2Mat(cfgname,path=0,save_dir='data_t14re'):
    config = Config()
    config.load_config(cfgname)

    targ_dirs = os.listdir(path=save_dir)
    if path==0:
        targ_dir = save_dir + "/" + targ_dirs[-1]
    else:
        targ_dir = path

    files = glob.glob(targ_dir + "*/*.bin")
    print("Target directory:",targ_dir)
    num_file = len(files)
    for i in range(num_file):
        print('\rfilenum: '+str(i+1)+'/'+str(num_file), end='')
        fname = files[i]
        fr = open(fname, 'rb')
        data = np.frombuffer(fr.read(),np.uint8)
        
        # data2 = data.reshape(FramesPerData,ChirpsetsPerFrame,NTx*NRx,SamplesPerChirp,2,2)
        data2 = data.reshape(config.FramesPerData,config.frameCfg['ChirpsetsPerFrame'],config.NTx*config.NRx,config.profileCfg['numAdcSamples'],2,2)
        
        I_data = (data2[:, :, :, :, 0, 0] + data2[:, :, :, :, 0, 1] * 256).astype(np.int16)
        Q_data = (data2[:, :, :, :, 1, 0] + data2[:, :, :, :, 1, 1] * 256).astype(np.int16)
        
        I_data[I_data >= 0x8000] -= 0x10000
        Q_data[Q_data >= 0x8000] -= 0x10000
        
        data_reshape = I_data + 1j * Q_data
        
        # save mat file
        scipy.io.savemat(fname.replace('.bin','') + ".mat", {'signals': data_reshape})
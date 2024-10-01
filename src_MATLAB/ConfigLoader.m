classdef ConfigLoader
    properties
        profileCfg = struct();
        frameCfg = struct();
        FramesPerData = 0;
        NumOfChirpsPerChirpset = 0;
        NRx = 0;
        NTx = 0;
        Fs = 0;
        TxOrder = [];
    end

    methods
        function self = load_config(self, cfgname)
            % get contents of config file
            fileID = fopen(cfgname, 'r');
            commands = textscan(fileID, '%s', 'Delimiter', '\n');
            commands = commands{1};
            fclose(fileID);

            for i = 1:length(commands)
                tmp_cfg = strsplit(commands{i});

                if strcmp(tmp_cfg{1},'profileCfg')
                    self.profileCfg.profileId = str2double(tmp_cfg{2});
                    self.profileCfg.startFreq = str2double(tmp_cfg{3}); % GHz
                    self.profileCfg.idleTime = str2double(tmp_cfg{4}); % us
                    self.profileCfg.adcStartTime = str2double(tmp_cfg{5}); % us
                    self.profileCfg.rampEndTime = str2double(tmp_cfg{6}); % us
                    self.profileCfg.txOutPower = str2double(tmp_cfg{7});
                    self.profileCfg.txPhaseShifter = str2double(tmp_cfg{8});
                    self.profileCfg.freqSlopeConst = str2double(tmp_cfg{9}); % MHz/us
                    self.profileCfg.txStartTime = str2double(tmp_cfg{10}); % us
                    self.profileCfg.numAdcSamples = str2double(tmp_cfg{11});
                    self.profileCfg.digOutSampleRate = str2double(tmp_cfg{12}); % kilo samples per sec

                    hpfCornerFreq1_idx = str2double(tmp_cfg{13});
                    hpfCornerFreq1 = [175, 235, 350, 700];
                    self.profileCfg.hpfCornerFreq1 = hpfCornerFreq1(hpfCornerFreq1_idx + 1);

                    hpfCornerFreq2_idx = str2double(tmp_cfg{14});
                    hpfCornerFreq2 = [350, 700, 1400, 2800];
                    self.profileCfg.hpfCornerFreq2 = hpfCornerFreq2(hpfCornerFreq2_idx + 1);

                    self.profileCfg.rxGain = str2double(tmp_cfg{15});

                elseif strcmp(tmp_cfg{1},'frameCfg')
                    self.frameCfg.ChirpStartIndex = str2double(tmp_cfg{2});
                    self.frameCfg.ChirpEndIndex = str2double(tmp_cfg{3});
                    self.frameCfg.ChirpsetsPerFrame = str2double(tmp_cfg{4}); % same as 'number of loops' in TI document
                    self.frameCfg.NumberOfFrames = str2double(tmp_cfg{5}); % 0 means infinite
                    self.frameCfg.Periodicity = str2double(tmp_cfg{6}); % ms
                    self.Fs = 1/(self.frameCfg.Periodicity*1e-3); % framerate
                    self.frameCfg.TriggerSelect = str2double(tmp_cfg{7}); % 1: software 2: hardware
                    self.frameCfg.FrameTriggerDelay = str2double(tmp_cfg{8}); % ms

                elseif strcmp(tmp_cfg{1},'channelCfg')
                    self.NRx = sum(dec2bin(str2double(tmp_cfg{2})) == '1');
                    self.NTx = sum(dec2bin(str2double(tmp_cfg{3})) == '1');

                elseif strcmp(tmp_cfg{1},'chirpCfg')
                    self.NumOfChirpsPerChirpset = self.NumOfChirpsPerChirpset + 1;
                    self.TxOrder = [self.TxOrder; str2double(tmp_cfg{end})];

                elseif contains(tmp_cfg{1}, 'framesPerData')
                    self.FramesPerData = str2double(tmp_cfg{end});
                end
            end
        end
    end
end

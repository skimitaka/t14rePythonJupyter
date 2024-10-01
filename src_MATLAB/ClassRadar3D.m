classdef ClassRadar3D
    %%% Needfiles
    % - mystft.m ()
    properties
        X
        Y
        Theta_deg
        Phi_deg
        noise_level = 30
        c = 3*1e8
        % f = 79*1e9
        f = (80.764+77.410)/2*1e9 % 厳密な中心周波数
        lmd % wave length
        d0 % element spacing
        % B = 3.6*1e9
        B = 3.354*1e9  % 厳密な帯域幅 (using high sampling config)
        dr % range resolution
        NTx
        NRx
        MS = 15
        N = 256
        dR = 0.0447
        range0
        fps = 100
        dt;
        ecg_freqmin = 0.8;
        ecg_freqmax = 12.5;

        syncTime
    end
    methods
        %% コンストラクタ
        function obj = ClassRadar3D(obj)
            obj.lmd = obj.c/obj.f;
            obj.d0 = obj.lmd/2;
            obj.dr = obj.c/(2*obj.B);

            obj.dt = 1/obj.fps;
        end

        %%
        function [obj, radar_data] = loadRadarBinData(obj, folderName, options)
            arguments (Input)
                obj
                folderName
                options.StartTime (1,1) double {mustBeNonnegative} = 0
                options.IntervalTime
                options.RangeCut
                options.NumFFTAngle
                options.NumFFTElevation
                options.FFTEnable=true
                options.configFile = "T14RE_3D_100fps.cfg"
            end

            % load config
            configFile = options.configFile;
            loader = ConfigLoader();
            loader = loader.load_config(configFile);
            matrix_size = [loader.FramesPerData, loader.frameCfg.ChirpsetsPerFrame, loader.NTx, loader.NRx, loader.profileCfg.numAdcSamples, 2, 2];

            radar_sample_rate = loader.Fs;
            FramesPerData = loader.FramesPerData;     % data sets per radar data
            ChirpsetsPerFrame = loader.frameCfg.ChirpsetsPerFrame;
            numAdcSamples = loader.profileCfg.numAdcSamples; % Range
            obj.NRx = loader.NRx; % Angle
            obj.NTx = loader.NTx; % Elevation
                              % the number of receive antennas
            if isfield(options, "NumFFTAngle")
                N_afft = options.NumFFTAngle;                      % angle dFt length
            else
                N_afft = obj.NRx*8;
            end
            if isfield(options, "NumFFTElevation")
                N_efft = options.NumFFTElevation;                      % angle dFt length
            else
                N_efft = obj.NTx*8;
            end
            Number_of_range = loader.profileCfg.numAdcSamples;
            N_rfft = 256; % range

            matFiles = fullfile(folderName, '*.bin');
            fileList = dir(matFiles);
            fileNum = size(fileList,1);

            if isfield(options, "RangeCut")
                delete_range = options.RangeCut;
            else
                delete_range = Number_of_range;
            end

            dr_pad = obj.dr * numAdcSamples / N_rfft;
            % range = (0:1:N_rfft-1)*dr_pad; %
            range = (0:delete_range-1)*dr_pad; %
            obj.range0 =  range;
            theta = asin(obj.lmd*(-N_afft/2:N_afft/2-1)./(N_afft*obj.d0)); %
            phi = asin(obj.lmd*(-N_efft/2:N_efft/2-1)./(N_efft*obj.d0)); %
            obj.Theta_deg = rad2deg( theta );
            obj.Phi_deg = rad2deg( phi );
            % [Range,Theta] = meshgrid(range,theta);
            % obj.X = Range.*sin(Theta);
            % obj.Y = Range.*cos(Theta);

            startFileId = floor( options.StartTime * (obj.fps / FramesPerData) + 1 );
            if isfield(options, "IntervalTime")
                loadFileNum = ceil( options.IntervalTime * (obj.fps / FramesPerData) );
                loadFileNum = min(loadFileNum, fileNum-startFileId+1);
            else
                loadFileNum = fileNum - startFileId + 1;
            end

            if options.FFTEnable
                radar_data = zeros(loadFileNum*FramesPerData, N_efft, N_afft, delete_range); % 配列の割り当て
            else
                radar_data = zeros(loadFileNum*FramesPerData, obj.NTx, obj.NRx, delete_range); % 配列の割り当て
            end

            strlen = 0;
            fprintf('Pleese wait ......... (startFileId=%d)\n', startFileId);
            fi=1;
            for FileID=startFileId:startFileId+loadFileNum-1
                % ログの表示
                Tmp = {'Current trial: %3d/%d\n', fi, loadFileNum};
                Tmp{1} = [ repmat(sprintf('\b'),[1 strlen]),  Tmp{1} ];
                Txt = sprintf(Tmp{1:3});
                fprintf(Txt);
                strlen = length(Txt) - strlen;

                filepath = fullfile(folderName, fileList(FileID).name);
                signals= obj.radardataconvert(char(filepath), matrix_size);

                signals = squeeze(mean(signals,2)); % coherent integration

                if options.FFTEnable
                    sig_rfft = fft(signals,N_rfft,4);
                    sig_rfft_afft = fftshift(fft(sig_rfft, N_afft,3), 3);
                    sig_rfft_afft_efft = fftshift(fft(sig_rfft_afft, N_efft,2), 2);
                else
                    sig_rfft_afft_efft = signals;
                end

                radar_data((fi-1)*FramesPerData+1:fi*FramesPerData, :, :,:) = sig_rfft_afft_efft(:, :, :,1:delete_range);

                fi = fi+1;
            end
            fprintf('................ done\n');
        end
        %%
        function [obj, radar_data] = loadRadarData(obj, folderName, options)
            arguments (Input)
                obj
                folderName
                options.StartTime (1,1) double {mustBeNonnegative} = 0
                options.IntervalTime
            end
            matFiles = fullfile(folderName, '*.mat');
            fileInfo = dir(matFiles);
            fileNum = size(fileInfo,1);

            N_rfft = 256;
            % N_afft = 92;
            N_afft = 48;
            % N_afft = 12;
            FS = 16;
            % delete_range = 72; % 3m
            delete_range = 256; % 3m

            I = 1;
            filepath = fullfile(folderName, fileInfo(I).name);
            load(filepath);
            [FramesPerData,ChirpsetsPerFrame,N,numAdcSamples] = size(signals);
            dr_pad = obj.dr * numAdcSamples / N_rfft;
            % range = (0:1:N_rfft-1)*dr_pad; %
            range = (0:1:delete_range-1)*dr_pad; %
            theta = asin(lmd*(-N_afft/2:N_afft/2-1)./(N_afft*d0)); %
            % theta = asin(linspace(1,-1,N_afft));
            [Range,Theta] = meshgrid(range,theta);
            obj.X = Range.*sin(Theta);
            obj.Y = Range.*cos(Theta);

            startFileId = floor( options.StartTime * (obj.fps / FramesPerData) + 1 );
            if isfield(options, "IntervalTime")
                loadFileNum = ceil( options.IntervalTime * (obj.fps / FramesPerData) );
                loadFileNum = min(loadFileNum, fileNum);
            else
                loadFileNum = fileNum;
            end

            radar_data = zeros(loadFileNum*FramesPerData, N_afft, delete_range); % 配列の割り当て
            strlen = 0;
            fprintf('Pleese wait .........\n');
            fi=1;
            for FileID=startFileId:startFileId+loadFileNum-1
                % ログの表示
                Tmp = {'Current trial: %3d/%d\n', fi, loadFileNum};
                Tmp{1} = [ repmat(sprintf('\b'),[1 strlen]),  Tmp{1} ];
                Txt = sprintf(Tmp{1:3});
                fprintf(Txt);
                strlen = length(Txt) - strlen;

                filepath = fullfile(folderName, fileInfo(FileID).name);
                load(filepath);
                signals = squeeze(mean(signals,2)); % coherent integration
                sig_rfft = fft(signals,N_rfft,3);
                sig_rfft_afft = fftshift(fft(sig_rfft,N_afft,2), 2);
                radar_data((fi-1)*FramesPerData+1:fi*FramesPerData, :, :) = sig_rfft_afft(:, :, 1:delete_range);

                fi = fi+1;
            end
            fprintf('................ done\n');
        end
        %%
        function showRadarImage(obj, Z, options)
            arguments (Input)
                obj
                Z % radar_data power (angle_num, range_num)
                options.Figure_id=1
                options.NormalizedMode=true
                options.Title
                options.Subplot = false
                options.DinamicRange = 30
            end

            if options.Subplot
                figure(gcf);
                nexttile;
            else
                figure(options.Figure_id);set(gcf,'color','white', 'visible', 'on');
                clf(gcf);
            end

            if options.NormalizedMode
                Z=10*log10(Z./max(Z,[],"all"));
            else
                Z=10*log10(Z);
            end

            % clf(gcf);
            surf(obj.X, obj.Y, Z, 'EdgeColor','none');
            view([0 0 1])
            % surfObj = pcolor(obj.X, obj.Y, Z);
            % surfObj.EdgeColor='none';
            daspect([1 1 1]);
            xlim([min(min(obj.X)),max(max(obj.X))]);
            ylim([0, obj.range0(end)]);
            xlabel('X (m)');
            ylabel('Y (m)');
            cb = colorbar;
            cb.Label.String = 'Intensity (dB)';

            if options.NormalizedMode
                clim([-options.DinamicRange, 0]);
            else
                clim([-options.DinamicRange, 0] + max(max(Z)));
            end


            colormap("jet");

            if isfield(options, "Title")
                title(options.Title);
            end
        end

        function [Z, v, idx] = showRadarIndex(obj, Z, options)
            arguments (Input)
                obj
                Z
                options.Figure_id=1
                options.Title
                options.Subplot=false
                options.DinamicRange = 30
            end

            if options.Subplot
                figure(gcf);
                nexttile;
            else
                figure(options.Figure_id);set(gcf,'color','white', 'visible', 'on');
                clf(gcf);
            end

            Z=10*log10(Z./max(Z,[],"all"));
            surf(Z.','EdgeColor','none');
            xlabel('axis index');
            ylabel('range index');
            cb = colorbar;
            clim([-options.DinamicRange, 0]);
            pbaspect([1 1 1]);
            view([0 0 1])
            colormap("jet");
            if isfield(options, "Title")
                title(options.Title);
            end

            % 電力最大点
            [v,idx] = obj.max_2d(Z);
            txt = sprintf('max (%d, %d)', idx(1), idx(2));
            text(1, 1.01, txt, 'Units','normalized')
        end
        %%
        function showMeanProcess(obj, radar_data, options)
            arguments (Input)
                obj
                radar_data
                options.Figure_id=1
            end
            figure_i=options.Figure_id;
            figure(figure_i);clf(gcf);set(gcf,'color','white', 'visible', 'on');
            mean_data = mean(radar_data, 1);
            radar_data_deaverage = radar_data - mean_data;

            %　時間平均削除後
            Z = squeeze(mean(abs(radar_data_deaverage),1));
            title_txt = "時間平均削除後";
            obj.showRadarImage(Z, Figure_id=figure_i, Title=title_txt, Subplot=true);
            obj.showRadarIndex(Z, Figure_id=figure_i, Title=title_txt, Subplot=true);

            %　時間平均削除前
            Z = squeeze(mean(abs(radar_data),1));
            title_txt = "時間平均削除前";
            obj.showRadarImage(Z, Figure_id=figure_i, Title=title_txt, Subplot=true);
            obj.showRadarIndex(Z, Figure_id=figure_i, Title=title_txt, Subplot=true);

            % %　時間平均画像
            % Z = 20*log10(squeeze(mean(abs(radar_data),1)));
            % title_txt = "時間平均";
            % obj.showRadarImage(Z, Figure_id=figure_i, Title=title_txt, Subplot=true);
            % obj.showRadarIndex(Z, Figure_id=figure_i, Title=title_txt, Subplot=true);
        end
        %%
        function showMovie(obj, radar_data, drawing_speed, options)
            arguments (Input)
                obj
                radar_data
                drawing_speed
                options.Figure_id=1
                options.Title
                options.StartTime=0
                options.IntervalTime=size(radar_data,1)
                options.SaveMovieFn

            end
            start_time_index = ceil(options.StartTime * obj.fps) + 1;
            end_time_index = start_time_index + ceil(options.IntervalTime * obj.fps) - 1;
            save_movie_flg = isfield(options, "SaveMovieFn");
            fprintf('(Start) %.2fsec - %.2fsec(End)\n', start_time_index/obj.fps, end_time_index/obj.fps);
            if save_movie_flg % ビデオの記録
                fprintf('SaveMovieMode On\n');
                v = VideoWriter(options.SaveMovieFn,'MPEG-4');
                cleanupObj = onCleanup(@() obj.cleanupFunction(v));
                open(v);
            end

            figObj = figure(options.Figure_id);clf(gcf);set(gcf,'Visible','on');
            axObj = gca;
            ti = 1;
            Z = 20*log10(abs(squeeze(radar_data(ti,:,:))));
            surfObj = surf(obj.X, obj.Y, Z,'EdgeColor','none');
            surfObj.ZDataSource = 'Z';
            xlabel('x (m)');
            ylabel('y (m)');
            cb = colorbar;
            % colormap parula;
            colormap jet;
            cb.Label.String = 'Intensity (dB)';
            xlim([min(min(obj.X)),max(max(obj.X))]);
            ylim([0, obj.range0(end)]);
            clim([-50,0]+max(max(Z)));
            daspect([1,1,1]);
            view(2);
            for ti=start_time_index:drawing_speed:end_time_index
                Z = 20*log10(abs(squeeze(radar_data(ti,:,:))));
                title(axObj, sprintf('%.1f (s) [ID %d]', ti/obj.fps, ti));
                refreshdata(surfObj, 'caller')
                drawnow

                if save_movie_flg % ビデオの記録
                    frame = getframe(gcf);
                    writeVideo(v,frame);
                end

            end

            if save_movie_flg % ビデオの記録
                close(v);
            end
        end
        %%
        function showRespECG(obj, id_x, id_y, iq_wave, options)
            arguments (Input)
                obj
                id_x
                id_y
                iq_wave
                options.Subplot
                options.Figure_id=1
                options.Title
            end
            if isfield(options, "Subplot")
                subplots = options.Subplot;
                figure(options.Figure_id);set(gcf,'color','white', 'visible', 'on');
                subplot(subplots(1), subplots(2), subplots(3));
            else
                figure(options.Figure_id);clf;set(gcf,'color','white', 'visible', 'on');
            end

            disp_wave = unwrap(angle(iq_wave));
            disp_wave_heart = bandpass(disp_wave, [obj.ecg_freqmin obj.ecg_freqmax], obj.fps);

            disp_resp = obj.lmd/(4*pi)*disp_wave * 1000; %
            disp_heart = obj.lmd/(4*pi)*disp_wave_heart * 1000;

            x_time = (0:numel(iq_wave)-1) * obj.dt;
            LineWidth=1;
            yyaxis left
            plot(x_time, disp_resp, 'LineWidth',LineWidth, 'Color',"#0072BD", 'DisplayName','Respiration');%hold on
            xlabel('Time (s)');
            ylabel('Displacement (mm)');
            yyaxis right
            plot(x_time, disp_heart, 'LineWidth',LineWidth, 'Color',"#D95319",  'DisplayName','Heart rate');
            xlim([x_time(1) x_time(end)])

            legend('location', 'bestoutside');

            if isfield(options, "Title")
                title(options.Title);
            end

            txt = sprintf('(x, y)=(%d, %d)', id_x, id_y);
            text(0.95, 1.03, txt, 'Units','normalized')
        end

        %%
        function showIQMovie(obj, id_x, id_y, radar_data, drawing_speed, options)
            arguments
                obj
                id_x
                id_y
                radar_data
                drawing_speed
                options.StartTime=0
                options.IntervalTime=size(radar_data,1)/obj.fps
                options.SaveMovieFn
                options.FigureID = 2
            end


            start_time_index = ceil(options.StartTime * obj.fps) + 1;
            end_time_index = start_time_index + ceil(options.IntervalTime * obj.fps) - 1;
            save_movie_flg = isfield(options, "SaveMovieFn");
            fprintf('(Start) %.2fsec - %.2fsec(End)\n', start_time_index/obj.fps, end_time_index/obj.fps);
            if save_movie_flg % ビデオの記録
                fprintf('SaveMovieMode On\n');
                v = VideoWriter(options.SaveMovieFn,'MPEG-4');
                cleanupObj = onCleanup(@() obj.cleanupFunction(v));
                open(v);
            end

            iq_wave = squeeze(radar_data(:, id_x, id_y));
            len_iq_wave = length(iq_wave);
            iq_wave_r = real(iq_wave);
            iq_wave_i = imag(iq_wave);
            disp_wave = unwrap(angle(iq_wave));
            disp_resp = obj.lmd/(4*pi)*disp_wave * 1000;


            figure(options.FigureID);clf(gcf);set(gcf, 'visible', 'on');
            
            subplot(2,1,1);          
            i = 1;
            p = plot(iq_wave_r(1:i), iq_wave_i(1:i), 'Color','#CCCCCC'); hold on;
            point_zero = plot(0,0,'.b','MarkerSize',15); hold on;
            m = scatter(iq_wave_r(i),iq_wave_i(i),15, 'r', 'filled'); hold off;
            daspect([1 1 1]);

            subplot(2,1,2);
            axObj = gca;
            x_time = (0:i-1) * obj.dt;
            p2 = plot(x_time, disp_resp(1:i));
            xlim([options.StartTime, options.StartTime+options.IntervalTime]);
            xlabel('Time (sec)');
            ylabel('Displacement (mm)');

            for ti = start_time_index:drawing_speed:end_time_index
                p.XData = iq_wave_r(1:ti);
                p.YData = iq_wave_i(1:ti);
                m.XData = iq_wave_r(ti);
                m.YData = iq_wave_i(ti);

                p2.XData = (0:ti-1) * obj.dt;
                p2.YData = disp_resp(1:ti);
                % refreshdata(surfObj, 'caller')
                title(axObj, sprintf('%.1f (s) [ID %d]', ti/obj.fps, ti));
                drawnow

                % pause(0.1)

                % Saving the figure
                if save_movie_flg % ビデオの記録
                    frame = getframe(gcf);
                    writeVideo(v,frame);
                end
            end

            if save_movie_flg % ビデオの記録
                close(v);
            end

        end

        function [u0, v0] = showIQImage(obj, id_x, id_y, radar_data, options)
            arguments (Input)
                obj
                id_x
                id_y
                radar_data
                options.FigureMode=false
                options.Figure_id=91
                options.Title
            end
            iq_wave = squeeze(radar_data(:, id_x, id_y));
            len_iq_wave = length(iq_wave);
            iq_wave_r = real(iq_wave);
            iq_wave_i = imag(iq_wave);
            disp_wave = unwrap(angle(iq_wave));
            disp_resp = obj.lmd/(4*pi)*disp_wave * 1000;

            % 中心推定
            func = @(x)sum( ( sqrt((iq_wave_r-x(1)).^2 + (iq_wave_i-x(2)).^2) - x(3) ).^2 );
            [x,fval] = fminunc(func, [-1,-1,-1]);
            u0 = x(1); v0 = x(2);
            if options.FigureMode
                figure(options.Figure_id);clf(gcf);set(gcf, 'visible', 'on');
                if isfield(options, "Title")
                    title(options.Title);
                end
                subplot(2,1,1);

                p = plot(iq_wave_r, iq_wave_i, 'Color','#CCCCCC'); hold on;

                scatter(u0, v0, 15, 'r', 'filled'); hold on;
                point_zero = plot(0,0,'.b','MarkerSize',15); hold off;
                daspect([1 1 1]);
                xlabel('I');
                ylabel('Q');

                subplot(2,1,2);
                x_time = (0:len_iq_wave-1) * obj.dt;
                p2 = plot(x_time, disp_resp);
                xlim([0 len_iq_wave*obj.dt]);
                xlabel('Time (sec)');
                ylabel('Displacement (mm)');
            end
        end

        function [u0, v0] = showIndexDisp(obj, id_x, id_y, radar_data, options)
            arguments (Input)
                obj
                id_x
                id_y
                radar_data
                options.FigureMode=false
                options.Figure_id=91
                options.Title
            end
            iq_wave = squeeze(radar_data(:, id_x, id_y));
            len_iq_wave = length(iq_wave);
            iq_wave_r = real(iq_wave);
            iq_wave_i = imag(iq_wave);
            disp_wave = unwrap(angle(iq_wave));
            disp_resp = obj.lmd/(4*pi)*disp_wave * 1000;

            % 中心推定
            func = @(x)sum( ( sqrt((iq_wave_r-x(1)).^2 + (iq_wave_i-x(2)).^2) - x(3) ).^2 );
            [x,fval] = fminunc(func, [-1,-1,-1]);
            u0 = x(1); v0 = x(2);


            figure(options.Figure_id);clf(gcf);set(gcf, 'visible', 'on');
            tiledlayout(2,1);
            nexttile;

            p = plot(iq_wave_r, iq_wave_i, 'Color','#CCCCCC'); hold on;

            scatter(u0, v0, 15, 'r', 'filled'); hold on;
            point_zero = plot(0,0,'.b','MarkerSize',15); hold off;
            daspect([1 1 1]);
            xlabel('I');
            ylabel('Q');

            subplot(2,1,2);
            x_time = (0:len_iq_wave-1) * obj.dt;
            p2 = plot(x_time, disp_resp);
            xlim([0 len_iq_wave*obj.dt]);
            xlabel('Time (sec)');
            ylabel('Displacement (mm)');

        end

        % Doppler of Single IF Signal
        function [Time,Freq,Y_stft,P_stft] = showDoppler(obj, IFSignal, options)
            arguments
                obj 
                IFSignal (:,1) % 1D-complex signal
                options.windowName
                options.Nfft=1024
                options.windowWidth=0.5 % 0.5 s
                options.FigureMode=true
                options.FigureID=10
            end
            [Time,Freq,Y_stft,P_stft] = my_stft(IFSignal, obj.fps, options.Nfft, options.windowWidth, options.windowName, isFFTShift=true);

            if options.FigureMode
                skip_slowt = 10;
                skip_freq = 10;
                Time_skip = Time(1:skip_slowt:end, 1:skip_freq:end);
                 Freq_skip = Freq(1:skip_slowt:end, 1:skip_freq:end);
                  P_stft_skip = P_stft(1:skip_slowt:end, 1:skip_freq:end);
                figure(options.FigureID);clf(gcf);
                % surface(Time,Freq,P_stft,P_stft,'EdgeColor','none');
                pcolorObj = pcolor(Time_skip, Freq_skip, P_stft_skip);
                pcolorObj.EdgeColor = 'None';

                colorbar;

                ylabel("Freq. (Hz)");xlabel("Time (s)")
            end

  
        end


        %% Utils
        % slice
        function out_radar_data = sliceTime(obj, radar_data, options)
            arguments (Input)
                obj
                radar_data (:,1) {mustBeNumeric}
                options.StartTime=0
                options.IntervalTime=length(radar_data)
            end
            start_time_index = ceil(options.StartTime * obj.fps) + 1;
            end_time_index = start_time_index + ceil(options.IntervalTime * obj.fps) - 1;

            start_time_index = max(start_time_index, 1);
            end_time_index = min(end_time_index, length(radar_data));
            out_radar_data = radar_data(start_time_index:end_time_index);
        end
        % figure x time
        function time_index = x_time(obj, length_signal, options)
            arguments (Input)
                obj
                length_signal % size n×1
                options.StartTime=0
            end
            time_index = (0 : length_signal-1)/obj.fps + options.StartTime;
        end
        % 
        function [M,I] = max_2d(~, a)
            % find max value of 2d matrix
            % Input
            %   a: matrix
            % Output
            %   M: max value
            %   I: [row,col] index
            if ~(length(size(a))==2)
                error('The dimension of the matrix must be 2.');
            end
            [tmp,aa] = max(a,[],1);
            [M,bb] = max(tmp);
            I = [aa(bb),bb];
        end

        %%% ユーザーが中断した場合に実行されるクリーンアップ関数
        function cleanupFunction(~, videoReader)
            disp('関数が終了（または強制終了）したよ');
            close(videoReader);
        end

        function radar_data = radardataconvert(obj, filename, matrix_size, options)
            % matrix_size is based on config files. before use this function, please load config
            arguments
                obj
                filename {mustBeText}
                matrix_size
                options.saveconvertfiles(1,1) {mustBeNumericOrLogical} = 0
            end

            % open files
            fileID = fopen(filename, 'r');
            % error check
            if fileID == -1
                error('Cannot open file: %s', filename);
            end
            % loading files
            data = fread(fileID, prod(matrix_size), '*uint8');

            % C language -- Row major same as numpy, but matlab load as column major
            % then load below
            matrix_size = matrix_size(end:-1:1);
            % change matrix size
            reshaped_data = reshape(data, matrix_size);
            % swap dimension
            reshaped_data = permute(reshaped_data, [7 6 5 4 3 2 1]);
            % change type
            reshaped_data = uint16(reshaped_data);

            I_data = (reshaped_data(:,:,:,:,:,1,1) + reshaped_data(:,:,:,:,:,1,2)*256);
            I_data = obj.uint16ToInt16(I_data);

            Q_data = (reshaped_data(:,:,:,:,:,2,1) + reshaped_data(:,:,:,:,:,2,2)*256);
            Q_data = obj.uint16ToInt16(Q_data);

            I_data(I_data >= 0x8000) = I_data(I_data >= 0x8000) - int16(0x10000);
            Q_data(Q_data >= 0x8000) = Q_data(Q_data >= 0x8000) - int16(0x10000);

            % make complex signal
            radar_data =  double(complex(I_data, Q_data));

            % close file
            fclose(fileID);

            if options.saveconvertfiles ~= 0
                % Use fileparts to get path, name, and extension
                [pathstr, name, ~] = fileparts(filename);
                filename_save = fullfile(pathstr, name);
                save(filename_save + ".mat", 'radar_data');
            end
        end

        % numpy used a wraparound method for data conversion when overflow occurred, while matlab filled in with the maximum value
        % Functions for type conversion in wrap-around
        function output = uint16ToInt16(~, input)
            tmp = int32(input);
            output = int16(input);

            % Mapping the upper half of uint16 (32768 - 65535) to the negative range of int16 (-32768 - -1)
            output(input > intmax('int16')) = int16(tmp(input > intmax('int16')) - 65536);
        end
    end
end





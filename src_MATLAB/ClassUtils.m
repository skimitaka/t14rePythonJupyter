classdef ClassUtils
    methods
        function [P1,f]=showFFT(obj, signal, fps, options)
            arguments (Input)
                obj
                signal
                fps
                options.Figure_id=10
                options.XLim
                options.LogMode=false
                options.Title='FFT'
                options.Subplot=false
                options.WindowFunction (1,1) string {mustBeMember(options.WindowFunction, ["rectangular", "hamming", "hann", "kaiser"])} = "rectangular"
            end

            len_signal = length(signal);

            switch options.WindowFunction
                case "rectangular"
                    signal_windowed = signal;
                case "hann"
                    signal_windowed = signal.*hann(len_signal, 'periodic');
                case "hamming"
                    signal_windowed = signal.*hamming(len_signal, 'periodic');
                case "kaiser"
                    signal_windowed = signal.*kaiser(len_signal, 30);
            end

            fft_signal = fft(signal_windowed);
            P2 = abs(fft_signal/len_signal);
            P1 = P2(1:len_signal/2+1);
            P1(2:end-1) = 2*P1(2:end-1);
            if options.LogMode
                P1 = 20*log10(P1);
            end
            f = fps*(0:(len_signal/2))/len_signal;

            if options.Subplot
                figure(gcf);
                nexttile;
            else
                figure(options.Figure_id);set(gcf,'color','white', 'visible', 'on');
                clf(gcf);
            end
            plot(f,P1, 'k', 'LineWidth',1);
            title(options.Title);
            if options.LogMode
                ylabel('Amplitude (dB)');
            else
                ylabel('Amplitude (μV)');
            end
            xlabel('Frequency (Hz)');

            % disp_mean=16.6;disp_sigma=2.8;
            % xregion((disp_mean-disp_sigma*3)/60,(disp_mean+disp_sigma*3)/60, FaceColor="#0072BD", FaceAlpha=0.1);
            % xregion(60/60, 100/60, FaceColor="#D95319", FaceAlpha=0.1);

            if isfield(options, "XLim")
                xlim(options.XLim);
            end
        end

        function showPSD(obj, signal, fps, options)
            arguments (Input)
                obj
                signal
                fps
                options.Figure_id=191
                options.XLim=[0 5]
                options.LogMode=true
                options.Normalized=true
                options.smooth (1,1) {mustBeNumeric}  = 0
                options.Subplot=false
                options.TwoSidedMode = false
                options.WindowFunction (1,1) string {mustBeMember(options.WindowFunction, ["rectangular", "hamming", "hann", "kaiser", "tailor"])} = "tailor"
            end

            len_signal = length(signal);
            switch options.WindowFunction
                case "rectangular"
                    signal_windowed = signal;
                case "hann"
                    signal_windowed = signal.*hann(len_signal, 'periodic');
                case "hamming"
                    signal_windowed = signal.*hamming(len_signal, 'periodic');
                case "kaiser"
                    signal_windowed = signal.*kaiser(len_signal, 10);
                case "tailor"
                    signal_windowed = signal.*taylorwin(len_signal);
            end

            fft_signal = fft(signal_windowed);
            if options.TwoSidedMode
                freq = (-len_signal/2:len_signal/2-1)/len_signal*fps;
                freq = freq(:);
                fft_signal = fftshift(fft_signal);
                psdx = (1/(fps*len_signal)) * abs(fft_signal).^2;

            else
                freq = fps*(0:(len_signal/2))/len_signal;
                xdft = fft_signal(1:len_signal/2+1);
                psdx = (1/(fps*len_signal)) * abs(xdft).^2;
                psdx(2:end-1) = 2*psdx(2:end-1);
            end

            if options.Normalized
                psdx = psdx./max(psdx);
            end
            if options.LogMode
                psdx = pow2db(psdx); % =10*log10(psdx);
            end
            if options.smooth ~= 0
                psdx = smooth(psdx, options.smooth/mean(diff(freq)));
            end

            if options.Subplot
                figure(gcf);
                nexttile;
            else
                figure(options.Figure_id);set(gcf,'color','white', 'visible', 'on');
                clf(gcf);
            end
            plot(freq, psdx, 'k', 'LineWidth',1);

            xlabel('Frequency (Hz)');
            if options.LogMode
                ylabel('Amplitude (dB)');
            else
                ylabel('Amplitude (μV)');
            end
            if isfield(options, "XLim")
                xlim(options.XLim);
            end
        end

        %% Figure
        function figure_obj = Figure(~, options)
            arguments (Input)
                ~
                options.ID
            end

            if isfield(options, "ID")
                figure_obj = figure(options.ID);clf(gcf);set(gcf, 'visible', 'on');
            else
                figure_obj = figure('Name','tmp');clf(gcf);set(gcf, 'visible', 'on');
            end
        end
    end
end
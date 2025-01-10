% CRISTIAN HERNANDEZ
% EE 399 CAPSTONE PROJECT
% NORTHWESTERN UNIVERSITY

% MATLAB Simulation for Acoustic Roughness Analysis

% Parameters
fs = 44100; % Sampling rate
duration = 1; % Duration in seconds
t = linspace(0, duration, fs * duration);
fc = 1000; % Center frequency
fm = 70; % Modulation frequency

% Generate initial pAM tone
phi = 0;
initial_pAM = generate_pAM_tone(phi, t, fc, fm);

% Parameters for the filter bank
numFilters = 10;
lowFreq = 100;
highFreq = 8000;
centerFreqs = logspace(log10(lowFreq), log10(highFreq), numFilters);

% Design bandpass filters
filters = cell(numFilters, 1);
for k = 1:numFilters
    f_low = centerFreqs(k) / sqrt(2);
    f_high = centerFreqs(k) * sqrt(2);
    [b, a] = butter(2, [f_low, f_high] / (fs / 2), 'bandpass');
    filters{k} = {b, a};
end

% Apply the filter bank to the initial pAM tone
filteredSignals = zeros(numFilters, length(initial_pAM));
for k = 1:numFilters
    [b, a] = filters{k}{:};
    filteredSignals(k, :) = filter(b, a, initial_pAM);
end

% Calculate roughness for each filtered signal
roughness_values = zeros(1, numFilters);
for k = 1:numFilters
    roughness_values(k) = calculate_roughness(filteredSignals(k, :));
end

% Total roughness
total_roughness = sum(roughness_values);
fprintf('Estimated Roughness for initial phase %f: %f\n', phi, total_roughness);

% Play the initial pAM tone
sound(initial_pAM, fs);

% Create a figure window
hFig = figure('Name', 'Acoustic Roughness Simulation', 'NumberTitle', 'off', ...
              'Position', [100, 100, 1000, 700], 'Color', [0.94, 0.94, 0.94]);

% Create panels
controlPanel = uipanel('Title', 'Controls', 'FontSize', 12, ...
                       'BackgroundColor', 'white', 'Position', [0.05, 0.05, 0.9, 0.25]);

plotPanel = uipanel('Title', 'Plots', 'FontSize', 12, ...
                    'BackgroundColor', 'white', 'Position', [0.05, 0.35, 0.9, 0.6]);

% Create sliders for phase, modulation frequency, and duration within the control panel
phaseSlider = uicontrol('Parent', controlPanel, 'Style', 'slider', 'Min', -pi, 'Max', pi, ...
                        'Value', 0, 'Position', [150, 70, 400, 20]);

fmSlider = uicontrol('Parent', controlPanel, 'Style', 'slider', 'Min', 20, 'Max', 200, ...
                     'Value', 70, 'Position', [150, 40, 400, 20]);

durationSlider = uicontrol('Parent', controlPanel, 'Style', 'slider', 'Min', 0.5, 'Max', 5, ...
                           'Value', 1, 'Position', [150, 10, 400, 20]);

% Create text labels for sliders within the control panel
uicontrol('Parent', controlPanel, 'Style', 'text', 'Position', [50, 65, 80, 30], ...
          'String', 'Phase', 'FontSize', 12, 'BackgroundColor', 'white');

uicontrol('Parent', controlPanel, 'Style', 'text', 'Position', [50, 35, 80, 30], ...
          'String', 'Mod Freq', 'FontSize', 12, 'BackgroundColor', 'white');

uicontrol('Parent', controlPanel, 'Style', 'text', 'Position', [50, 5, 80, 30], ...
          'String', 'Duration', 'FontSize', 12, 'BackgroundColor', 'white');

% Create dynamic text to display the current values of the sliders
phaseLabel = uicontrol('Parent', controlPanel, 'Style', 'text', 'Position', [580, 65, 50, 30], ...
                       'String', '0', 'FontSize', 12, 'BackgroundColor', 'white');

fmLabel = uicontrol('Parent', controlPanel, 'Style', 'text', 'Position', [580, 35, 50, 30], ...
                    'String', '70', 'FontSize', 12, 'BackgroundColor', 'white');

durationLabel = uicontrol('Parent', controlPanel, 'Style', 'text', 'Position', [580, 5, 50, 30], ...
                          'String', '1', 'FontSize', 12, 'BackgroundColor', 'white');

% Plot area within the plot panel
waveformAxes = axes('Parent', plotPanel, 'Position', [0.1, 0.55, 0.35, 0.4]);
waveformPlot = plot(waveformAxes, t, initial_pAM);
title(waveformAxes, 'Waveform');
xlabel(waveformAxes, 'Time (s)');
ylabel(waveformAxes, 'Amplitude');

roughnessAxes = axes('Parent', plotPanel, 'Position', [0.55, 0.55, 0.35, 0.4]);
roughnessPlot = bar(roughnessAxes, roughness_values, 'r');
title(roughnessAxes, 'Roughness');
xlabel(roughnessAxes, 'Filter Channel');
ylabel(roughnessAxes, 'Roughness');

spectrogramAxes = axes('Parent', plotPanel, 'Position', [0.1, 0.1, 0.8, 0.35]);
spectrogram(initial_pAM, 256, [], [], fs, 'yaxis');
title(spectrogramAxes, 'Spectrogram');

% Store UI elements and parameters in guidata
handles = struct('phaseSlider', phaseSlider, 'fmSlider', fmSlider, 'durationSlider', durationSlider, ...
                 'waveformPlot', waveformPlot, 'roughnessPlot', roughnessPlot, 'spectrogramAxes', spectrogramAxes, ...
                 't', t, 'fc', fc, 'fs', fs, 'filters', {filters}, 'numFilters', numFilters, ...
                 'phaseLabel', phaseLabel, 'fmLabel', fmLabel, 'durationLabel', durationLabel);
guidata(hFig, handles);

% Add listeners to the sliders
addlistener(phaseSlider, 'ContinuousValueChange', @(src, evt)updateSimulation(hFig));
addlistener(fmSlider, 'ContinuousValueChange', @(src, evt)updateSimulation(hFig));
addlistener(durationSlider, 'ContinuousValueChange', @(src, evt)updateSimulation(hFig));

% Update function
function updateSimulation(hFig)
    handles = guidata(hFig);
    
    % Get slider values
    phase = get(handles.phaseSlider, 'Value');
    fm = get(handles.fmSlider, 'Value');
    duration = get(handles.durationSlider, 'Value');
    
    % Update slider labels
    set(handles.phaseLabel, 'String', num2str(phase, '%.2f'));
    set(handles.fmLabel, 'String', num2str(fm, '%.0f'));
    set(handles.durationLabel, 'String', num2str(duration, '%.2f'));
    
    % Update time vector based on new duration
    handles.t = linspace(0, duration, handles.fs * duration);
    
    % Generate pAM tone with new parameters
    new_pAM = generate_pAM_tone(phase, handles.t, handles.fc, fm);
    
    % Play the new pAM tone
    sound(new_pAM, handles.fs);
    
    % Update filtered signals
    filteredSignals = zeros(handles.numFilters, length(new_pAM));
    for k = 1:handles.numFilters
        [b, a] = handles.filters{k}{:};
        filteredSignals(k, :) = filter(b, a, new_pAM);
    end
    
    % Update roughness values
    roughness_values = zeros(1, handles.numFilters);
    for k = 1:handles.numFilters
        roughness_values(k) = calculate_roughness(filteredSignals(k, :));
    end
    total_roughness = sum(roughness_values);
    
    % Update plots
    set(handles.waveformPlot, 'YData', new_pAM, 'XData', handles.t);
    set(handles.roughnessPlot, 'YData', roughness_values);
    
    % Plot spectrogram on spectrogramAxes
    cla(handles.spectrogramAxes); % Clear previous spectrogram
    axes(handles.spectrogramAxes); % Set current axes to spectrogramAxes
    spectrogram(new_pAM, 256, [], [], handles.fs, 'yaxis');
    title(handles.spectrogramAxes, 'Spectrogram');
    
    fprintf('Estimated Roughness for phase %f, mod freq %f, and duration %f: %f\n', phase, fm, duration, total_roughness);
end

% Define the pAM tone generation function at the end
function tone = generate_pAM_tone(phi, t, fc, fm)
    tone = 0.5 * cos(2 * pi * (fc - fm) * t) + ...
           cos(2 * pi * fc * t + phi) + ...
           0.5 * cos(2 * pi * (fc + fm) * t);
end

% Calculate the roughness of a signal
function roughness = calculate_roughness(signal)
    roughness = rms(abs(hilbert(signal)));
end

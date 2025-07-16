clear all; close all; clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script Name : SpCord_retro_Plots
% Author      : Shreyas Suryanarayana
% Affiliation : Duke University School of Medicine
% Email       : shreyas.suryanarayana@duke.edu
% Date        : 07/14/25
%
% Description:
% Plotting and statistical testing of GFP and RFP cell counts from lightsheet imaging data. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Specify the file path
disp('Select Br-1/Br-2 GFP Excel File');
[filename1, filepath1] = uigetfile('*.xlsx', 'Select GFP Excel File');
% disp('Select Br-1 RFP Excel File');
% [filename2, filepath2] = uigetfile('*.csv', 'Select RFP Excel File');
disp('Select Br-1/Br-2 RFP Excel File');
[filename2, filepath2] = uigetfile('*.xlsx', 'Select GFP Excel File');
% disp('Select Br-2 RFP Excel File');
% [filename5, filepath5] = uigetfile('*.csv', 'Select GFP Excel File');


fullFileName1 = fullfile(filepath1, filename1);
fullFileName2 = fullfile(filepath2, filename2);


disp(['Selected file: ', fullFileName1]);
disp(['Selected file: ', fullFileName2]);
%%
% Read the data into a table
dataTable1 = readtable(filename1);
dataTable2 = readtable(filename2);
dataTable3 = readtable('/Users/sm824/Desktop/Anatomy Analysis/Retrograde Analysis/AreaIDS_Retro');
dataTable4 = readtable('/Users/sm824/Desktop/Anatomy Analysis/Retrograde Analysis/AreaIDS_Lev1_Retro');
dataTable5 = readtable('/Users/sm824/Desktop/Anatomy Analysis/Retrograde Analysis/AreaIDS_Lev2_Retro');
dataTable6 = readtable('/Users/sm824/Desktop/Anatomy Analysis/Retrograde Analysis/SpCorP04/Isocortex_Comp_Areas');
dataTable7 = readtable('/Users/sm824/Desktop/Anatomy Analysis/Retrograde Analysis/AreaIDS_Lev3_Retro');
%% convert tables to arrays
dataTable_GFP(:,[1,2,3,4]) = table2array(dataTable1(:,[1,4,5,6]));
dataTable_RFP(:,[1,2,3,4]) = table2array(dataTable2(:,[1,4,5,6]));


% extracting selected areas and their names %
dataTable_Areas = table2array(dataTable3(:,1));
dataTable_Area_names = cellfun(@(x) x(1:(length(x)-2)), table2array(dataTable3(:,3)), 'UniformOutput', false);

dataTable_Areas_Lev1 = table2array(dataTable4(:,1));
dataTable_Area_Lev1_names = cellfun(@(x) x(1:(length(x)-2)), table2array(dataTable4(:,3)), 'UniformOutput', false);

dataTable_Areas_comp = table2array(dataTable6(:,1));
dataTable_Areas_comp_names = cellfun(@(x) x(1:(length(x)-2)), table2array(dataTable6(:,3)), 'UniformOutput', false);

dataTable_Areas_Lev3 = table2array(dataTable7(:,1));
dataTable_Area_Lev3_names = cellfun(@(x) x(1:(length(x)-2)), table2array(dataTable7(:,3)), 'UniformOutput', false);

% saving %
writecell(dataTable_Area_names, fullfile(filepath1, 'Area_Names.xlsx'));
writecell(dataTable_Area_Lev1_names, fullfile(filepath1, 'Area_Lev1_Names.xlsx'));

%% Getting GFP/RFP data for all selected areas
% ipsi
for i = 1:size((dataTable_Areas),1)
    for j = 1:size((dataTable_GFP),1)
        if dataTable_Areas(i,1) == dataTable_GFP(j,1)
            dataTable_sel_ipsi_GFP(i,1) = dataTable_GFP(j,1);
            dataTable_sel_ipsi_GFP(i,2) = dataTable_GFP(j,2);
            dataTable_sel_ipsi_GFP(i,3) = dataTable_GFP(j,3);
            dataTable_sel_ipsi_GFP(i,4) = dataTable_GFP(j,4);

        end
        if dataTable_Areas(i,1) == dataTable_RFP(j,1)
            dataTable_sel_ipsi_RFP(i,1) = dataTable_RFP(j,1);
            dataTable_sel_ipsi_RFP(i,2) = dataTable_RFP(j,2);
            dataTable_sel_ipsi_RFP(i,3) = dataTable_RFP(j,3);
            dataTable_sel_ipsi_RFP(i,4) = dataTable_RFP(j,4);
        end
    end
end

% contra %
for i = 1:size((dataTable_Areas),1)
    for j = 1:size((dataTable_GFP),1)
        if (dataTable_Areas(i,1) - 10000) == dataTable_GFP(j,1)
            dataTable_sel_contra_GFP(i,1) = dataTable_GFP(j,1);
            dataTable_sel_contra_GFP(i,2) = dataTable_GFP(j,2);
            dataTable_sel_contra_GFP(i,3) = dataTable_GFP(j,3);
            dataTable_sel_contra_GFP(i,4) = dataTable_GFP(j,4);
        end
        if (dataTable_Areas(i,1) - 10000) == dataTable_RFP(j,1)
            dataTable_sel_contra_RFP(i,1) = dataTable_RFP(j,1);
            dataTable_sel_contra_RFP(i,2) = dataTable_RFP(j,2);
            dataTable_sel_contra_RFP(i,3) = dataTable_RFP(j,3);
            dataTable_sel_contra_RFP(i,4) = dataTable_RFP(j,4);
        end
    end
end
%% Getting GFP/RFP data for Level 1 selected areas
% ipsi
for i = 1:size((dataTable_Areas_Lev1),1)
    for j = 1:size((dataTable_GFP),1)
        if dataTable_Areas_Lev1(i,1) == dataTable_GFP(j,1)
            dataTable_Lev1_ipsi_GFP(i,1) = dataTable_GFP(j,1);
            dataTable_Lev1_ipsi_GFP(i,2) = dataTable_GFP(j,2);
            dataTable_Lev1_ipsi_GFP(i,3) = dataTable_GFP(j,3);
            dataTable_Lev1_ipsi_GFP(i,4) = dataTable_GFP(j,4);
        end
        if dataTable_Areas_Lev1(i,1) == dataTable_RFP(j,1)
            dataTable_Lev1_ipsi_RFP(i,1) = dataTable_RFP(j,1);
            dataTable_Lev1_ipsi_RFP(i,2) = dataTable_RFP(j,2);
            dataTable_Lev1_ipsi_RFP(i,3) = dataTable_RFP(j,3);
            dataTable_Lev1_ipsi_RFP(i,4) = dataTable_RFP(j,4);
        end
    end
end

% contra %
for i = 1:size((dataTable_Areas_Lev1),1)
    for j = 1:size((dataTable_GFP),1)
        if (dataTable_Areas_Lev1(i,1) - 10000) == dataTable_GFP(j,1)
            dataTable_Lev1_contra_GFP(i,1) = dataTable_GFP(j,1);
            dataTable_Lev1_contra_GFP(i,2) = dataTable_GFP(j,2);
            dataTable_Lev1_contra_GFP(i,3) = dataTable_GFP(j,3);
            dataTable_Lev1_contra_GFP(i,4) = dataTable_GFP(j,4);
        end
        if (dataTable_Areas_Lev1(i,1) - 10000) == dataTable_RFP(j,1)
            dataTable_Lev1_contra_RFP(i,1) = dataTable_RFP(j,1);
            dataTable_Lev1_contra_RFP(i,2) = dataTable_RFP(j,2);
            dataTable_Lev1_contra_RFP(i,3) = dataTable_RFP(j,3);
            dataTable_Lev1_contra_RFP(i,4) = dataTable_RFP(j,4);
        end
    end
end
%% Getting GFP/RFP data for Level 3 selected areas
% ipsi
for i = 1:size((dataTable_Areas_Lev3),1)
    for j = 1:size((dataTable_GFP),1)
        if dataTable_Areas_Lev3(i,1) == dataTable_GFP(j,1)
            dataTable_Lev3_ipsi_GFP(i,1) = dataTable_GFP(j,1);
            dataTable_Lev3_ipsi_GFP(i,2) = dataTable_GFP(j,2);
            dataTable_Lev3_ipsi_GFP(i,3) = dataTable_GFP(j,3);
            dataTable_Lev3_ipsi_GFP(i,4) = dataTable_GFP(j,4);
        end
        if dataTable_Areas_Lev3(i,1) == dataTable_RFP(j,1)
            dataTable_Lev3_ipsi_RFP(i,1) = dataTable_RFP(j,1);
            dataTable_Lev3_ipsi_RFP(i,2) = dataTable_RFP(j,2);
            dataTable_Lev3_ipsi_RFP(i,3) = dataTable_RFP(j,3);
            dataTable_Lev3_ipsi_RFP(i,4) = dataTable_RFP(j,4);
        end
    end
end

% contra %
for i = 1:size((dataTable_Areas_Lev3),1)
    for j = 1:size((dataTable_GFP),1)
        if (dataTable_Areas_Lev3(i,1) - 10000) == dataTable_GFP(j,1)
            dataTable_Lev3_contra_GFP(i,1) = dataTable_GFP(j,1);
            dataTable_Lev3_contra_GFP(i,2) = dataTable_GFP(j,2);
            dataTable_Lev3_contra_GFP(i,3) = dataTable_GFP(j,3);
            dataTable_Lev3_contra_GFP(i,4) = dataTable_GFP(j,4);
        end
        if (dataTable_Areas_Lev3(i,1) - 10000) == dataTable_RFP(j,1)
            dataTable_Lev3_contra_RFP(i,1) = dataTable_RFP(j,1);
            dataTable_Lev3_contra_RFP(i,2) = dataTable_RFP(j,2);
            dataTable_Lev3_contra_RFP(i,3) = dataTable_RFP(j,3);
            dataTable_Lev3_contra_RFP(i,4) = dataTable_RFP(j,4);
        end
    end
end
%% Calculate mean %%
ipsi_Lev1_GFP_mean = mean(dataTable_Lev1_ipsi_GFP(:,[2,3,4]), 2);
contra_Lev1_GFP_mean = mean(dataTable_Lev1_contra_GFP(:,[2,3,4]), 2);

ipsi_Lev3_GFP_mean = mean(dataTable_Lev3_ipsi_GFP(:,[2,3,4]), 2);
contra_Lev3_GFP_mean = mean(dataTable_Lev3_contra_GFP(:,[2,3,4]), 2);

ipsi_Lev1_RFP_mean = mean(dataTable_Lev1_ipsi_RFP(:,[2,3,4]), 2);
contra_Lev1_RFP_mean = mean(dataTable_Lev1_contra_RFP(:,[2,3,4]), 2);

ipsi_Lev3_RFP_mean = mean(dataTable_Lev3_ipsi_RFP(:,[2,3,4]), 2);
contra_Lev3_RFP_mean = mean(dataTable_Lev3_contra_RFP(:,[2,3,4]), 2);
%% Create the horizontal bar graph ipsi-contra %%
% Create the horizontal bar graph for ipsilateral and contralateral data
figure('Units', 'inches', 'Position', [2, 2, 50/25.4, 200/25.4]); % Set figure size to 100mm x 150mm (converted to inches)
hold on; % Hold on to the current axes

% Define y-axis positions for ipsilateral
y1_ipsi = [size(ipsi_Lev1_GFP_mean, 1):-1:1]; % Main y positions
y2_ipsi = y1_ipsi + 0.4; % Offset positions

% Bar plots for ipsilateral
barh(y1_ipsi, ipsi_Lev1_RFP_mean, 0.4, 'FaceColor', 'y', 'EdgeColor', 'y', 'DisplayName', 'Ipsi RFP Mean'); % Ipsi RFP
barh(y2_ipsi, ipsi_Lev1_GFP_mean, 0.4, 'FaceColor', 'm', 'EdgeColor', 'm', 'DisplayName', 'Ipsi GFP Mean'); % Ipsi GFP

% Overlay individual data points and connect with lines for ipsilateral
for i = 1:length(y1_ipsi)
    % RFP points and line
    x_rfp = dataTable_Lev1_ipsi_RFP(i, 2:4); % Extract values from columns 2, 3, and 4
    y_rfp = repmat(y1_ipsi(i), size(x_rfp)); % Ensure y-coordinates match x-coordinates
    plot(x_rfp, y_rfp, '-k', 'LineWidth', 0.5); % Connect points with a line
    scatter(x_rfp, y_rfp, 4, 'k', 'filled'); % Yellow points for RFP

    % GFP points and line
    x_gfp = dataTable_Lev1_ipsi_GFP(i, 2:4);
    y_gfp = repmat(y2_ipsi(i), size(x_gfp));
    plot(x_gfp, y_gfp, '-k', 'LineWidth', 0.5);
    scatter(x_gfp, y_gfp, 4, 'k', 'filled'); % Magenta points for GFP
end

% Define y-axis positions for contralateral
y1_contra = [size(contra_Lev1_GFP_mean, 1):-1:1]; % Main y positions
y2_contra = y1_contra + 0.4; % Offset positions

% Bar plots for contralateral
barh(y1_contra, -contra_Lev1_RFP_mean, 0.4, 'FaceColor', 'y', 'EdgeColor', 'y', 'DisplayName', 'Contra RFP Mean'); % Contra RFP
barh(y2_contra, -contra_Lev1_GFP_mean, 0.4, 'FaceColor', 'm', 'EdgeColor', 'm', 'DisplayName', 'Contra GFP Mean'); % Contra GFP

% Overlay individual data points and connect with lines for contralateral
for i = 1:length(y1_contra)
    % RFP points and line
    x_rfp = -dataTable_Lev1_contra_RFP(i, 2:4); % Extract values and negate them
    y_rfp = repmat(y1_contra(i), size(x_rfp));
    plot(x_rfp, y_rfp, '-k', 'LineWidth', 0.5);
    scatter(x_rfp, y_rfp, 4, 'k', 'filled'); % Yellow points for RFP

    % GFP points and line
    x_gfp = -dataTable_Lev1_contra_GFP(i, 2:4);
    y_gfp = repmat(y2_contra(i), size(x_gfp));
    plot(x_gfp, y_gfp, '-k', 'LineWidth', 0.5);
    scatter(x_gfp, y_gfp, 4, 'k', 'filled'); % Magenta points for GFP
end

% Add labels and title
ylabel('Isocortical Areas', 'FontSize', 5); % Set font size to 5
xlabel('Number of Labeled Neurons', 'FontName', 'Arial', 'FontWeight', 'bold', 'FontSize', 5); % Set font size to 5
title('Corticospinal P04 Contralateral/Ipsilateral Projections', 'FontName', 'Arial', 'FontWeight', 'bold', 'FontSize', 5); % Set font size to 5

% Customize y-ticks
yy = 1:size(contra_Lev1_GFP_mean);
yticks(yy);
yticklabels(flip(dataTable_Area_Lev1_names));

% Adjust axes properties
ax = gca; % Get the current axes
ax.YAxis.FontName = 'Arial'; % Set font to Arial
ax.YAxis.FontWeight = 'bold'; % Set font weight to bold
ax.YAxis.FontSize = 5; % Set font size to 5
ax.XAxis.FontName = 'Arial'; % Set font to Arial
ax.XAxis.FontWeight = 'bold'; % Set font weight to bold
ax.XAxis.FontSize = 5; % Set font size to 5

% Adjust x-ticks to show absolute values
xtick_values = xticks; % Get current x-tick values
xticklabels(arrayfun(@num2str, abs(xtick_values), 'UniformOutput', false)); % Make all labels positive

hold off;
% save figure %
saveas(gcf, fullfile(filepath1, 'Cor_spip4_ipsicontra'), 'epsc');

%% Create the horizontal bar graph ipsi-contra (log scale) %%
% Calculate mean %
ipsi_Lev1_GFP_mean = mean(dataTable_Lev1_ipsi_GFP(:,[2,3,4]), 2);
contra_Lev1_GFP_mean = mean(dataTable_Lev1_contra_GFP(:,[2,3,4]), 2);

ipsi_Lev1_RFP_mean = mean(dataTable_Lev1_ipsi_RFP(:,[2,3,4]), 2);
contra_Lev1_RFP_mean = mean(dataTable_Lev1_contra_RFP(:,[2,3,4]), 2);

% Calculate SEM across replicates (3 samples)
n_samples = 3;
ipsi_Lev1_GFP_sd = std(dataTable_Lev1_ipsi_GFP(:, [2,3,4]), 0, 2);
ipsi_Lev1_GFP_sem = ipsi_Lev1_GFP_sd / sqrt(n_samples);

contra_Lev1_GFP_sd = std(dataTable_Lev1_contra_GFP(:, [2,3,4]), 0, 2);
contra_Lev1_GFP_sem = contra_Lev1_GFP_sd / sqrt(n_samples);

ipsi_Lev1_RFP_sd = std(dataTable_Lev1_ipsi_RFP(:, [2,3,4]), 0, 2);
ipsi_Lev1_RFP_sem = ipsi_Lev1_RFP_sd / sqrt(n_samples);

contra_Lev1_RFP_sd = std(dataTable_Lev1_contra_RFP(:, [2,3,4]), 0, 2);
contra_Lev1_RFP_sem = contra_Lev1_RFP_sd / sqrt(n_samples);

% Create the horizontal bar graph ipsi-contra (log scale) %%
figure('Units', 'inches', 'Position', [2, 2, 50/25.4, 200/25.4]);
hold on;

% Define y-axis positions for ipsilateral
y1_ipsi = [size(ipsi_Lev1_GFP_mean, 1):-1:1];
y2_ipsi = y1_ipsi + 0.3;

% Apply log10 to ipsilateral means
ipsi_rfp_log = log10(ipsi_Lev1_RFP_mean);
ipsi_gfp_log = log10(ipsi_Lev1_GFP_mean);

% Propagate SEM to log10 scale: SEM_log = SEM / (mean * ln(10))
ipsi_rfp_sem_log = ipsi_Lev1_RFP_sem ./ (ipsi_Lev1_RFP_mean * log(10));
ipsi_gfp_sem_log = ipsi_Lev1_GFP_sem ./ (ipsi_Lev1_GFP_mean * log(10));

% Bar plots for ipsilateral
barh(y1_ipsi, ipsi_rfp_log, 0.3, 'FaceColor', 'y', 'EdgeColor', 'y');
barh(y2_ipsi, ipsi_gfp_log, 0.3, 'FaceColor', 'm', 'EdgeColor', 'm');

% Add SEM error bars for ipsilateral
errorbar(ipsi_rfp_log, y1_ipsi, ipsi_rfp_sem_log, 'horizontal', 'k', 'LineStyle', 'none', 'CapSize', 3);
errorbar(ipsi_gfp_log, y2_ipsi, ipsi_gfp_sem_log, 'horizontal', 'k', 'LineStyle', 'none', 'CapSize', 3);

% Define y-axis positions for contralateral
y1_contra = [size(contra_Lev1_GFP_mean, 1):-1:1];
y2_contra = y1_contra + 0.3;

% Apply log10 to contralateral means and negate
contra_rfp_log = -log10(contra_Lev1_RFP_mean);
contra_gfp_log = -log10(contra_Lev1_GFP_mean);

% Propagate SEM to log10 scale for contralateral
contra_rfp_sem_log = contra_Lev1_RFP_sem ./ (contra_Lev1_RFP_mean * log(10));
contra_gfp_sem_log = contra_Lev1_GFP_sem ./ (contra_Lev1_GFP_mean * log(10));

% Bar plots for contralateral
barh(y1_contra, contra_rfp_log, 0.3, 'FaceColor', 'y', 'EdgeColor', 'y');
barh(y2_contra, contra_gfp_log, 0.3, 'FaceColor', 'm', 'EdgeColor', 'm');

% Add SEM error bars for contralateral (negate errors for mirrored plot)
errorbar(contra_rfp_log, y1_contra, contra_rfp_sem_log, 'horizontal', 'k', 'LineStyle', 'none', 'CapSize', 3);
errorbar(contra_gfp_log, y2_contra, contra_gfp_sem_log, 'horizontal', 'k', 'LineStyle', 'none', 'CapSize', 3);

% Add labels and title
ylabel('Isocortical Areas', 'FontSize', 5);
xlabel('log_{10}(Number of Labeled Neurons)', 'FontName', 'Arial', 'FontWeight', 'bold', 'FontSize', 5);
title('Corticopontine Contralateral/Ipsilateral Projections (log scale)', 'FontName', 'Arial', 'FontWeight', 'bold', 'FontSize', 5);

% Customize y-ticks
yy = 1:size(contra_Lev1_GFP_mean);
yticks(yy);
yticklabels(flip(dataTable_Area_Lev1_names));

% Adjust axes properties
ax = gca;
ax.YAxis.FontName = 'Arial';
ax.YAxis.FontWeight = 'bold';
ax.YAxis.FontSize = 5;

ax.XAxis.FontName = 'Arial';
ax.XAxis.FontWeight = 'bold';
ax.XAxis.FontSize = 5;

% Set x-tick labels to show antilog values (approximate)
xticks_vals = xticks;
xticklabels(arrayfun(@(x) sprintf('%g', 10.^abs(x)), xticks_vals, 'UniformOutput', false));

hold off;

% Save figure
saveas(gcf, fullfile(filepath1, 'Cor_ipsicontra_log_sem'), 'epsc');



%% Correlation Analysis %%
% Data for Mouse 1
Ind1_Cond1 = [dataTable_sel_ipsi_GFP(:,2); dataTable_sel_contra_GFP(:,2)];  % Condition 1
Ind1_Cond2 = [dataTable_sel_ipsi_RFP(:,2); dataTable_sel_contra_RFP(:,2)];  % Condition 2

% Data for Mouse 2
Ind2_Cond1 = [dataTable_sel_ipsi_GFP(:,3); dataTable_sel_contra_GFP(:,3)];  % Condition 1
Ind2_Cond2 = [dataTable_sel_ipsi_RFP(:,3); dataTable_sel_contra_RFP(:,3)];  % Condition 2

% Data for Mouse 1
Ind3_Cond1 = [dataTable_sel_ipsi_GFP(:,4); dataTable_sel_contra_GFP(:,4)];  % Condition 1
Ind3_Cond2 = [dataTable_sel_ipsi_RFP(:,4); dataTable_sel_contra_RFP(:,4)];  % Condition 2
% Combine all data into a matrix (rows: observations, columns: conditions/individuals)
combined_data = [
    Ind1_Cond1(:), Ind1_Cond2(:), ...
    Ind2_Cond1(:), Ind2_Cond2(:), ...
    Ind3_Cond1(:), Ind3_Cond2(:), ...
    ];

% Compute the correlation matrix
corr_matrix = corrcoef(combined_data);

% Display the correlation matrix
disp('Direct and Indirect Correlation Matrix:');
disp(corr_matrix);

% Plot the correlation matrix
figure;
heatmap(corr_matrix, ...
    'XData', {'M1 Direct', 'M1 Indirect', 'M2 Direct', 'M2 Indirect', 'M3 Direct', 'M3 Indirect'}, ...
    'YData', {'M1 Direct', 'M1 Indirect', 'M2 Direct', 'M2 Indirect', 'M3 Direct', 'M3 Indirect'}, ...
    'Colormap', jet, 'ColorLimits', [-1, 1]);
title('Corticospinal P04 Direct and Indirect (across isocortex)');

% Save figure
saveas(gcf, fullfile(filepath1, 'Cort_spip4_Correlation_Mouse'), 'epsc');
%% total isocortex %%
% Data arrays
y1_GFP = [dataTable_sel_ipsi_GFP(4, [2, 3, 4]); dataTable_sel_contra_GFP(4, [2, 3, 4])];
y1_RFP = [dataTable_sel_ipsi_RFP(4, [2, 3, 4]); dataTable_sel_contra_RFP(4, [2, 3, 4])];

% Saving data
writematrix(y1_GFP, fullfile(filepath1, 'Cor_Spi_tot_GFP.xlsx'));
writematrix(y1_RFP, fullfile(filepath1, 'Cor_Spi_tot_RFP.xlsx'));

% Compute means and standard deviations
mean_GFP = mean(y1_GFP(:));
mean_RFP = mean(y1_RFP(:));
std_GFP = std(y1_GFP(:));
std_RFP = std(y1_RFP(:));

% Create figure
figure('Units', 'inches', 'Position', [2, 2, 6, 10]);
hold on;

% Bar heights
bar_data = [mean_GFP, mean_RFP];

% Bar colors
bar_colors = [1, 0, 1; 1, 1, 0];
b = bar([1, 2], bar_data, 'FaceColor', 'flat', 'BarWidth', 0.5);

% Assign colors
for k = 1:2
    b.CData(k, :) = bar_colors(k, :);
end

% Add vertical error bars
errorbar([1, 2], bar_data, [std_GFP, std_RFP], 'k.', 'CapSize', 10, 'LineStyle', 'none');

% Add data points
x_positions = [1, 2];
for i = 1:2
    if i == 1
        data_points = y1_GFP;
        color = [0, 0, 0];
    else
        data_points = y1_RFP;
        color = [0, 0, 0];
    end

    % Plot row 1 (filled)
    scatter(repmat(x_positions(i), 1, size(data_points, 2)), data_points(1, :), ...
        100, color, 'filled');

    % Plot row 2 (unfilled)
    scatter(repmat(x_positions(i), 1, size(data_points, 2)), data_points(2, :), ...
        100, 'o', 'MarkerEdgeColor', color);
end

% Statistical test
[h, p] = ttest2(y1_GFP(:), y1_RFP(:));

% Annotate plot
text(mean([x_positions]), max(bar_data) + 0.05 * max(bar_data), ...
    sprintf('p = %.10f', p), 'HorizontalAlignment', 'center', ...
    'FontSize', 12, 'FontWeight', 'bold');

% Title and labels
title('Corticospinal P04 neurons', 'FontName', 'Arial', 'FontWeight', 'bold', 'FontSize', 16);
set(gca, 'XTick', x_positions, 'XTickLabel', {'Indirect', 'Direct'}, 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Number of neurons', 'FontSize', 14, 'FontWeight', 'bold');

% Aesthetics
box on;
grid on;
ylim([0, max(bar_data) + 0.2 * max(bar_data)]);

hold off;

% Save figure
saveas(gcf, fullfile(filepath1, 'Cor_spip4_stat'), 'epsc');


%% Pie Chart %%
% Data arrays
y1_GFP = [dataTable_sel_ipsi_GFP(4, [2, 3, 4]); dataTable_sel_contra_GFP(4, [2, 3, 4])];
y1_RFP = [dataTable_sel_ipsi_RFP(4, [2, 3, 4]); dataTable_sel_contra_RFP(4, [2, 3, 4])];

% Compute sum of values across the rows for GFP and RFP
sum_GFP = sum(y1_GFP(:)); % Sum of all values in y1_GFP
sum_RFP = sum(y1_RFP(:)); % Sum of all values in y1_RFP

% Total sum for both groups
total_sum = sum_GFP + sum_RFP;

% Create pie chart
figure('Units', 'inches', 'Position', [2, 2, 6, 6]); % Adjust figure size
pie_data = [sum_GFP, sum_RFP];

% Labels with percentages
labels = {sprintf('Indirect (GFP)\n%.1f%%', (sum_GFP / total_sum) * 100), ...
    sprintf('Direct (RFP)\n%.1f%%', (sum_RFP / total_sum) * 100)};

% Pie chart with color
h = pie(pie_data, labels);
colormap([1, 0, 1; 1, 1, 0]); % Magenta for GFP, Yellow for RFP

% Set title for the pie chart
title_handle = title('Corticospinal P04 Neurons', 'FontName', 'Arial', 'FontWeight', 'bold', 'FontSize', 14);

% Adjust the title position (move it up)
title_handle.Position(2) = title_handle.Position(2) + 0.05; % Adjust the y-position slightly higher

% Set font to Arial and bold for all pie chart labels
for k = 1:length(h)
    if isprop(h(k), 'Text')  % Only set the font for text properties
        h(k).FontName = 'Arial';
        h(k).FontWeight = 'bold';
        h(k).FontSize = 12;
    end
end

% Display the pie chart
set(gca, 'FontSize', 12, 'FontWeight', 'bold');
% save figure %
saveas(gcf,fullfile(filepath1, 'Cor_spip4_Pie'), 'epsc');
saveas(gcf,fullfile(filepath1, 'Cor_spip4_Pie'), 'jpeg');
%% Extracting rows for layers 5 and 6 %%
% Pre-allocate the arrays for L5 and L6
dataTable_Areas_L5 = []; % Will store the rows copied directly from GFP
dataTable_Areas_L6 = []; % Will store the summed rows from GFP

% Initialize empty tables
dataTable_Layer5 = [];  % For rows 1, 4, 7, etc.
dataTable_Layer6a = []; % For rows 2, 5, 8, etc.
dataTable_Layer6b = []; % For rows 3, 6, 9, etc.

% Loop through each row of dataTable5
for i = 1:size(dataTable5, 1)
    if mod(i, 3) == 1
        % For rows where the index is 1 modulo 3 (1st, 4th, 7th, etc.), copy to Layer5
        dataTable_Layer5 = [dataTable_Layer5; dataTable5(i, :)];
    elseif mod(i, 3) == 2
        % For rows where the index is 2 modulo 3 (2nd, 5th, 8th, etc.), copy to Layer6a
        dataTable_Layer6a = [dataTable_Layer6a; dataTable5(i, :)];
    else
        % For rows where the index is 0 modulo 3 (3rd, 6th, 9th, etc.), copy to Layer6b
        dataTable_Layer6b = [dataTable_Layer6b; dataTable5(i, :)];
    end
end
%% Extracting layer areas and their names %

dataTable_Areas_Layer5 = table2array(dataTable_Layer5(:,1));
dataTable_Area_names_L5 = cellfun(@(x) x(1:(length(x)-2)), table2array(dataTable_Layer5(:,3)), 'UniformOutput', false);

dataTable_Areas_Layer6a = table2array(dataTable_Layer6a(:,1));
dataTable_Area_names_L6a = cellfun(@(x) x(1:(length(x)-2)), table2array(dataTable_Layer6a(:,3)), 'UniformOutput', false);

dataTable_Areas_Layer6b = table2array(dataTable_Layer6b(:,1));
dataTable_Area_names_L6b = cellfun(@(x) x(1:(length(x)-2)), table2array(dataTable_Layer6b(:,3)), 'UniformOutput', false);

%% Getting GFP/RFP data for selected layers each area
% Layer 5 %
% ipsi %
for i = 1:size(dataTable_Areas_Layer5,1)
    for j = 1:size(dataTable_GFP,1)
        if dataTable_Areas_Layer5(i,1) == dataTable_GFP(j,1)
            dataTable_L5_ipsi_GFP(i,1) = dataTable_GFP(j,1);
            dataTable_L5_ipsi_GFP(i,2) = dataTable_GFP(j,2);
            dataTable_L5_ipsi_GFP(i,3) = dataTable_GFP(j,3);
            dataTable_L5_ipsi_GFP(i,4) = dataTable_GFP(j,4);
        end
        if dataTable_Areas_Layer5(i,1) == dataTable_RFP(j,1)
            dataTable_L5_ipsi_RFP(i,1) = dataTable_RFP(j,1);
            dataTable_L5_ipsi_RFP(i,2) = dataTable_RFP(j,2);
            dataTable_L5_ipsi_RFP(i,3) = dataTable_RFP(j,3);
            dataTable_L5_ipsi_RFP(i,4) = dataTable_RFP(j,4);
        end
    end
end

% contra %
for i = 1:size(dataTable_Areas_Layer5,1)
    for j = 1:size(dataTable_GFP,1)
        if (dataTable_Areas_Layer5(i,1)-10000) == dataTable_GFP(j,1)
            dataTable_L5_contra_GFP(i,1) = dataTable_GFP(j,1);
            dataTable_L5_contra_GFP(i,2) = dataTable_GFP(j,2);
            dataTable_L5_contra_GFP(i,3) = dataTable_GFP(j,3);
            dataTable_L5_contra_GFP(i,4) = dataTable_GFP(j,4);
        end
        if (dataTable_Areas_Layer5(i,1)-10000) == dataTable_RFP(j,1)
            dataTable_L5_contra_RFP(i,1) = dataTable_RFP(j,1);
            dataTable_L5_contra_RFP(i,2) = dataTable_RFP(j,2);
            dataTable_L5_contra_RFP(i,3) = dataTable_RFP(j,3);
            dataTable_L5_contra_RFP(i,4) = dataTable_RFP(j,4);
        end
    end
end

%% Layer 6a %%
% ipsi %
for i = 1:size(dataTable_Areas_Layer6a,1)
    for j = 1:size(dataTable_GFP,1)
        if dataTable_Areas_Layer6a(i,1) == dataTable_GFP(j,1)
            dataTable_L6a_ipsi_GFP(i,1) = dataTable_GFP(j,1);
            dataTable_L6a_ipsi_GFP(i,2) = dataTable_GFP(j,2);
            dataTable_L6a_ipsi_GFP(i,3) = dataTable_GFP(j,3);
            dataTable_L6a_ipsi_GFP(i,4) = dataTable_GFP(j,4);
        end
        if dataTable_Areas_Layer6a(i,1) == dataTable_RFP(j,1)
            dataTable_L6a_ipsi_RFP(i,1) = dataTable_RFP(j,1);
            dataTable_L6a_ipsi_RFP(i,2) = dataTable_RFP(j,2);
            dataTable_L6a_ipsi_RFP(i,3) = dataTable_RFP(j,3);
            dataTable_L6a_ipsi_RFP(i,4) = dataTable_RFP(j,4);
        end
    end
end

% contra %
for i = 1:size(dataTable_Areas_Layer6a,1)
    for j = 1:size(dataTable_GFP,1)
        if (dataTable_Areas_Layer6a(i,1)-10000) == dataTable_GFP(j,1)
            dataTable_L6a_contra_GFP(i,1) = dataTable_GFP(j,1);
            dataTable_L6a_contra_GFP(i,2) = dataTable_GFP(j,2);
            dataTable_L6a_contra_GFP(i,3) = dataTable_GFP(j,3);
            dataTable_L6a_contra_GFP(i,4) = dataTable_GFP(j,4);
        end
        if (dataTable_Areas_Layer6a(i,1)-10000) == dataTable_RFP(j,1)
            dataTable_L6a_contra_RFP(i,1) = dataTable_RFP(j,1);
            dataTable_L6a_contra_RFP(i,2) = dataTable_RFP(j,2);
            dataTable_L6a_contra_RFP(i,3) = dataTable_RFP(j,3);
            dataTable_L6a_contra_RFP(i,4) = dataTable_RFP(j,4);
        end
    end
end

%% Layer 6b %%
% ipsi %
for i = 1:size(dataTable_Areas_Layer6b,1)
    for j = 1:size(dataTable_GFP,1)
        if dataTable_Areas_Layer6b(i,1) == dataTable_GFP(j,1)
            dataTable_L6b_ipsi_GFP(i,1) = dataTable_GFP(j,1);
            dataTable_L6b_ipsi_GFP(i,2) = dataTable_GFP(j,2);
            dataTable_L6b_ipsi_GFP(i,3) = dataTable_GFP(j,3);
            dataTable_L6b_ipsi_GFP(i,4) = dataTable_GFP(j,4);
        end
        if dataTable_Areas_Layer6b(i,1) == dataTable_RFP(j,1)
            dataTable_L6b_ipsi_RFP(i,1) = dataTable_RFP(j,1);
            dataTable_L6b_ipsi_RFP(i,2) = dataTable_RFP(j,2);
            dataTable_L6b_ipsi_RFP(i,3) = dataTable_RFP(j,3);
            dataTable_L6b_ipsi_RFP(i,4) = dataTable_RFP(j,4);
        end
    end
end

% contra %
for i = 1:size(dataTable_Areas_Layer6b,1)
    for j = 1:size(dataTable_GFP,1)
        if (dataTable_Areas_Layer6b(i,1)-10000) == dataTable_GFP(j,1)
            dataTable_L6b_contra_GFP(i,1) = dataTable_GFP(j,1);
            dataTable_L6b_contra_GFP(i,2) = dataTable_GFP(j,2);
            dataTable_L6b_contra_GFP(i,3) = dataTable_GFP(j,3);
            dataTable_L6b_contra_GFP(i,4) = dataTable_GFP(j,4);
        end
        if (dataTable_Areas_Layer6b(i,1)-10000) == dataTable_RFP(j,1)
            dataTable_L6b_contra_RFP(i,1) = dataTable_RFP(j,1);
            dataTable_L6b_contra_RFP(i,2) = dataTable_RFP(j,2);
            dataTable_L6b_contra_RFP(i,3) = dataTable_RFP(j,3);
            dataTable_L6b_contra_RFP(i,4) = dataTable_RFP(j,4);
        end
    end
end

%%
% Data arrays
L5_GFP = [dataTable_L5_ipsi_GFP(:, [2, 3, 4]);  dataTable_L5_contra_GFP(:, [2, 3, 4])];
L5_RFP = [dataTable_L5_ipsi_RFP(:, [2, 3, 4]); dataTable_L5_contra_RFP(:, [2, 3, 4])];


L6a_GFP = [dataTable_L6a_ipsi_GFP(:, [2, 3, 4]);  dataTable_L6a_contra_GFP(:, [2, 3, 4])];
L6a_RFP = [dataTable_L6a_ipsi_RFP(:, [2, 3, 4]); dataTable_L6a_contra_RFP(:, [2, 3, 4])];



L6b_GFP = [dataTable_L6b_ipsi_GFP(:, [2, 3, 4]);  dataTable_L6b_contra_GFP(:, [2, 3, 4])];
L6b_RFP = [dataTable_L6b_ipsi_RFP(:, [2, 3, 4]); dataTable_L6b_contra_RFP(:, [2, 3, 4])];


% Summing the corresponding rows of L6a_GFP and L6b_GFP
L6_GFP = L6a_GFP + L6b_GFP;

% Summing the corresponding rows of L6a_RFP and L6b_RFP
L6_RFP = L6a_RFP + L6b_RFP;



%% Pie chart  for layer 5 %%
% Calculate the total of L5_GFP and L5_RFP
total_L5_GFP = sum(L5_GFP(:));  % Sum of all values in L5_GFP
total_L5_RFP = sum(L5_RFP(:));  % Sum of all values in L5_RFP

% Prepare the data for the pie chart
data = [total_L5_GFP, total_L5_RFP];

total_sum_L5 = total_L5_GFP + total_L5_RFP;

% Labels with percentages
labels = {sprintf('Indirect (GFP)\n%.1f%%', (total_L5_GFP / total_sum_L5) * 100), ...
    sprintf('Direct (RFP)\n%.1f%%', (total_L5_RFP / total_sum_L5) * 100)};
% Create the pie chart
figure;
h = pie(data, labels);
colormap([1, 0, 1; 1, 1, 0]); % Magenta for GFP, Yellow for RFP

% Set title for the pie chart
title_handle = title(' L5 Corticospinal P04 Neurons', 'FontName', 'Arial', 'FontWeight', 'bold', 'FontSize', 14);

% Adjust the title position (move it up)
title_handle.Position(2) = title_handle.Position(2) + 0.05; % Adjust the y-position slightly higher

% Set font to Arial and bold for all pie chart labels
for k = 1:length(h)
    if isprop(h(k), 'Text')  % Only set the font for text properties
        h(k).FontName = 'Arial';
        h(k).FontWeight = 'bold';
        h(k).FontSize = 12;
    end
end

% Display the pie chart
set(gca, 'FontSize', 12, 'FontWeight', 'bold');
% save figure %
saveas(gcf,fullfile(filepath1, 'Cor_spip4_pie_L5'), 'epsc');
saveas(gcf,fullfile(filepath1, 'Cor_spip4_pie_L5'), 'jpeg');
%% Pie chart  for layer 6 %%
% Calculate the total of L6_GFP and L6_RFP
total_L6_GFP = sum(L6_GFP(:));  % Sum of all values in L5_GFP
total_L6_RFP = sum(L6_RFP(:));  % Sum of all values in L5_RFP

% Prepare the data for the pie chart
data = [total_L6_GFP, total_L6_RFP];

total_sum_L6 = total_L6_GFP + total_L6_RFP;

% Labels with percentages
labels = {sprintf('Indirect (GFP)\n%.1f%%', (total_L6_GFP / total_sum_L6) * 100), ...
    sprintf('Direct (RFP)\n%.1f%%', (total_L6_RFP / total_sum_L6) * 100)};
% Create the pie chart
figure;
h = pie(data, labels);

% Change the colors of the pie chart slices (patches)
h(1).FaceColor = 'magenta';  % Set color for L5 GFP (first slice)
h(3).FaceColor = 'yellow';   % Set color for L5 RFP (second slice)


% Set title for the pie chart
title_handle = title(' L6 Corticospinal P04 Neurons', 'FontName', 'Arial', 'FontWeight', 'bold', 'FontSize', 14);

% Adjust the title position (move it up)
title_handle.Position(2) = title_handle.Position(2) + 0.05; % Adjust the y-position slightly higher

% Set font to Arial and bold for all pie chart labels
for k = 1:length(h)
    if isprop(h(k), 'Text')  % Only set the font for text properties
        h(k).FontName = 'Arial';
        h(k).FontWeight = 'bold';
        h(k).FontSize = 12;
    end
end

% Display the pie chart
set(gca, 'FontSize', 12, 'FontWeight', 'bold');
% save figure %
saveas(gcf,fullfile(filepath1, 'Cor_spip4_pie_L6'), 'epsc');
saveas(gcf,fullfile(filepath1, 'Cor_spip4_pie_L6'), 'jpeg');
%% Getting GFP/RFP data for comparision of selected areas %%
% ipsi
for i = 1:size((dataTable_Areas_comp),1)
    for j = 1:size((dataTable_GFP),1)
        if dataTable_Areas_comp(i,1) == dataTable_GFP(j,1)
            dataTable_comp_ipsi_GFP(i,1) = dataTable_GFP(j,1);
            dataTable_comp_ipsi_GFP(i,2) = dataTable_GFP(j,2);
            dataTable_comp_ipsi_GFP(i,3) = dataTable_GFP(j,3);
            dataTable_comp_ipsi_GFP(i,4) = dataTable_GFP(j,4);
        end
        if dataTable_Areas_comp(i,1) == dataTable_RFP(j,1)
            dataTable_comp_ipsi_RFP(i,1) = dataTable_RFP(j,1);
            dataTable_comp_ipsi_RFP(i,2) = dataTable_RFP(j,2);
            dataTable_comp_ipsi_RFP(i,3) = dataTable_RFP(j,3);
            dataTable_comp_ipsi_RFP(i,4) = dataTable_RFP(j,4);
        end
    end
end

% contra %
for i = 1:size((dataTable_Areas_comp),1)
    for j = 1:size((dataTable_GFP),1)
        if (dataTable_Areas_comp(i,1) - 10000) == dataTable_GFP(j,1)
            dataTable_comp_contra_GFP(i,1) = dataTable_GFP(j,1);
            dataTable_comp_contra_GFP(i,2) = dataTable_GFP(j,2);
            dataTable_comp_contra_GFP(i,3) = dataTable_GFP(j,3);
            dataTable_comp_contra_GFP(i,4) = dataTable_GFP(j,4);
        end
        if (dataTable_Areas_comp(i,1) - 10000) == dataTable_RFP(j,1)
            dataTable_comp_contra_RFP(i,1) = dataTable_RFP(j,1);
            dataTable_comp_contra_RFP(i,2) = dataTable_RFP(j,2);
            dataTable_comp_contra_RFP(i,3) = dataTable_RFP(j,3);
            dataTable_comp_contra_RFP(i,4) = dataTable_RFP(j,4);
        end
    end
end
writematrix(dataTable_comp_ipsi_GFP, fullfile(filepath1, 'dataTable_comp_ipsi_GFP.xlsx'));
writematrix(dataTable_comp_ipsi_RFP, fullfile(filepath1, 'dataTable_comp_ipsi_RFP.xlsx'));
writematrix(dataTable_comp_contra_GFP, fullfile(filepath1, 'dataTable_comp_contra_GFP.xlsx'));
writematrix(dataTable_comp_contra_RFP, fullfile(filepath1, 'dataTable_comp_contra_RFP.xlsx'));
%% Zoomed bar graph for ipsilateral: Highlight RFP and truncate GFP %%
figure('Units', 'inches', 'Position', [2, 2, 50/25.4, 200/25.4]);
hold on;

% Combine values
ipsi_combined = [ipsi_Lev3_RFP_mean, min(ipsi_Lev3_GFP_mean, 1.05 * max(ipsi_Lev3_RFP_mean))]; % Truncate GFP
y_ipsi = size(ipsi_combined, 1):-1:1;

% Grouped bar plot
bh = barh(y_ipsi, ipsi_combined, 0.8, 'grouped', 'EdgeColor', 'none');
bh(1).FaceColor = [1 1 0]; % RFP - yellow
bh(2).FaceColor = [1 0 1]; % GFP - magenta

% Overlay truncated markers for GFP
for i = 1:length(y_ipsi)
    if ipsi_Lev3_GFP_mean(i) > 1.05 * max(ipsi_Lev3_RFP_mean)
        % Add cap mark at truncation limit
        line([1.05 * max(ipsi_Lev3_RFP_mean), 1.05 * max(ipsi_Lev3_RFP_mean)], ...
            [y_ipsi(i)-0.3, y_ipsi(i)+0.3], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 0.5);
    end
end

% Overlay individual data points
for i = 1:length(y_ipsi)
    x_rfp = dataTable_Lev3_ipsi_RFP(i, 2:4);
    x_gfp = min(dataTable_Lev3_ipsi_GFP(i, 2:4), 1.05 * max(ipsi_Lev3_RFP_mean)); % Truncate GFP
    y_rfp = repmat(y_ipsi(i) + 0.15, size(x_rfp));
    y_gfp = repmat(y_ipsi(i) - 0.15, size(x_gfp));

    plot(x_rfp, y_rfp, '-k', 'LineWidth', 0.5);
    scatter(x_rfp, y_rfp, 4, 'k', 'filled');

    plot(x_gfp, y_gfp, '-k', 'LineWidth', 0.5);
    scatter(x_gfp, y_gfp, 4, 'k', 'filled');
end

% Value labels for RFP
for i = 1:length(y_ipsi)
    text(ipsi_Lev3_RFP_mean(i) + 0.1, y_ipsi(i) + 0.25, ...
        num2str(ipsi_Lev3_RFP_mean(i), '%.1f'), 'FontSize', 4, ...
        'Color', 'k', 'HorizontalAlignment', 'left');
end

% Label truncated GFPs
for i = 1:length(y_ipsi)
    if ipsi_Lev3_GFP_mean(i) > 1.05 * max(ipsi_Lev3_RFP_mean)
        text(1.05 * max(ipsi_Lev3_RFP_mean) + 0.1, y_ipsi(i) - 0.35, ...
            sprintf('>%.1f', 1.05 * max(ipsi_Lev3_RFP_mean)), ...
            'FontSize', 8, 'Color', 'k', 'HorizontalAlignment', 'left');
    else
        text(ipsi_Lev3_GFP_mean(i) + 0.1, y_ipsi(i) - 0.35, ...
            num2str(ipsi_Lev3_GFP_mean(i), '%.1f'), ...
            'FontSize', 8, 'Color', 'k', 'HorizontalAlignment', 'left');
    end
end

% Axis labels and title
ylabel('Isocortical Areas', 'FontSize', 5);
xlabel('Number of Labeled Neurons', 'FontName', 'Arial', 'FontWeight', 'bold', 'FontSize', 5);
title('Corticotectal Ipsilateral Projections (RFP Zoomed)', ...
    'FontName', 'Arial', 'FontWeight', 'bold', 'FontSize', 5);

% Y-ticks and labels
yticks(1:size(ipsi_Lev3_GFP_mean, 1));
yticklabels(flip(dataTable_Area_Lev3_names));

% Axes style
ax = gca;
ax.YAxis.FontName = 'Arial';
ax.YAxis.FontWeight = 'bold';
ax.YAxis.FontSize = 5;
ax.XAxis.FontName = 'Arial';
ax.XAxis.FontWeight = 'bold';
ax.XAxis.FontSize = 5;

% Set x-axis limit to zoom on RFP
xlim([0, 1.2 * max(ipsi_Lev3_RFP_mean)]);

hold off;

% Save figure
saveas(gcf, fullfile(filepath1, 'Cor_tec_ipsi_only_zoomRFP'), 'epsc');
%% Zoomed bar graph for contralateral: Highlight RFP and truncate GFP %%
figure('Units', 'inches', 'Position', [2, 2, 50/25.4, 200/25.4]);
hold on;

% Negate for contralateral, truncate GFP at RFP zoom level
rfp_max = max(contra_Lev3_RFP_mean);
gfp_trunc_limit = 1.05 * rfp_max;

contra_combined = [-contra_Lev3_RFP_mean, -min(contra_Lev3_GFP_mean, gfp_trunc_limit)];
y_contra = size(contra_combined, 1):-1:1;

% Grouped bar plot
bh = barh(y_contra, contra_combined, 0.8, 'grouped', 'EdgeColor', 'none');
bh(1).FaceColor = [1 1 0]; % RFP - yellow
bh(2).FaceColor = [1 0 1]; % GFP - magenta

% Overlay truncation marks for GFP
for i = 1:length(y_contra)
    if contra_Lev3_GFP_mean(i) > gfp_trunc_limit
        % Truncation line at the left cutoff
        x_trunc = -gfp_trunc_limit;
        line([x_trunc, x_trunc], [y_contra(i)-0.3, y_contra(i)+0.3], ...
            'Color', 'k', 'LineStyle', '--', 'LineWidth', 0.5);
    end
end

% Overlay individual data points
for i = 1:length(y_contra)
    x_rfp = -dataTable_Lev3_contra_RFP(i, 2:4);
    x_gfp = -min(dataTable_Lev3_contra_GFP(i, 2:4), gfp_trunc_limit); % Truncate GFP
    y_rfp = repmat(y_contra(i) + 0.15, size(x_rfp));
    y_gfp = repmat(y_contra(i) - 0.15, size(x_gfp));

    plot(x_rfp, y_rfp, '-k', 'LineWidth', 0.5);
    scatter(x_rfp, y_rfp, 4, 'k', 'filled');

    plot(x_gfp, y_gfp, '-k', 'LineWidth', 0.5);
    scatter(x_gfp, y_gfp, 4, 'k', 'filled');
end

% Add text labels
for i = 1:length(y_contra)
    % RFP label
    text(-contra_Lev3_RFP_mean(i) - 0.1, y_contra(i) + 0.25, ...
        num2str(contra_Lev3_RFP_mean(i), '%.1f'), ...
        'FontSize', 4, 'Color', 'k', 'HorizontalAlignment', 'right');

    % GFP label (truncated or exact)
    if contra_Lev3_GFP_mean(i) > gfp_trunc_limit
        text(-gfp_trunc_limit - 0.1, y_contra(i) - 0.35, ...
            sprintf('>%.1f', gfp_trunc_limit), ...
            'FontSize', 4, 'Color', 'k', 'HorizontalAlignment', 'right');
    else
        text(-contra_Lev3_GFP_mean(i) - 0.1, y_contra(i) - 0.35, ...
            num2str(contra_Lev3_GFP_mean(i), '%.1f'), ...
            'FontSize', 4, 'Color', 'k', 'HorizontalAlignment', 'right');
    end
end

% Axis labels and title
ylabel('Isocortical Areas', 'FontSize', 5);
xlabel('Number of Labeled Neurons', 'FontName', 'Arial', 'FontWeight', 'bold', 'FontSize', 5);
title('Corticotectal Contralateral Projections (RFP Zoomed)', ...
    'FontName', 'Arial', 'FontWeight', 'bold', 'FontSize', 5);

% Y-ticks and labels
yticks(1:size(contra_Lev3_GFP_mean, 1));
yticklabels(flip(dataTable_Area_Lev3_names));

% Axes formatting
ax = gca;
ax.YAxis.FontName = 'Arial';
ax.YAxis.FontWeight = 'bold';
ax.YAxis.FontSize = 5;
ax.XAxis.FontName = 'Arial';
ax.XAxis.FontWeight = 'bold';
ax.XAxis.FontSize = 5;

% Set x-axis limit to zoom in on RFP
xlim([-1.2 * rfp_max, 0]);

hold off;

% Save figure
saveas(gcf, fullfile(filepath1, 'Cor_tec_contra_only_zoomRFP'), 'epsc');

%% Zero plots %%
% Shared Axis Bar Plot (Log Scale, Skipping NaNs)

% Find Level 1 areas with RFP mean = 0
zero_idx_Lev1_ipsi = find(ipsi_Lev1_RFP_mean == 0);
zero_idx_Lev1_contra = find(contra_Lev1_RFP_mean == 0);

% Extract raw GFP values for these areas
gfp_zero_ipsi_raw = floor(ipsi_Lev1_GFP_mean(zero_idx_Lev1_ipsi));
gfp_zero_ipsi_raw(gfp_zero_ipsi_raw == 1) = 0;
gfp_zero_contra_raw = floor(contra_Lev1_GFP_mean(zero_idx_Lev1_contra));
gfp_zero_contra_raw(gfp_zero_contra_raw == 1) = 0;


% Replace zeros with NaN to skip in log10
gfp_zero_ipsi_raw(gfp_zero_ipsi_raw == 0) = NaN;
gfp_zero_contra_raw(gfp_zero_contra_raw == 0) = NaN;

% Compute log10
gfp_zero_ipsi = log10(gfp_zero_ipsi_raw);
gfp_zero_contra = log10(gfp_zero_contra_raw);

% Get area names
area_Lev1_ipsi_zero_names = dataTable_Area_Lev1_names(zero_idx_Lev1_ipsi);
area_Lev1_contra_zero_names = dataTable_Area_Lev1_names(zero_idx_Lev1_contra);

% Remove entries with NaNs (i.e., originally 0 values)
valid_idx_ipsi = ~isnan(gfp_zero_ipsi);
valid_idx_contra = ~isnan(gfp_zero_contra);

gfp_zero_ipsi = gfp_zero_ipsi(valid_idx_ipsi);
area_Lev1_ipsi_zero_names = area_Lev1_ipsi_zero_names(valid_idx_ipsi);

gfp_zero_contra = gfp_zero_contra(valid_idx_contra);
area_Lev1_contra_zero_names = area_Lev1_contra_zero_names(valid_idx_contra);

% Combine and reverse area names for plot
all_names = [area_Lev1_ipsi_zero_names; area_Lev1_contra_zero_names];
all_names = unique(all_names, 'stable');
all_names = flipud(all_names);
y = 1:length(all_names);

% Match indices
[~, idx_ipsi] = ismember(area_Lev1_ipsi_zero_names, all_names);
[~, idx_contra] = ismember(area_Lev1_contra_zero_names, all_names);

% Fill log values into full list (initialize with zeros)
gfp_ipsi = zeros(size(all_names));
gfp_contra = zeros(size(all_names));

gfp_ipsi(idx_ipsi) = gfp_zero_ipsi;
gfp_contra(idx_contra) = -gfp_zero_contra; % mirror for contra

% Plot
figure('Name', 'Level 1 Zero RFP Areas - Log10 Horizontal Bar Plot', ...
    'Color', 'w', 'Position', [100 100 500 500]);
hold on;

barh(y, gfp_contra, 0.4, 'FaceColor', [1 0 1], 'EdgeColor', 'k');  % Contra: cyan
barh(y, gfp_ipsi, 0.4, 'FaceColor', [1 0 1], 'EdgeColor', 'k');    % Ipsi: magenta

xline(0, '--k');
yticks(y);
yticklabels(all_names);
xlabel('Signal (log_{10} scale, Contra Negative | Ipsi Positive)', 'FontSize', 12);
title('Level 1 Areas with RFP = 0 | Mirrored Log10 Plot', 'FontSize', 14);
grid on;
box on;

% Axis formatting
ax = gca;
ax.YAxis.FontName = 'Arial';
ax.YAxis.FontWeight = 'bold';
ax.YAxis.FontSize = 8;

ax.XAxis.FontName = 'Arial';
ax.XAxis.FontWeight = 'bold';
ax.XAxis.FontSize = 8;

% Set X-limits slightly padded beyond max abs log10 value
xlim([-max(abs([gfp_ipsi; gfp_contra]))*1.2, max(abs([gfp_ipsi; gfp_contra]))*1.2]);

% Custom xticks and labels with antilog format
xticks_vals = -3:1:3;
xticks(xticks_vals);
xticklabels(arrayfun(@(x) sprintf('10^{%d}', abs(x)), xticks_vals, 'UniformOutput', false));

% Save figure (ensure filepath1 is defined)
saveas(gcf, fullfile(filepath1, 'Cor_tec_Zero_Log10'), 'epsc');

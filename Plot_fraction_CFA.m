clear all; close all; clc
% Specify the file path
disp('Select Br-1/Br-2 GFP Excel File');
[filename1, filepath1] = uigetfile('*.xlsx', 'Select GFP Excel File');
% disp('Select Br-1 RFP Excel File');
% [filename2, filepath2] = uigetfile('*.csv', 'Select RFP Excel File');
disp('Select Br-1/Br-2 RFP Excel File');
[filename2, filepath2] = uigetfile('*.xlsx', 'Select GFP Excel File');
% disp('Select Br-2 RFP Excel File');
% [filename5, filepath5] = uigetfile('*.csv', 'Select GFP Excel File');
disp('Select Area ID Excel File');
[filename3, filepath3] = uigetfile('*.xlsx', 'Select Area ID Excel File');

fullFileName1 = fullfile(filepath1, filename1);
fullFileName2 = fullfile(filepath2, filename2);
% fullFileName4 = fullfile(filepath4, filename4);
% fullFileName5 = fullfile(filepath5, filename5);
fullFileName3 = fullfile(filepath3, filename3);

disp(['Selected file: ', fullFileName1]);
disp(['Selected file: ', fullFileName2]);
% disp(['Selected file: ', fullFileName3]);
% disp(['Selected file: ', fullFileName4]);
% disp(['Selected file: ', fullFileName5]);
%%
% Read the data into a table
dataTable1 = readtable(filename1);
dataTable2 = readtable(filename2);
dataTable3 = readtable(filename3);
% dataTable4 = readtable(filename4);
% dataTable5 = readtable(filename5);

%% convert tables to arrays
% using fraction occupation %
dataTable_GFP(:,[1,2]) = table2array(dataTable1(:,[1,4]));
dataTable_GFP(:,3) = table2array(dataTable1(:,5));
dataTable_RFP(:,[1,2]) = table2array(dataTable2(:,[1,4]));
dataTable_RFP(:,3) = table2array(dataTable2(:,5));
% extracting selected areas and their names %
dataTable_Areas = table2array(dataTable3(:,1));
dataTable_Area_names = cellfun(@(x) x(1:(length(x)-2)), table2array(dataTable3(:,3)), 'UniformOutput', false);
% saving %
writecell(dataTable_Area_names, fullfile(filepath1, 'CFA_Area_Names.xlsx'));

%% Getting GFP/RFP data for selected areas
% ipsi
for i = 1:size((dataTable_Areas),1)
for j = 1:size((dataTable_GFP),1)
if dataTable_Areas(i,1) == dataTable_GFP(j,1)
dataTable_sel_ipsi_GFP(i,1) = dataTable_GFP(j,1);
dataTable_sel_ipsi_GFP(i,2) = dataTable_GFP(j,2);
dataTable_sel_ipsi_GFP(i,3) = dataTable_GFP(j,3);
end
if dataTable_Areas(i,1) == dataTable_RFP(j,1)
dataTable_sel_ipsi_RFP(i,1) = dataTable_RFP(j,1);
dataTable_sel_ipsi_RFP(i,2) = dataTable_RFP(j,2);
dataTable_sel_ipsi_RFP(i,3) = dataTable_RFP(j,3);
end
end
end

%% contra
for i = 1:size((dataTable_Areas),1)
for j = 1:size((dataTable_GFP),1)
if (dataTable_Areas(i,1) - 10000) == dataTable_GFP(j,1)
dataTable_sel_contra_GFP(i,1) = dataTable_GFP(j,1);
dataTable_sel_contra_GFP(i,2) = dataTable_GFP(j,2);
dataTable_sel_contra_GFP(i,3) = dataTable_GFP(j,3);
end
if (dataTable_Areas(i,1) - 10000) == dataTable_RFP(j,1)
dataTable_sel_contra_RFP(i,1) = dataTable_RFP(j,1);
dataTable_sel_contra_RFP(i,2) = dataTable_RFP(j,2);
dataTable_sel_contra_RFP(i,3) = dataTable_RFP(j,3);
end
end
end
% saving %
writematrix(dataTable_sel_contra_GFP , fullfile(filepath1, 'CFA_dataTable_sel_contra_GFP.xlsx'));
writematrix(dataTable_sel_contra_RFP, fullfile(filepath1, 'CFA_dataTable_sel_contra_RFP.xlsx'));
writematrix(dataTable_sel_ipsi_GFP , fullfile(filepath1, 'CFA_dataTable_sel_ipsi_GFP.xlsx'));
writematrix(dataTable_sel_ipsi_RFP , fullfile(filepath1, 'CFA_dataTable_sel_ipsi_RFP.xlsx'));

%% Using fraction occupation %%
%
dataTable_sel_ipsi_Norm_GFP(:,1) = dataTable_sel_ipsi_GFP(:,1);
dataTable_sel_contra_Norm_GFP(:,1) = dataTable_sel_contra_GFP(:,1);
dataTable_sel_ipsi_Norm_GFP(:,2) = dataTable_sel_ipsi_GFP(:,2);
dataTable_sel_contra_Norm_GFP(:,2) = dataTable_sel_contra_GFP(:,2);
dataTable_sel_ipsi_Norm_GFP(:,3) = dataTable_sel_ipsi_GFP(:,3);
dataTable_sel_contra_Norm_GFP(:,3) = dataTable_sel_contra_GFP(:,3);

dataTable_sel_ipsi_Norm_RFP(:,1) = dataTable_sel_ipsi_RFP(:,1);
dataTable_sel_contra_Norm_RFP(:,1) = dataTable_sel_contra_RFP(:,1);
dataTable_sel_ipsi_Norm_RFP(:,2) = dataTable_sel_ipsi_RFP(:,2);
dataTable_sel_contra_Norm_RFP(:,2) = dataTable_sel_contra_RFP(:,2);
dataTable_sel_ipsi_Norm_RFP(:,3) = dataTable_sel_ipsi_RFP(:,3);
dataTable_sel_contra_Norm_RFP(:,3) = dataTable_sel_contra_RFP(:,3);

%% Calculate mean %%
ipsi_Norm_GFP_mean = mean(dataTable_sel_ipsi_Norm_GFP(:,[2,3]), 2);
contra_Norm_GFP_mean = mean(dataTable_sel_contra_Norm_GFP(:,[2,3]), 2);


ipsi_Norm_RFP_mean = mean(dataTable_sel_ipsi_Norm_RFP(:,[2,3]), 2);
contra_Norm_RFP_mean = mean(dataTable_sel_contra_Norm_RFP(:,[2,3]), 2);




%% Min-Max normalization global mean
%ipsi
combined_ipsi = [ipsi_Norm_GFP_mean, ipsi_Norm_RFP_mean];
global_min_ipsi = min(min(combined_ipsi));
global_max_ipsi = max(max(combined_ipsi));

ipsi_Norm_GFP_mean_Norm = (ipsi_Norm_GFP_mean - global_min_ipsi)/(global_max_ipsi - global_min_ipsi);
ipsi_Norm_RFP_mean_Norm = (ipsi_Norm_RFP_mean - global_min_ipsi)/(global_max_ipsi - global_min_ipsi);

%contra
combined_contra = [contra_Norm_GFP_mean, contra_Norm_RFP_mean];
global_min_contra = min(min(combined_contra));
global_max_contra = max(max(combined_contra));

contra_Norm_GFP_mean_Norm = (contra_Norm_GFP_mean - global_min_contra)/(global_max_contra - global_min_contra);
contra_Norm_RFP_mean_Norm = (contra_Norm_RFP_mean - global_min_contra)/(global_max_contra - global_min_contra);

% saving %
writematrix(ipsi_Norm_GFP_mean_Norm , fullfile(filepath1, 'CFA_Norm_mean_GFP_ipsi.xlsx'));
writematrix(ipsi_Norm_RFP_mean_Norm , fullfile(filepath1, 'CFA_Norm_mean_RFP_ipsi.xlsx'));
writematrix(ipsi_Norm_GFP_mean_Norm , fullfile(filepath1, 'CFA_Norm_mean_GFP_contra.xlsx'));
writematrix(ipsi_Norm_RFP_mean_Norm , fullfile(filepath1, 'CFA_Norm_mean_RFP_contra.xlsx'));
%%
% Create the horizontal bar graph ipsi %
figure('Units', 'inches', 'Position', [100,100,4,8]);
ye = [0.8500 0.3250 0.0980];
hold on; % Hold on to the current axes
y1 = size((ipsi_Norm_GFP_mean_Norm),1):-1:1;
y2 = (size((ipsi_Norm_GFP_mean_Norm),1):-1:1) + 0.3;

barh(y1, ipsi_Norm_GFP_mean_Norm, 0.3, 'FaceColor', 'y', 'DisplayName', 'Set 1'); % First set in green
barh(y2, ipsi_Norm_RFP_mean_Norm, 0.3, 'FaceColor', 'm', 'DisplayName', 'Set 2'); % Second set in red
% Add labels and title
ylabel('Y Axis Label');
xlabel('Normalized percent occupation','FontName', 'Arial', 'FontWeight', 'bold','FontSize', 5);
title('Ipsilateral Projections','FontName', 'Arial', 'FontWeight', 'bold','FontSize', 5);

% % Customize y-ticks
yy = 1:size(ipsi_Norm_GFP_mean_Norm);
yticks(yy);
yticklabels(flip(dataTable_Area_names));
% Customize y-axis tick labels
ax = gca; % Get the current axes
ax.YAxis.FontName = 'Arial'; % Set font to Arial
ax.YAxis.FontWeight = 'bold'; % Set font weight to bold
ax.YAxis.FontSize = 5; % Set font size to 12
ax.XAxis.FontName = 'Arial'; % Set font to Arial
ax.XAxis.FontWeight = 'bold'; % Set font weight to bold
ax.XAxis.FontSize = 5; % Set font size to 12
hold off

% save figure %
saveas(gcf,fullfile(filepath1, 'Ipsi_bar_proj'), 'epsc');

%% Create the horizontal bar graph ipsi-contra %%
figure('Units', 'inches', 'Position', [100,100,4,8]);
hold on; % Hold on to the current axes

% %ipsi %
y1 = [size((ipsi_Norm_GFP_mean_Norm),1):-1:1];
y2 = [size((ipsi_Norm_GFP_mean_Norm),1):-1:1] + 0.35;

barh(y1, ipsi_Norm_GFP_mean_Norm, 0.3, 'FaceColor', 'y', 'DisplayName', 'Set 1'); % First set in green
barh(y2, ipsi_Norm_RFP_mean_Norm, 0.3, 'FaceColor', 'm', 'DisplayName', 'Set 2'); % Second set in red

hold on
y1 = [size((contra_Norm_GFP_mean_Norm),1):-1:1];
y2 = [size((contra_Norm_GFP_mean_Norm),1):-1:1] + 0.35;

% contra %
barh(y1, -contra_Norm_GFP_mean_Norm, 0.3, 'FaceColor', 'y', 'DisplayName', 'Set 1'); % First set in green
barh(y2, -contra_Norm_RFP_mean_Norm, 0.3, 'FaceColor', 'm', 'DisplayName', 'Set 2'); % Second set in red
% Add labels and title
ylabel('Y Axis Label');
xlabel('Normalized percent occupation','FontName', 'Arial', 'FontWeight', 'bold','FontSize', 7);
title('Contralateral/Ipsilateral Projections','FontName', 'Arial', 'FontWeight', 'bold','FontSize', 7);
hold on
% % Customize y-ticks
yy = 1:size(contra_Norm_GFP_mean_Norm);
yticks(yy);
yticklabels(flip(dataTable_Area_names));
ax = gca; % Get the current axes
ax.YAxis.FontName = 'Arial'; % Set font to Arial
ax.YAxis.FontWeight = 'bold'; % Set font weight to bold
ax.YAxis.FontSize = 5; % Set font size to 12
ax.XAxis.FontName = 'Arial'; % Set font to Arial
ax.XAxis.FontWeight = 'bold'; % Set font weight to bold
ax.XAxis.FontSize = 5; % Set font size to 12
xtick_values = xticks; % Get current x-tick values
xticklabels(arrayfun(@num2str, abs(xtick_values), 'UniformOutput', false)); % Make all labels positive
hold off
% save figure %
saveas(gcf,fullfile(filepath1, 'Ipsi_contra_bar_proj'), 'epsc');

%% colormap plots %%
figure('Units', 'inches', 'Position', [100,100,4,6]);

load('Copper_edited.mat');
% ipsi %
imagesc([ipsi_Norm_GFP_mean_Norm,ipsi_Norm_RFP_mean_Norm]);
colormap(Copper_edited)
colorbar;
title('Ipsilateral Projections','FontName', 'Arial', 'FontWeight', 'bold','FontSize', 7);
yy = 1:size(ipsi_Norm_GFP_mean_Norm);
yticks(yy);
yticklabels((dataTable_Area_names));
% Customize y-axis tick labels
ax = gca; % Get the current axes
ax.YAxis.FontName = 'Arial'; % Set font to Arial
ax.YAxis.FontWeight = 'bold'; % Set font weight to bold
ax.YAxis.FontSize = 5; % Set font size to 12
ax.XAxis.FontName = 'Arial'; % Set font to Arial
ax.XAxis.FontWeight = 'bold'; % Set font weight to bold
ax.XAxis.FontSize = 5; % Set font size to 12
hold off
% save figure %
saveas(gcf,fullfile(filepath1, 'Colormap_ipsi_proj'), 'epsc');
% contra %

figure('Units', 'inches', 'Position', [100,100,4,6]);
imagesc([contra_Norm_GFP_mean_Norm,contra_Norm_RFP_mean_Norm]);
colormap(Copper_edited)
colorbar;
title('Contralateral Projections','FontName', 'Arial', 'FontWeight', 'bold','FontSize', 7);
yy = 1:size(ipsi_Norm_GFP_mean_Norm);
yticks(yy);
yticklabels(dataTable_Area_names);
% Customize y-axis tick labels
ax = gca; % Get the current axes
ax.YAxis.FontName = 'Arial'; % Set font to Arial
ax.YAxis.FontWeight = 'bold'; % Set font weight to bold
ax.YAxis.FontSize = 5; % Set font size to 12
ax.XAxis.FontName = 'Arial'; % Set font to Arial
ax.XAxis.FontWeight = 'bold'; % Set font weight to bold
ax.XAxis.FontSize = 5; % Set font size to 12
hold off
% save figure %
saveas(gcf,fullfile(filepath1, 'Colormap_contra_proj'), 'epsc');
%% Regression %%
%
% % Fit linear regression model: y = m*x + b
% p_gen = polyfit(ipsi_Norm_GFP_mean, ipsi_Norm_RFP_mean, 1);  % p(1) = slope, p(2) = intercept
%
% % Generate the regression line for plotting
% Y_fit_gen = polyval(p_gen, ipsi_Norm_GFP_mean);
%
% % Different plots %
% % Generic %
% figure ();
% scatter(ipsi_Norm_GFP_mean, ipsi_Norm_RFP_mean, 50,'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'none','LineWidth', 2);
% hold on;
% % Plot the regression line
% plot(ipsi_Norm_GFP_mean, Y_fit_gen, 'k-', 'LineWidth', 2);  % Regression line
% % Add labels from the list of names to each data point
% for i = 1:length(contra_Norm_GFP_mean)
%     text(ipsi_Norm_GFP_mean(i), ipsi_Norm_RFP_mean(i), dataTable_Area_names{i}, ...
%         'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
% end
%
% % Add labels and title
% xlabel('Direct Projections Fraction Occupation');
% ylabel('Indirect Projections Fraction Occupation');
% title('CFA Projections');
% hold off;

%% Thalamus/Brainstem Regression %%

TG = [dataTable_sel_ipsi_GFP(32:47,2);dataTable_sel_contra_GFP(32:47,2);dataTable_sel_ipsi_GFP(32:47,3);dataTable_sel_contra_GFP(32:47,3)];
TR = [dataTable_sel_ipsi_RFP(32:47,2);dataTable_sel_contra_RFP(32:47,2);dataTable_sel_ipsi_RFP(32:47,3);dataTable_sel_contra_RFP(32:47,3)];
BG = [dataTable_sel_ipsi_GFP(48:end,2);dataTable_sel_contra_GFP(48:end,2);dataTable_sel_ipsi_GFP(48:end,3);dataTable_sel_contra_GFP(48:end,3)];
BR = [dataTable_sel_ipsi_RFP(48:end,2);dataTable_sel_contra_RFP(48:end,2);dataTable_sel_ipsi_RFP(48:end,3);dataTable_sel_contra_RFP(48:end,3)];


% Thalamus %

Xt = TG;
Yt = TR;
Zt = [dataTable_Area_names([32:47]);dataTable_Area_names([32:47]);dataTable_Area_names([32:47]);dataTable_Area_names([32:47])];
figure('Units', 'inches', 'Position', [100,100,3,2]);
% Fit linear regression model: y = m*x + b
pt = polyfit(Xt, Yt, 1);  % p(1) = slope, p(2) = intercept

% Generate the regression line for plotting
Y_fitt = polyval(pt, Xt);

% Plot the regression line
plot(Xt, Y_fitt, 'c-', 'LineWidth', 1);  % Regression line
hold on

% Number of parts to split into (4 parts)
numParts = 4;

% Loop through each 1/4th segment of Xt and Yt
for i = 0:numParts-1
% Define the start and end indices for each 1/4th of the data
startIdx = floor(i * length(Xt) / numParts) + 1;
endIdx = floor((i + 1) * length(Xt) / numParts);
if i == 0 || i == 2
% Plot the current 1/4th segment with cyan color
scatter(Xt(startIdx:endIdx), Yt(startIdx:endIdx), 10, 'MarkerEdgeColor', 'cyan', 'MarkerFaceColor', 'cyan', 'LineWidth', 1);
hold on;
else
scatter(Xt(startIdx:endIdx), Yt(startIdx:endIdx), 10, 'MarkerEdgeColor', 'cyan', 'MarkerFaceColor', 'w', 'LineWidth', 1);
hold on;
end
end

% Add labels from the list of names to each data point
for i = 1:length(Xt)
text(Xt(i), Yt(i), Zt{i}, ...
'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right','FontName', 'Arial','FontWeight', 'bold', 'FontSize', 5);
end
hold on
% Brainstem %

Xb = BG;
Yb = BR
Zb = [dataTable_Area_names([48:end]);dataTable_Area_names([48:end]);dataTable_Area_names([48:end]);dataTable_Area_names([48:end])];

% saving %
writematrix(TG, fullfile(filepath1, 'Thalamus_GFP_CFA.xlsx'));
writematrix(TR, fullfile(filepath1, 'Thalamus_RFP_CFA.xlsx'));
writematrix(BG, fullfile(filepath1, 'Brainstem_GFP_CFA.xlsx'));
writematrix(BR, fullfile(filepath1, 'Brainstem_RFP_CFA.xlsx'));
writecell(Zt, fullfile(filepath1, 'Thalamus_Areas_CFA.xlsx'));
writecell(Zb, fullfile(filepath1, 'Brainstem_Areas_CFA.xlsx'));

% Fit linear regression model: y = m*x + b
pb = polyfit(Xb, Yb, 1);  % p(1) = slope, p(2) = intercept

% Generate the regression line for plotting
Y_fitb = polyval(pb, Xb);

% Plot the regression line
plot(Xb, Y_fitb, 'b-', 'LineWidth', 1);  % Regression line
hold on

% Number of parts to split into (4 parts)
numParts = 4;

% Loop through each 1/4th segment of Xb and Yb
for i = 0:numParts-1
% Define the start and end indices for each 1/4th of the data
startIdx = floor(i * length(Xb) / numParts) + 1;
endIdx = floor((i + 1) * length(Xb) / numParts);
if i == 0 || i == 2
% Plot the current 1/4th segment with blue color
scatter(Xb(startIdx:endIdx), Yb(startIdx:endIdx), 10, 'MarkerEdgeColor', 'blue', 'MarkerFaceColor', 'blue', 'LineWidth', 1);
hold on;
else
scatter(Xb(startIdx:endIdx), Yb(startIdx:endIdx), 10, 'MarkerEdgeColor', 'blue', 'MarkerFaceColor', 'w', 'LineWidth', 1);
hold on;
end
end

% Add labels from the list of names to each data point

for i = 1:length(Xb)
text(Xb(i), Yb(i), Zb{i}, ...
'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', ...
'FontName', 'Arial', 'FontSize', 5, 'FontWeight', 'bold', 'Color', 'k');  % 'k' is black
end

% Add labels and title
set(gca, 'FontName', 'Arial', 'FontSize', 5, 'FontWeight', 'bold');
xlabel('Direct Projections Fraction Occupation [%]', 'FontName', 'Arial','FontWeight', 'bold', 'FontSize', 7);
ylabel('Indirect Projections Fraction Occupation [%]', 'FontName', 'Arial', 'FontWeight', 'bold','FontSize', 7);
title('CFA Projections: \color{cyan}Thalamus \color{black}vs \color{blue}Brainstem','FontName', 'Arial', 'FontWeight', 'bold','FontSize', 7);
hold off;
% save figure %
saveas(gcf,fullfile(filepath1, 'Th_vs_BS_reg'), 'epsc');
%% statistical tests including ipsi and contra projections %%

% TG = [dataTable_sel_ipsi_GFP(32:45,2);dataTable_sel_ipsi_GFP(32:45,3)];
% TR = [dataTable_sel_ipsi_RFP(32:45,2);dataTable_sel_ipsi_RFP(32:45,3)];
% BG = [dataTable_sel_ipsi_GFP(46:end,2);dataTable_sel_ipsi_GFP(46:end,2)];
% BR = [dataTable_sel_ipsi_RFP(46:end,2);dataTable_sel_ipsi_RFP(46:end,2)];

% Normalize all data together
all_data = [TG; TR; BG; BR];

mean_val = mean(all_data);  % Global mean
std_val = std(all_data);    % Global standard deviation

TG_norm = (TG - mean_val) / std_val;
TR_norm = (TR - mean_val) / std_val;
BG_norm = (BG - mean_val) / std_val;
BR_norm = (BR - mean_val) / std_val;
%% raincloud plot %%
raincloud_plot(TG, TR, BG, BR, 'TG', 'TR', 'BG', 'BR');


%% Difference %%

diff_bs = [BR_norm-BG_norm];
diff_th = [TR_norm-TG_norm];
%% Indivisually Thalamus and Brainstem  tests%%
% Brainstem %
% [h, ~] = kstest((diff_bs - mean(diff_bs)) / std(diff_bs)); % Kolmogorov-Smirnov test for normality
% if h == 0
%     disp('Data follows normal distribution (p > 0.05).');
% else
%     disp('Data does not follow normal distribution (p <= 0.05).');
% end


% % Shapiro–Wilk test
[h, pValue, W] = swtest(diff_bs, 0.05);  % 0.05 is the significance level

if h == 0
disp('Data follows normal distribution (p > 0.05)');
else
disp('Data does NOT follow normal distribution (p <= 0.05)');
end

% Perform Wilcoxon Signed-Rank Test
pbs = signrank(diff_bs);
disp(['p-value: ', num2str(pbs)]);

% Calculate the median and quartiles
median_val = median(diff_bs);
q1 = prctile(diff_bs, 25); % First quartile
q3 = prctile(diff_bs, 75); % Third quartile

% Create the plot
figure('Units', 'inches', 'Position', [100,100,3,2]);
hold on;

% Plot the box manually
x_center = 1; % X position for diff_bs
box_width = 0.1;
rectangle('Position', [x_center - box_width / 2, q1, box_width, q3 - q1], ...
'EdgeColor', 'blue', 'FaceColor', [1, 1, 1]); % blue box for IQR
plot([x_center - box_width / 2, x_center + box_width / 2], [median_val, median_val], ...
'b-', 'LineWidth', 2); % Line for median

% Overlay individual data points in blue
scatter(x_center * ones(size(diff_bs)), diff_bs, 10, 'MarkerFaceColor', 'blue', ...
'MarkerEdgeColor', 'blue', 'LineWidth', 1.2);

% Add p-value
y_pos = q3 + 0.2 * range(diff_bs); % Position for p-value
text(x_center - box_width, y_pos, sprintf('p = %.20f', pbs), ...
'HorizontalAlignment', 'center', 'FontSize', 5, 'FontWeight', 'bold', 'Color', 'black');

% Beautify plot
xlim([0.5, 1.5]);
xticks(x_center);
xticklabels({'diff\_bs'});
set(gca, 'FontName', 'Arial', 'FontSize', 5, 'FontWeight', 'bold');
ylabel('Values','FontName', 'Arial', 'FontWeight', 'bold','FontSize', 7);
title('CFA:Brainstem','FontName', 'Arial', 'FontWeight', 'bold','FontSize', 7);
grid on;
hold off;
% save figure %
saveas(gcf,fullfile(filepath1, 'BS_stat'), 'epsc');
%% Thalamus %%

% [h, ~] = kstest((diff_th - mean(diff_th)) / std(diff_th)); % Kolmogorov-Smirnov test for normality
% if h == 0
%     disp('Data follows normal distribution (p > 0.05).');
% else
%     disp('Data does not follow normal distribution (p <= 0.05).');
% end


% % Shapiro–Wilk test
[h, pValue, W] = swtest(diff_th, 0.05);  % 0.05 is the significance level

if h == 0
disp('Data follows normal distribution (p > 0.05)');
else
disp('Data does NOT follow normal distribution (p <= 0.05)');
end

% Perform Wilcoxon Signed-Rank Test
pth = signrank(diff_th);
disp(['p-value: ', num2str(pth)]);

% Calculate the median and quartiles
median_val = median(diff_th);
q1 = prctile(diff_th, 25); % First quartile
q3 = prctile(diff_th, 75); % Third quartile

% Create the plot
figure('Units', 'inches', 'Position', [100,100,3,2]);
hold on;

% Plot the box manually
x_center = 1; % X position for diff_bs
box_width = 0.1;
rectangle('Position', [x_center - box_width / 2, q1, box_width, q3 - q1], ...
'EdgeColor', 'cyan', 'FaceColor', [1, 1, 1]); % Cyan box for IQR
plot([x_center - box_width / 2, x_center + box_width / 2], [median_val, median_val], ...
'c-', 'LineWidth', 2); % Line for median

% Overlay individual data points in cyan
scatter(x_center * ones(size(diff_th)), diff_th, 10, 'MarkerFaceColor', 'cyan', ...
'MarkerEdgeColor', 'c', 'LineWidth', 1);

% Add p-value
y_pos = q3 + 0.2 * range(diff_bs); % Position for p-value
text(x_center - box_width, y_pos, sprintf('p = %.5f', pth), ...
'HorizontalAlignment', 'center', 'FontSize', 7, 'FontWeight', 'bold', 'Color', 'black');

% Beautify plot
xlim([0.5, 1.5]);
xticks(x_center);
xticklabels({'\Delta Th'});
set(gca, 'FontName', 'Arial', 'FontSize', 5, 'FontWeight', 'bold');
ylabel('Values','FontName', 'Arial', 'FontWeight', 'bold','FontSize', 7);
title('CFA:Thalamus','FontName', 'Arial', 'FontWeight', 'bold','FontSize', 7);
grid on;
hold off;
% save figure %
saveas(gcf,fullfile(filepath1, 'Th_stat'), 'epsc');


%% Thalamus/Midbrain Regression %%
MG = [dataTable_sel_ipsi_GFP([50:56],2);dataTable_sel_contra_GFP([50:56],2);dataTable_sel_ipsi_GFP([50:56],3);dataTable_sel_contra_GFP([50:56],3)];
MR = [dataTable_sel_ipsi_RFP([50:56],2);dataTable_sel_contra_RFP([50:56],2);dataTable_sel_ipsi_RFP([50:56],3);dataTable_sel_contra_RFP([50:56],3)];

% Thalamus %
figure('Units', 'inches', 'Position', [100,100,3,2]);
% % Fit linear regression model: y = m*x + b
% pt = polyfit(Xt, Yt, 1);  % p(1) = slope, p(2) = intercept
% 
% % Generate the regression line for plotting
% Y_fitt = polyval(pt, Xt);
% 
% % Plot the regression line
% plot(Xt, Y_fitt, 'c-', 'LineWidth', 1);  % Regression line
% hold on
% 
% scatter(Xt, Yt, 10,'MarkerEdgeColor', 'c', 'MarkerFaceColor', 'c','LineWidth', 1);
% hold on;
% % Add labels from the list of names to each data point
% for i = 1:length(Xt)
%     text(Xt(i), Yt(i), Zt{i}, ...
%         'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right','FontName', 'Arial', 'FontSize', 5, 'FontWeight', 'bold', 'Color', 'k');
% end
% hold on

% Midbrain %

Xm = MG;
Ym = MR;
Zm = [dataTable_Area_names([50:56]); dataTable_Area_names([50:56]); dataTable_Area_names([50:56]); dataTable_Area_names([50:56])];

% saving %
writematrix(MG, fullfile(filepath1, 'MB_GFP_PFC.xlsx'));
writematrix(MR, fullfile(filepath1, 'MB_RFP_PFC.xlsx'));
writecell(Zm, fullfile(filepath1, 'MB_Areas_PFC.xlsx'));



% Fit linear regression model: y = m*x + b
m_fit = polyfit(Xm, Ym, 1);  % m_fit(1) = slope, m_fit(2) = intercept
% Define the color 'pink' as RGB values
pink = [1.0000 0.4118 0.7059];

% Generate the regression line for plotting
Y_fitm = polyval(m_fit, Xm);

% Plot the regression line
plot(Xm, Y_fitm, '-', 'Color', pink, 'LineWidth', 1);  % Regression line
hold on
% Version 6: Xm, Ym, color = pink (RGB values)

% Loop through each 1/4th segment of Xm and Ym
for i = 0:numParts-1
% Define the start and end indices for each 1/4th of the data
startIdx = floor(i * length(Xm) / numParts) + 1;
endIdx = floor((i + 1) * length(Xm) / numParts);
if i == 0 || i == 2
% Plot the current 1/4th segment with pink color
scatter(Xm(startIdx:endIdx), Ym(startIdx:endIdx), 20, 'MarkerEdgeColor', pink, 'MarkerFaceColor', pink, 'LineWidth', 1);
hold on;
else
scatter(Xm(startIdx:endIdx), Ym(startIdx:endIdx), 20, 'MarkerEdgeColor', pink, 'MarkerFaceColor', 'w', 'LineWidth', 1);
hold on;
end
end

% Add labels from the list of names to each data point
for i = 1:length(Xm)
text(Xm(i), Ym(i), Zm{i}, ...
'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', ...
'FontName', 'Arial', 'FontSize', 10, 'FontWeight', 'bold', 'Color', 'k');  % 'k' is black
end

% Add labels and title
set(gca, 'FontName', 'Arial', 'FontSize', 5, 'FontWeight', 'bold');
xlabel('Direct Projections Fraction Occupation [%]', 'FontName', 'Arial', 'FontWeight', 'bold', 'FontSize', 7);
ylabel('Indirect Projections Fraction Occupation [%]', 'FontName', 'Arial', 'FontWeight', 'bold', 'FontSize', 7);
title('PFC Projections: \color{cyan}Thalamus \color{black}vs \color[rgb]{1.0000, 0.4118, 0.7059}Midbrain', 'FontName', 'Arial', 'FontWeight', 'bold', 'FontSize', 7);
hold off;
% Save figure %
saveas(gcf, fullfile(filepath1, 'Th_vs_Midbrain_reg'), 'epsc');



%% Statistical tests including ipsi and contra projections %%

% Normalize all data together
all_data = [TG; TR; MG; MR];

mean_val = mean(all_data);  % Global mean
std_val = std(all_data);    % Global standard deviation

MG_norm = (MG - mean_val) / std_val;
MR_norm = (MR - mean_val) / std_val;
TG_norm = (TG - mean_val) / std_val;
TR_norm = (TR - mean_val) / std_val;
%% raincloud plot %%
raincloud_plot(TG, TR, MG, MR, 'TG', 'TR', 'MG', 'MR');
% Difference %

diff_m = [MR_norm - MG_norm];
diff_th = [TR_norm - TG_norm];
%% Individually Thalamus and Midbrain tests %%
% Midbrain %
[h, ~] = kstest((diff_m - mean(diff_m)) / std(diff_m)); % Kolmogorov-Smirnov test for normality
if h == 0
disp('Data follows normal distribution (p > 0.05).');
else
disp('Data does not follow normal distribution (p <= 0.05).');
end

% % Perform Wilcoxon Signed-Rank Test
% m_p = signrank(diff_m);
% disp(['p-value: ', num2str(m_p)]);


% Perform one-sample t-test
[fo_t, m_p] = ttest(diff_m, 0);  % Testing if the mean of diff_ho is significantly different from 0

disp(['t-statistic: ', num2str(fo_t)]);
disp(['p-value: ', num2str(m_p)]);


% Calculate the median and quartiles
median_val = median(diff_m);
q1 = prctile(diff_m, 25); % First quartile
q3 = prctile(diff_m, 75); % Third quartile

% Create the plot
figure('Units', 'inches', 'Position', [100,100,3,2]);
hold on;

% Plot the box manually
x_center = 1; % X position for diff_m
box_width = 0.1;
rectangle('Position', [x_center - box_width / 2, q1, box_width, q3 - q1], ...
'EdgeColor', [1.0000 0.4118 0.7059], 'FaceColor', [1, 1, 1]); % Pink box for IQR
plot([x_center - box_width / 2, x_center + box_width / 2], [median_val, median_val], ...
'-', 'Color', [1.0000 0.4118 0.7059], 'LineWidth', 1); % Line for median

% Overlay individual data points in pink
scatter(x_center * ones(size(diff_m)), diff_m, 10, 'MarkerFaceColor', [1.0000 0.4118 0.7059], ...
'MarkerEdgeColor', [1.0000 0.0000 0.0000], 'LineWidth', 1);

% Add p-value
y_pos = q3 + 0.2 * range(diff_m); % Position for p-value
text(x_center - box_width, y_pos, sprintf('p = %.5f', m_p), ...
'HorizontalAlignment', 'center', 'FontSize', 7, 'FontWeight', 'bold', 'Color', [1.0000 0.4118 0.7059]);

% Beautify plot
xlim([0.5, 1.5]);
xticks(x_center);
xticklabels({'diff\_m'});
set(gca, 'FontName', 'Arial', 'FontSize', 5, 'FontWeight', 'bold');
ylabel('Values', 'FontName', 'Arial', 'FontWeight', 'bold', 'FontSize', 7);
title('CFA:Midbrain', 'FontName', 'Arial', 'FontWeight', 'bold', 'FontSize', 7);
grid on;
hold off;

% Save figure %
saveas(gcf, fullfile(filepath1, 'Midbrain_stat'), 'epsc');
% Thalamus not needed, already done %


%% Thalamus/PONS Regression %%


PoG = [dataTable_sel_ipsi_GFP(57:63, 2); dataTable_sel_contra_GFP(57:63, 2); dataTable_sel_ipsi_GFP(57:63, 3); dataTable_sel_contra_GFP(57:63, 3)];
PoR = [dataTable_sel_ipsi_RFP(57:63, 2); dataTable_sel_contra_RFP(57:63, 2); dataTable_sel_ipsi_RFP(57:63, 3); dataTable_sel_contra_RFP(57:63, 3)];


pu = [1, 0.5, 0]; 

% Thalamus %
figure('Units', 'inches', 'Position', [100,100,3,2]);
% % Fit linear regression model: y = m*x + b
% pt = polyfit(Xt, Yt, 1);  % p(1) = slope, p(2) = intercept
% 
% % Generate the regression line for plotting
% Y_fitt = polyval(pt, Xt);
% 
% % Plot the regression line
% plot(Xt, Y_fitt, 'c-', 'LineWidth', 1);  % Regression line
% hold on
% 
% scatter(Xt, Yt, 10,'MarkerEdgeColor', 'c', 'MarkerFaceColor', 'c','LineWidth', 1);
% hold on;
% % Add labels from the list of names to each data point
% for i = 1:length(Xt)
%     text(Xt(i), Yt(i), Zt{i}, ...
%         'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right','FontName', 'Arial', 'FontSize', 5, 'FontWeight', 'bold', 'Color', 'k');
% end
% hold on
% 
% % Add labels from the list of names to each data point
% for i = 1:length(Xt)
%     % Convert Zc to string format using num2str
%     text(Xt(i), Yt(i), Zt{i}, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontName', 'Arial', 'FontSize', 5, 'FontWeight', 'bold', 'Color', 'k');
% end
% hold on;


% PONS %

Xpo = PoG;  % Replaced PG with CoG
Ypo = PoR;  % Replaced PR with CoR
Zpo = [dataTable_Area_names([57:63]); dataTable_Area_names([57:63]); dataTable_Area_names([57:63]); dataTable_Area_names([57:63])];

% saving %
writematrix(PoG, fullfile(filepath1, 'PONS_PFC.xlsx'));
writematrix(PoR, fullfile(filepath1, 'PONS_PFC.xlsx'));
writecell(Zpo, fullfile(filepath1, 'PONS_Areas_PFC.xlsx'));

% Fit linear regression model: y = m*x + b
po = polyfit(Xpo, Ypo, 1);  % po(1) = slope, po(2) = intercept

% Generate the regression line for plotting
Y_fith_po = polyval(po, Xpo);

% Plot the regression line with new color
% Define the color 'pu' as RGB values
pu = [1, 0.5, 0];
plot(Xpo, Y_fith_po, 'Color', pu, 'LineWidth', 1);  % Regression line
hold on;

% Number of parts to split into (4 parts)
numParts = 4;

% Loop through each 1/4th segment of Xpo and Ypo
for i = 0:numParts-1
% Define the start and end indices for each 1/4th of the data
startIdx = floor(i * length(Xpo) / numParts) + 1;
endIdx = floor((i + 1) * length(Xpo) / numParts);
if i == 0 || i == 2
% Plot the current 1/4th segment with the 'pu' color
scatter(Xpo(startIdx:endIdx), Ypo(startIdx:endIdx), 20, 'MarkerEdgeColor', pu, 'MarkerFaceColor', pu, 'LineWidth', 1);
hold on;
else
scatter(Xpo(startIdx:endIdx), Ypo(startIdx:endIdx), 20, 'MarkerEdgeColor', pu, 'MarkerFaceColor', 'w', 'LineWidth', 1);
hold on;
end
end


% Add labels from the list of names to each data point
for i = 1:length(Xpo)
text(Xpo(i), Ypo(i), Zpo{i}, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontName', 'Arial', 'FontSize', 5, 'FontWeight', 'bold', 'Color', 'k');
end

% Add labels and title
set(gca, 'FontName', 'Arial', 'FontSize', 5, 'FontWeight', 'bold');
xlabel('Direct Projections Fraction Occupation [%]', 'FontName', 'Arial', 'FontWeight', 'bold', 'FontSize', 7);
ylabel('Indirect Projections Fraction Occupation [%]', 'FontName', 'Arial', 'FontWeight', 'bold', 'FontSize', 7);
title(['CFA Projections: \color{cyan}Thalamus \color[rgb]{black} vs \color[rgb]{0.4940, 0.1840, 0.5560}PONS'], ...
'FontName', 'Arial', 'FontWeight', 'bold', 'FontSize', 7);

hold off;

% save figure %
saveas(gcf,fullfile(filepath1, 'Th_vs_Pons_reg'), 'epsc');

%% Individually Thalamus/PONS tests %%
% Normalize all data together
all_data = [TG; TR; PoG; PoR];

mean_val = mean(all_data);  % Global mean
std_val = std(all_data);    % Global standard deviation

TG_norm = (TG - mean_val) / std_val;
TR_norm = (TR - mean_val) / std_val;
PoG_norm = (PoG - mean_val) / std_val;
PoR_norm = (PoR - mean_val) / std_val;
%% raincloud plot %%
raincloud_plot(TG, TR, PoG, PoR, 'TG', 'TR', 'PoG', 'PoR');

% Difference %
diff_po = [PoR_norm - PoG_norm];
diff_th = [TR_norm - TG_norm];
%% PONS %%
[po, ~] = kstest((diff_po - mean(diff_po)) / std(diff_po)); % Kolmogorov-Smirnov test for normality
if po == 0
disp('Data follows normal distribution (po > 0.05).');
else
disp('Data does not follow normal distribution (po <= 0.05).');
end

% % Perform Wilcoxon Signed-Rank Test for paired data
% [~, p_val] = signrank(PoR_norm, PoG_norm);  % Second output is p-value
% disp(['Wilcoxon Signed-Rank Test p-value: ', num2str(p_val)]);  % Display p-value


% Perform one-sample t-test
[fo_t, p_val] = ttest(diff_po, 0);  % Testing if the mean of diff_ho is significantly different from 0

disp(['t-statistic: ', num2str(fo_t)]);
disp(['p-value: ', num2str(p_val)]);

% Calculate the median and quartiles
median_val = median(diff_po);
q1 = prctile(diff_po, 25); % First quartile
q3 = prctile(diff_po, 75); % Third quartile

% Create the plot
figure('Units', 'inches', 'Position', [100,100,3,2]);
hold on;
pu = [1, 0.5, 0];  % Purple color for the plot

% Plot the box manually
x_center = 1; % X position for diff_po
box_width = 0.1;
rectangle('Position', [x_center - box_width / 2, q1, box_width, q3 - q1], ...
'EdgeColor', pu, 'FaceColor', [1, 1, 1]); % white box for IQR
plot([x_center - box_width / 2, x_center + box_width / 2], [median_val, median_val], ...
'-', 'LineWidth', 1, 'Color', pu); % Line for median

% Overlay individual data points in pu
scatter(x_center * ones(size(diff_po)), diff_po, 10, 'MarkerFaceColor', pu, ...
'MarkerEdgeColor', pu, 'LineWidth', 1);

% Add p-value as text
y_pos = q3 + 0.2 * range(diff_po); % Position for p-value text
text(x_center - box_width, y_pos, sprintf('p = %.5f', p_val), ...
'HorizontalAlignment', 'center', 'FontSize', 7, 'FontWeight', 'bold', 'Color', pu);

% Beautify plot
xlim([0.5, 1.5]);
xticks(x_center);
xticklabels({'diff\_po'});
set(gca, 'FontName', 'Arial', 'FontSize', 5, 'FontWeight', 'bold');
ylabel('Values','FontName', 'Arial', 'FontWeight', 'bold','FontSize', 7);
title('PFC:PONS','FontName', 'Arial', 'FontWeight', 'bold','FontSize', 7);
grid on;
hold off;

% Save figure
saveas(gcf, fullfile(filepath1, 'PONS_stat'), 'epsc');

%% Thalamus/Medulla Regression %%
MeG = [dataTable_sel_ipsi_GFP([64:end],2);dataTable_sel_contra_GFP([64:end],2);dataTable_sel_ipsi_GFP([64:end],3);dataTable_sel_contra_GFP([64:end],3)];
MeR = [dataTable_sel_ipsi_RFP([64:end],2);dataTable_sel_contra_RFP([64:end],2);dataTable_sel_ipsi_RFP([64:end],3);dataTable_sel_contra_RFP([64:end],3)];

% Thalamus %
figure('Units', 'inches', 'Position', [100,100,3,2]);
% % Fit linear regression model: y = m*x + b
% pt = polyfit(Xt, Yt, 1);  % p(1) = slope, p(2) = intercept
% 
% % Generate the regression line for plotting
% Y_fitt = polyval(pt, Xt);
% 
% % Plot the regression line
% plot(Xt, Y_fitt, 'c-', 'LineWidth', 1);  % Regression line
% hold on
% 
% scatter(Xt, Yt, 10,'MarkerEdgeColor', 'c', 'MarkerFaceColor', 'c','LineWidth', 1);
% hold on;
% % Add labels from the list of names to each data point
% for i = 1:length(Xt)
%     text(Xt(i), Yt(i), Zt{i}, ...
%         'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right','FontName', 'Arial', 'FontSize', 5, 'FontWeight', 'bold', 'Color', 'k');
% end
% hold on

% Medulla %

Xme = MeG;
Yme = MeR;
Zme = [dataTable_Area_names([64:end]); dataTable_Area_names([64:end]); dataTable_Area_names([64:end]); dataTable_Area_names([64:end])];

% saving %
writematrix(MeG, fullfile(filepath1, 'Me_GFP_PFC.xlsx'));
writematrix(MeR, fullfile(filepath1, 'Me_RFP_PFC.xlsx'));
writecell(Zme, fullfile(filepath1, 'Me_Areas_PFC.xlsx'));


% Fit linear regression model: y = m*x + b
me_fit = polyfit(Xme, Yme, 1);  % me_fit(1) = slope, me_fit(2) = intercept
% Define the color 'violet' as RGB values
violet = [0.5804 0.0000 0.8275];

% Generate the regression line for plotting
Y_fitme = polyval(me_fit, Xme);

% Plot the regression line
plot(Xme, Y_fitme, '-', 'Color', violet, 'LineWidth', 1);  % Regression line
hold on
% Version 6: Xme, Yme, color = violet (RGB values)

% Loop through each 1/4th segment of Xme and Yme
for i = 0:numParts-1
% Define the start and end indices for each 1/4th of the data
startIdx = floor(i * length(Xme) / numParts) + 1;
endIdx = floor((i + 1) * length(Xme) / numParts);
if i == 0 || i == 2
% Plot the current 1/4th segment with violet color
scatter(Xme(startIdx:endIdx), Yme(startIdx:endIdx), 20, 'MarkerEdgeColor', violet, 'MarkerFaceColor', violet, 'LineWidth', 1);
hold on;
else
scatter(Xme(startIdx:endIdx), Yme(startIdx:endIdx), 20, 'MarkerEdgeColor', violet, 'MarkerFaceColor', 'w', 'LineWidth', 1);
hold on;
end
end

% Add labels from the list of names to each data point
for i = 1:length(Xme)
text(Xme(i), Yme(i), Zme{i}, ...
'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', ...
'FontName', 'Arial', 'FontSize', 5, 'FontWeight', 'bold', 'Color', 'k');  % 'k' is black
end

% Add labels and title
set(gca, 'FontName', 'Arial', 'FontSize', 5, 'FontWeight', 'bold');
xlabel('Direct Projections Fraction Occupation [%]', 'FontName', 'Arial', 'FontWeight', 'bold', 'FontSize', 7);
ylabel('Indirect Projections Fraction Occupation [%]', 'FontName', 'Arial', 'FontWeight', 'bold', 'FontSize', 7);
title('CFA Projections: \color{cyan}Thalamus \color{black}vs \color[rgb]{0.5804, 0.0000, 0.8275}Medulla', 'FontName', 'Arial', 'FontWeight', 'bold', 'FontSize', 7);
hold off;
% Save figure %
saveas(gcf, fullfile(filepath1, 'Th_vs_Medulla_reg'), 'epsc');
%% Statistical tests including ipsi and contra projections %%

% Normalize all data together
all_data = [TG; TR; MeG; MeR];

mean_val = mean(all_data);  % Global mean
std_val = std(all_data);    % Global standard deviation

MeG_norm = (MeG - mean_val) / std_val;
MeR_norm = (MeR - mean_val) / std_val;
TG_norm = (TG - mean_val) / std_val;
TR_norm = (TR - mean_val) / std_val;
%% raincloud plot %%
raincloud_plot(TG, TR, MeG, MeR, 'TG', 'IR', 'MeG', 'MeR');
% Difference %
diff_me = [MeR_norm - MeG_norm];
diff_th = [TR_norm - TG_norm];
%% Individually Thalamus and Medulla tests %%
% Medulla %
[h, ~] = kstest((diff_me - mean(diff_me)) / std(diff_me)); % Kolmogorov-Smirnov test for normality
if h == 0
disp('Data follows normal distribution (p > 0.05).');
else
disp('Data does not follow normal distribution (p <= 0.05).');
end

% % Perform Wilcoxon Signed-Rank Test
% me_p = signrank(diff_me);
% disp(['p-value: ', num2str(me_p)]);


% Perform one-sample t-test
[fo_t, me_p] = ttest(diff_me, 0);  % Testing if the mean of diff_ho is significantly different from 0

disp(['t-statistic: ', num2str(fo_t)]);
disp(['p-value: ', num2str(me_p)]);



% Calculate the median and quartiles
median_val = median(diff_me);
q1 = prctile(diff_me, 25); % First quartile
q3 = prctile(diff_me, 75); % Third quartile

% Create the plot
figure('Units', 'inches', 'Position', [100,100,3,2]);
hold on;

% Plot the box manually
x_center = 1; % X position for diff_me
box_width = 0.1;
rectangle('Position', [x_center - box_width / 2, q1, box_width, q3 - q1], ...
'EdgeColor', [0.5412 0.1686 0.8863], 'FaceColor', [1, 1, 1]); % Violet box for IQR
plot([x_center - box_width / 2, x_center + box_width / 2], [median_val, median_val], ...
'-', 'Color', [0.5412 0.1686 0.8863], 'LineWidth', 1); % Line for median

% Overlay individual data points in violet
scatter(x_center * ones(size(diff_me)), diff_me, 10, 'MarkerFaceColor', [0.5412 0.1686 0.8863], ...
'MarkerEdgeColor', [0.5412 0.1686 0.8863], 'LineWidth', 1);

% Add p-value
y_pos = q3 + 0.2 * range(diff_me); % Position for p-value
text(x_center - box_width, y_pos, sprintf('p = %.5f', me_p), ...
'HorizontalAlignment', 'center', 'FontSize', 7, 'FontWeight', 'bold', 'Color', [0.5412 0.1686 0.8863]);

% Beautify plot
xlim([0.5, 1.5]);
xticks(x_center);
xticklabels({'diff\_me'});
set(gca, 'FontName', 'Arial', 'FontSize', 5, 'FontWeight', 'bold');
ylabel('Values', 'FontName', 'Arial', 'FontWeight', 'bold', 'FontSize', 7);
title('CFA:Medulla', 'FontName', 'Arial', 'FontWeight', 'bold', 'FontSize', 7);
grid on;
hold off;

% Save figure %
saveas(gcf, fullfile(filepath1, 'Medulla_stat'), 'epsc');

%% FO vs HO regression %%

% First order thalamus %

FTG = [dataTable_sel_ipsi_GFP([34:35,39:40],2);dataTable_sel_contra_GFP([34:35,39:40],2);dataTable_sel_ipsi_GFP([34:35,39:40],3);dataTable_sel_contra_GFP([34:35,39:40],3)];
FTR = [dataTable_sel_ipsi_RFP([34:35,39:40],2);dataTable_sel_contra_RFP([34:35,39:40],2);dataTable_sel_ipsi_RFP([34:35,39:40],3);dataTable_sel_contra_RFP([34:35,39:40],3)];

Xfot = FTG;
Yfot = FTR;
Zfot = [dataTable_Area_names([34:35,39:40]); dataTable_Area_names([34:35,39:40]); dataTable_Area_names([34:35,39:40]); dataTable_Area_names([34:35,39:40])];

% saving %
writematrix(FTG, fullfile(filepath1, 'FT_GFP_PFC.xlsx'));
writematrix(FTR, fullfile(filepath1, 'FT_RFP_PFC.xlsx'));
writecell(Zfot, fullfile(filepath1, 'FT_Areas_PFC.xlsx'));

% Fit linear regression model: y = m*x + b
fot_fit = polyfit(Xfot, Yfot, 1);  % fot_fit(1) = slope, fot_fit(2) = intercept
% Define the color 'fo' as RGB values
fo = [0 0.4470 0.7410];

% Generate the regression line for plotting
Y_fitfot = polyval(fot_fit, Xfot);
figure('Units', 'inches', 'Position', [100,100,3,2]);
% Plot the regression line
plot(Xfot, Y_fitfot, '-', 'Color', fo, 'LineWidth', 1);  % Regression line
hold on
% Version 6: Xfot, Yfot, color = fo (RGB values)

% Loop through each 1/4th segment of Xfot and Yfot
for i = 0:numParts-1
% Define the start and end indices for each 1/4th of the data
startIdx = floor(i * length(Xfot) / numParts) + 1;
endIdx = floor((i + 1) * length(Xfot) / numParts);
if i == 0 || i == 2
% Plot the current 1/4th segment with fo color
scatter(Xfot(startIdx:endIdx), Yfot(startIdx:endIdx), 20, 'MarkerEdgeColor', fo, 'MarkerFaceColor', fo, 'LineWidth', 1);
hold on;
else
scatter(Xfot(startIdx:endIdx), Yfot(startIdx:endIdx), 20, 'MarkerEdgeColor', fo, 'MarkerFaceColor', 'w', 'LineWidth', 1);
hold on;
end
end

% Add labels from the list of names to each data point
for i = 1:length(Xfot)
text(Xfot(i), Yfot(i), Zfot{i}, ...
'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', ...
'FontName', 'Arial', 'FontSize', 5, 'FontWeight', 'bold', 'Color', 'k');  % 'k' is black
end

% Higher order thalamus %
HTG = [dataTable_sel_ipsi_GFP([32:33,36:38,41:47],2);dataTable_sel_contra_GFP([32:33,36:38,41:47],2);dataTable_sel_ipsi_GFP([32:33,36:38,41:47],3);dataTable_sel_contra_GFP([32:33,36:38,41:47],3)];
HTR = [dataTable_sel_ipsi_RFP([32:33,36:38,41:47],2);dataTable_sel_contra_RFP([32:33,36:38,41:47],2);dataTable_sel_ipsi_RFP([32:33,36:38,41:47],3);dataTable_sel_contra_RFP([32:33,36:38,41:47],3)];

Xhot = HTG;
Yhot = HTR;
Zhot = [dataTable_Area_names([32:33,36:38,41:47]); dataTable_Area_names([32:33,36:38,41:47]); dataTable_Area_names([32:33,36:38,41:47]); dataTable_Area_names([32:33,36:38,41:47])];

% saving %
writematrix(HTG, fullfile(filepath1, 'HT_GFP_PFC.xlsx'));
writematrix(HTR, fullfile(filepath1, 'HT_RFP_PFC.xlsx'));
writecell(Zhot, fullfile(filepath1, 'HT_Areas_PFC.xlsx'));

% Fit linear regression model: y = m*x + b
hot_fit = polyfit(Xhot, Yhot, 1);  % hot_fit(1) = slope, hot_fit(2) = intercept
% Define the color 'ho' as RGB values
ho = [0, 0.5, 0];

% Generate the regression line for plotting
Y_fithot = polyval(hot_fit, Xhot);

% Plot the regression line
plot(Xhot, Y_fithot, '-', 'Color', ho, 'LineWidth', 1);  % Regression line
hold on
% Version 6: Xhot, Yhot, color = ho (RGB values)

% Loop through each 1/4th segment of Xhot and Yhot
for i = 0:numParts-1
% Define the start and end indices for each 1/4th of the data
startIdx = floor(i * length(Xhot) / numParts) + 1;
endIdx = floor((i + 1) * length(Xhot) / numParts);
if i == 0 || i == 2
% Plot the current 1/4th segment with ho color
scatter(Xhot(startIdx:endIdx), Yhot(startIdx:endIdx), 20, 'MarkerEdgeColor', ho, 'MarkerFaceColor', ho, 'LineWidth', 1);
hold on;
else
scatter(Xhot(startIdx:endIdx), Yhot(startIdx:endIdx), 20, 'MarkerEdgeColor', ho, 'MarkerFaceColor', 'w', 'LineWidth', 1);
hold on;
end
end

% Add labels from the list of names to each data point
for i = 1:length(Xhot)
text(Xhot(i), Yhot(i), Zhot{i}, ...
'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', ...
'FontName', 'Arial', 'FontSize', 5, 'FontWeight', 'bold', 'Color', 'k');  % 'k' is black
end

% Add labels and title
set(gca, 'FontName', 'Arial', 'FontSize', 5, 'FontWeight', 'bold');
xlabel('Direct Projections Fraction Occupation [%]', 'FontName', 'Arial', 'FontWeight', 'bold', 'FontSize', 7);
ylabel('Indirect Projections Fraction Occupation [%]', 'FontName', 'Arial', 'FontWeight', 'bold', 'FontSize', 7);
title('CFA Projections: \color[rgb]{0, 0.4470, 0.7410}FO Thalamus \color{black}vs \color[rgb]{0.4940, 0.1840, 0.5560}HO Thalamus', ...
'FontName', 'Arial', 'FontWeight', 'bold', 'FontSize', 7);
hold off;
% Save figure %
saveas(gcf, fullfile(filepath1, 'FO_vs_HO_Thalamus_reg'), 'epsc');

%% Statistical tests including ipsi and contra projections %%

% Normalize all data together
all_data = [FTG; FTR; HTG; HTR];

mean_val = mean(all_data);  % Global mean
std_val = std(all_data);    % Global standard deviation

FTG_norm = (FTG - mean_val) / std_val;
FTR_norm = (FTR - mean_val) / std_val;
HTG_norm = (HTG - mean_val) / std_val;
HTR_norm = (HTR - mean_val) / std_val;
%% raincloud plot %%
raincloud_plot(FTG, FTR, HTG, HTR, 'FTG', 'FTR', 'HTG', 'HTR');
% Difference %
diff_fo = [FTR_norm - FTG_norm];
diff_ho = [HTR_norm - HTG_norm];

%% Individually FO Thalamus and HO Thalamus tests %%
% First order thalamus %
% [h, ~] = kstest((diff_fo - mean(diff_fo)) / std(diff_fo)); % Kolmogorov-Smirnov test for normality
% if h == 0
%     disp('Data follows normal distribution (p > 0.05).');
% else
%     disp('Data does not follow normal distribution (p <= 0.05).');
% end

% % Shapiro–Wilk test
[h, pValue, W] = swtest(diff_fo, 0.05);  % 0.05 is the significance level

if h == 0
disp('Data follows normal distribution (p > 0.05)');
else
disp('Data does NOT follow normal distribution (p <= 0.05)');
end

% Perform one-sample t-test
[fo_t, fo_p] = ttest(diff_fo, 0);  % Testing if the mean of diff_ho is significantly different from 0

disp(['t-statistic: ', num2str(fo_t)]);
disp(['p-value: ', num2str(fo_p)]);


% Perform Wilcoxon Signed-Rank Test
fo_p = signrank(diff_fo);
disp(['p-value: ', num2str(fo_p)]);
% Calculate the median and quartiles
median_val = median(diff_fo);
q1 = prctile(diff_fo, 25); % First quartile
q3 = prctile(diff_fo, 75); % Third quartile

% Create the plot
figure('Units', 'inches', 'Position', [100,100,3,2]);
hold on;

% Plot the box manually
x_center = 1; % X position for diff_fo
box_width = 0.1;
rectangle('Position', [x_center - box_width / 2, q1, box_width, q3 - q1], ...
'EdgeColor', [0 0.4470 0.7410], 'FaceColor', [1, 1, 1]); % fo color box for IQR
plot([x_center - box_width / 2, x_center + box_width / 2], [median_val, median_val], ...
'-', 'Color', [0 0.4470 0.7410], 'LineWidth', 1); % Line for median

% Overlay individual data points in fo color
scatter(x_center * ones(size(diff_fo)), diff_fo, 10, 'MarkerFaceColor', [0 0.4470 0.7410], ...
'MarkerEdgeColor', [0 0.4470 0.7410], 'LineWidth', 1);

% Add p-value
y_pos = q3 + 0.2 * range(diff_fo); % Position for p-value
text(x_center - box_width, y_pos, sprintf('p = %.5f', fo_p), ...
'HorizontalAlignment', 'center', 'FontSize', 5, 'FontWeight', 'bold', 'Color', [0 0.4470 0.7410]);

% Beautify plot
xlim([0.5, 1.5]);
xticks(x_center);
xticklabels({'diff\_fo'});
set(gca, 'FontName', 'Arial', 'FontSize', 5, 'FontWeight', 'bold');
ylabel('Values', 'FontName', 'Arial', 'FontWeight', 'bold', 'FontSize', 7);
title('CFA:FO Thalamus', 'FontName', 'Arial', 'FontWeight', 'bold', 'FontSize', 7);
grid on;
hold off;

% Save figure %
saveas(gcf, fullfile(filepath1, 'FO_stat'), 'epsc');

%% Higher order thalamus %%
[h, ~] = kstest((diff_ho - mean(diff_ho)) / std(diff_ho)); % Kolmogorov-Smirnov test for normality
if h == 0
disp('Data follows normal distribution (p > 0.05).');
else
disp('Data does not follow normal distribution (p <= 0.05).');
end

% % Shapiro–Wilk test
[h, pValue, W] = swtest(diff_ho, 0.05);  % 0.05 is the significance level

if h == 0
disp('Data follows normal distribution (p > 0.05)');
else
disp('Data does NOT follow normal distribution (p <= 0.05)');

end
% Perform Wilcoxon Signed-Rank Test
ho_p = signrank(diff_ho);
disp(['p-value: ', num2str(ho_p)]);
% Calculate the median and quartiles
median_val = median(diff_ho);
q1 = prctile(diff_ho, 25); % First quartile
q3 = prctile(diff_ho, 75); % Third quartile

% Create the plot
figure('Units', 'inches', 'Position', [100,100,3,2]);
hold on;

% Plot the box manually
x_center = 1; % X position for diff_ho
box_width = 0.1;
rectangle('Position', [x_center - box_width / 2, q1, box_width, q3 - q1], ...
'EdgeColor', [0 0.5 0], 'FaceColor', [1, 1, 1]); % ho color box for IQR
plot([x_center - box_width / 2, x_center + box_width / 2], [median_val, median_val], ...
'-', 'Color', [0 0.5 0], 'LineWidth', 1); % Line for median

% Overlay individual data points in ho color
scatter(x_center * ones(size(diff_ho)), diff_ho, 10, 'MarkerFaceColor', [0 0.5 0], ...
'MarkerEdgeColor', [0 0.5 0], 'LineWidth', 1);

% Add p-value
y_pos = q3 + 0.2 * range(diff_ho); % Position for p-value
text(x_center - box_width, y_pos, sprintf('p = %.5f', ho_p), ...
'HorizontalAlignment', 'center', 'FontSize', 5, 'FontWeight', 'bold', 'Color', [0 0.5 0]);

% Beautify plot
xlim([0.5, 1.5]);
xticks(x_center);
xticklabels({'diff\_ho'});
set(gca, 'FontName', 'Arial', 'FontSize', 5, 'FontWeight', 'bold');
ylabel('Values', 'FontName', 'Arial', 'FontWeight', 'bold', 'FontSize', 7);
title('CFA:HO Thalamus', 'FontName', 'Arial', 'FontWeight', 'bold', 'FontSize', 7);
grid on;
hold off;

% Save figure %
saveas(gcf, fullfile(filepath1, 'HO_stat'), 'epsc');
%% Comparing differences in FO with HO %%
% Perform Mann-Whitney U Test (Wilcoxon rank-sum test)
[pFH, h] = ranksum(diff_fo, diff_ho);

disp(['Mann-Whitney U Test p-value: ', num2str(pFH)]);

% Data for plotting
group_data = {diff_fo, diff_ho}; % Updated group data
group_names = {'$\delta$ FO', '$\delta$ HO'}; % Updated group labels

% Specify colors for box plots and individual points
box_colors = [0, 0.4470, 0.7410; 0, 0.5, 0];  % Colors for each box (Group 1: fo, Group 2: ho)
point_colors = [0, 0.4470, 0.7410; 0, 0.5, 0];  % Colors for individual points

% Check if group_data is a cell array and contains numeric vectors
if iscell(group_data) && all(cellfun(@isnumeric, group_data))
% Ensure that group_data{1} and group_data{2} are vectors (1D arrays)
group_data{1} = group_data{1}(:);  % Convert to column vector if not already
group_data{2} = group_data{2}(:);  % Convert to column vector if not already

% Combine the group data into a single vector for the boxplot
combined_data = [group_data{1}; group_data{2}]; % Concatenate vertically (column-wise)
group_ids = [ones(1, length(group_data{1})), 2*ones(1, length(group_data{2}))]; % Group IDs for the boxplot

% Create a figure
figure('Units', 'inches', 'Position', [100,100,3,2]);
hold on;

% Box plot for both groups
h = boxplot(combined_data, group_ids, 'Labels', group_names, 'Whisker', 1.5, 'Colors', [0, 0, 0]);

% Set font properties for axes
set(gca, 'FontName', 'Arial', 'FontSize', 5, 'FontWeight', 'bold');  % Font settings for the axes

% Customize box colors for each group
for i = 1:2
% Set box colors
set(h(:, i), 'Color', box_colors(i, :), 'LineWidth', 1); % Black outline for boxes
set(h(5, i), 'Color', box_colors(i, :));   % Color for box (5th row is the box color)
set(h(6, i), 'Color', box_colors(i, :));   % Color for the median line (6th row is median line)
end

% Overlay individual data points for each group
for i = 1:length(group_data)
scatter(ones(1, length(group_data{i})) * i, group_data{i}, 10, 'MarkerFaceColor', point_colors(i, :), 'MarkerEdgeColor', point_colors(i, :), 'LineWidth', 1);
end

% Calculate the position for displaying the p-value (above the box plot)
y_pos = max(combined_data) + 0.1 * range(combined_data); % Position for p-value above the boxes

% Display the p-value as stars
if pFH < 0.001
p_str = 'p < 0.001';
elseif pFH < 0.01
p_str = 'p < 0.01';
elseif pFH < 0.05
p_str = 'p < 0.05';
else
p_str = sprintf('p = %.5f', pFH); % Show exact p-value if it's above 0.05
end

% Add p-value text to the plot
text(1.5, y_pos, p_str, 'HorizontalAlignment', 'center', 'FontSize', 7, 'FontWeight', 'bold', 'Color', 'black');

% Beautify plot
% Set title with LaTeX interpreter
title('Comparison of Direct and Indirect Projections', 'FontName', 'Arial', 'FontSize', 7, 'FontWeight', 'bold', 'Interpreter', 'latex');

% Set ylabel with LaTeX interpreter
ylabel('Fraction Occupation', 'FontName', 'Arial', 'FontSize', 7, 'FontWeight', 'bold', 'Interpreter', 'latex');

% Set xtick labels with LaTeX interpreter
xticklabels({'$\delta$ FO', '$\delta$ HO'});  % Using LaTeX for delta symbol
% % % set(gca, 'TickLabelInterpreter', 'latex');  % Apply LaTeX to the tick labels

grid on;
hold off;
end

% Save figure
saveas(gcf, fullfile(filepath1, 'FO_vs_HO_stat'), 'epsc');
%% Thalamus/Basal Ganglia [GG] Regression %%
GG = [dataTable_sel_ipsi_GFP([27,30:31,48,50,55],2);dataTable_sel_contra_GFP([27,30:31,48,50,55],2);dataTable_sel_ipsi_GFP([27,30:31,48,50,55],3);dataTable_sel_contra_GFP([27,30:31,48,50,55],3)];
GR = [dataTable_sel_ipsi_RFP([27,30:31,48,50,55],2);dataTable_sel_contra_RFP([27,30:31,48,50,55],2);dataTable_sel_ipsi_RFP([27,30:31,48,50,55],3);dataTable_sel_contra_RFP([27,30:31,48,50,55],3)];

% Thalamus %
figure('Units', 'inches', 'Position', [100,100,3,2]);
% % Fit linear regression model: y = m*x + b
% pt = polyfit(Xt, Yt, 1);  % p(1) = slope, p(2) = intercept
% 
% % Generate the regression line for plotting
% Y_fitt = polyval(pt, Xt);
% 
% % Plot the regression line
% plot(Xt, Y_fitt, 'c-', 'LineWidth', 1);  % Regression line
% hold on
% % Number of parts to split into (4 parts)
% numParts = 4;
% % Loop through each 1/4th segment of Xt and Yt
% for i = 0:numParts-1
%     % Define the start and end indices for each 1/4th of the data
%     startIdx = floor(i * length(Xt) / numParts) + 1;
%     endIdx = floor((i + 1) * length(Xt) / numParts);
%     if i == 0 || i == 2
%         % Plot the current 1/4th segment with cyan color
%         scatter(Xt(startIdx:endIdx), Yt(startIdx:endIdx), 20, 'MarkerEdgeColor', 'cyan', 'MarkerFaceColor', 'cyan', 'LineWidth', 1);
%         hold on;
%     else
%         scatter(Xt(startIdx:endIdx), Yt(startIdx:endIdx), 20, 'MarkerEdgeColor', 'cyan', 'MarkerFaceColor', 'w', 'LineWidth', 1);
%         hold on;
%     end
% end
% % Add labels from the list of names to each data point
% for i = 1:length(Xt)
%     text(Xt(i), Yt(i), Zt{i}, ...
%         'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right','FontName', 'Arial', 'FontWeight', 'bold', 'FontSize', 5);
% end
% hold on


% Basal Ganglia %

Xbg = GG;
Ybg = GR;
Zbg = [dataTable_Area_names([27,30:31,48,50,55]); dataTable_Area_names([27,30:31,48,50,55]); dataTable_Area_names([27,30:31,48,50,55]); dataTable_Area_names([27,30:31,48,50,55])];

% saving %
writematrix(GG, fullfile(filepath1, 'BG_GFP_PFC.xlsx'));
writematrix(GR, fullfile(filepath1, 'BG_RFP_PFC.xlsx'));
writecell(Zbg, fullfile(filepath1, 'BG_Areas_PFC.xlsx'));

% Define the color 'mag' as RGB values
mag = [1.0000 0.0000 0000];
% Fit linear regression model: y = m*x + b
bg_fit = polyfit(Xbg, Ybg, 1);  % bg_fit(1) = slope, bg_fit(2) = intercept

% Generate the regression line for plotting
Y_fitbg = polyval(bg_fit, Xbg);

% Plot the regression line
plot(Xbg, Y_fitbg, 'm-', 'LineWidth', 1);  % Regression line
hold on
% Loop through each 1/4th segment of Xbg and Ybg
for i = 0:numParts-1
% Define the start and end indices for each 1/4th of the data
startIdx = floor(i * length(Xbg) / numParts) + 1;
endIdx = floor((i + 1) * length(Xbg) / numParts);
if i == 0 || i == 2
% Plot the current 1/4th segment with magenta color
scatter(Xbg(startIdx:endIdx), Ybg(startIdx:endIdx), 20, 'MarkerEdgeColor', mag, 'MarkerFaceColor', mag, 'LineWidth', 1);
hold on;
else
scatter(Xbg(startIdx:endIdx), Ybg(startIdx:endIdx), 20, 'MarkerEdgeColor', mag, 'MarkerFaceColor', 'w', 'LineWidth', 1);
hold on;
end
end

% Add labels from the list of names to each data point

for i = 1:length(Xbg)
text(Xbg(i), Ybg(i), Zbg{i}, ...
'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', ...
'FontName', 'Arial', 'FontSize', 5, 'FontWeight', 'bold', 'Color', 'k');  % 'k' is black
end

% Add labels and title
set(gca, 'FontName', 'Arial', 'FontSize', 5, 'FontWeight', 'bold');
xlabel('Direct Projections Fraction Occupation [%]', 'FontName', 'Arial', 'FontWeight', 'bold', 'FontSize', 7);
ylabel('Indirect Projections Fraction Occupation [%]', 'FontName', 'Arial', 'FontWeight', 'bold', 'FontSize', 7);
title('CFA Projections: \color{cyan}Thalamus \color{black}vs \color{magenta}BasGanG', 'FontName', 'Arial', 'FontWeight', 'bold', 'FontSize', 7);
hold off;
% Save figure %
saveas(gcf, fullfile(filepath1, 'Th_vs_BasGanG_reg'), 'epsc');


%% Statistical tests including ipsi and contra projections %%

% Normalize all data together
all_data = [TG; TR; GG; GR];

mean_val = mean(all_data);  % Global mean
std_val = std(all_data);    % Global standard deviation

GG_norm = (GG - mean_val) / std_val;
GR_norm = (GR - mean_val) / std_val;
TG_norm = (TG - mean_val) / std_val;
TR_norm = (TR - mean_val) / std_val;
%% raincloud plot %%
raincloud_plot(TG, TR, GG, GR, 'TG', 'TR', 'GG', 'GR');
% Difference %

diff_bg = [GR_norm - GG_norm];
diff_th = [TR_norm - TG_norm];
%% Individually Thalamus and BasGang tests %%
% BasGang %
[h, ~] = kstest((diff_bg - mean(diff_bg)) / std(diff_bg)); % Kolmogorov-Smirnov test for normality
if h == 0
disp('Data follows normal distribution (p > 0.05).');
else
disp('Data does not follow normal distribution (p <= 0.05).');
end
%
% % Perform Wilcoxon Signed-Rank Test
% bg_p = signrank(diff_bg);
% disp(['p-value: ', num2str(bg_p)]);

% Perform One-Sample t-test
[~, bg_p] = ttest(diff_bg, 0);
disp(['p-value: ', num2str(bg_p)]);


% Calculate the median and quartiles
median_val = median(diff_bg);
q1 = prctile(diff_bg, 25); % First quartile
q3 = prctile(diff_bg, 75); % Third quartile

% Create the plot
figure('Units', 'inches', 'Position', [100,100,3,2]);
hold on;

% Plot the box manually
x_center = 1; % X position for diff_bg
box_width = 0.1;
rectangle('Position', [x_center - box_width / 2, q1, box_width, q3 - q1], ...
'EdgeColor', 'red', 'FaceColor', [1, 1, 1]); % red box for IQR
plot([x_center - box_width / 2, x_center + box_width / 2], [median_val, median_val], ...
'r-', 'LineWidth', 1); % Line for median

% Overlay individual data points in red
scatter(x_center * ones(size(diff_bg)), diff_bg, 10, 'MarkerFaceColor', 'red', ...
'MarkerEdgeColor', 'red', 'LineWidth', 1);

% Add p-value
y_pos = q3 + 0.2 * range(diff_bg); % Position for p-value
text(x_center - box_width, y_pos, sprintf('p = %.5f', bg_p), ...
'HorizontalAlignment', 'center', 'FontSize', 7, 'FontWeight', 'bold', 'Color', 'red');

% Beautify plot
xlim([0.5, 1.5]);
xticks(x_center);
xticklabels({'diff\_bg'});
set(gca, 'FontName', 'Arial', 'FontSize', 5, 'FontWeight', 'bold');
ylabel('Values', 'FontName', 'Arial', 'FontWeight', 'bold', 'FontSize', 7);
title('CFA:BasGang', 'FontName', 'Arial', 'FontWeight', 'bold', 'FontSize', 7);
grid on;
hold off;

% Save figure %
saveas(gcf, fullfile(filepath1, 'BasGang_stat'), 'epsc');
% Thalamus not needed, already done  %

%% All together with thalamus multiple comparisions plot %%


% === Step 1: Global z-score normalization ===
all_GFP = [TG; BG; MG; PoG; MeG; GG];
all_RFP = [TR; BR; MR; PoR; MeR; GR];
all_values = [all_GFP; all_RFP];

mean_val = mean(all_values);
std_val = std(all_values);

TG_z = (TG - mean_val) / std_val;
TR_z = (TR - mean_val) / std_val;
BG_z = (BG - mean_val) / std_val;
BR_z = (BR - mean_val) / std_val;
MG_z = (MG - mean_val) / std_val;
MR_z = (MR - mean_val) / std_val;
PoG_z = (PoG - mean_val) / std_val;
PoR_z = (PoR - mean_val) / std_val;
MeG_z = (MeG - mean_val) / std_val;
MeR_z = (MeR - mean_val) / std_val;
GG_z = (GG - mean_val) / std_val;
GR_z = (GR - mean_val) / std_val;

% === Step 2: Compute differences ===
T_diff = TR_z - TG_z;
B_diff = BR_z - BG_z;
M_diff = MR_z - MG_z;
Po_diff = PoR_z - PoG_z;
Me_diff = MeR_z - MeG_z;
G_diff = GR_z - GG_z;

region_names = {'Thalamus', 'Brainstem', 'Midbrain', 'Pons', 'Medulla', 'Basal Ganglia'};
diff_groups = {T_diff, B_diff, M_diff, Po_diff, Me_diff, G_diff};

% === Step 3: Test for normality ===
is_normal = true(1, length(diff_groups));
alpha_normal = 0.05;

fprintf('Normality test results (Lilliefors):\n');
for i = 1:length(diff_groups)
[h, p] = lillietest(diff_groups{i});
if h == 0
status = 'Normal';
else
status = 'Non-normal';
end
fprintf('%s: p = %.3f -> %s\n', region_names{i}, p, status);
is_normal(i) = (h == 0);
end

% === Step 4: Combine data and labels ===
all_diff = [];
group_labels = [];
for i = 1:length(diff_groups)
all_diff = [all_diff; diff_groups{i}(:)];
group_labels = [group_labels; repmat(region_names(i), length(diff_groups{i}), 1)];
end

th_idx = find(strcmp(region_names, 'Thalamus'));

% === Step 5: Run test based on normality ===
if all(is_normal)
% One-way ANOVA
fprintf('\nAll groups normal. Running One-way ANOVA...\n');
[p_anova, tbl, stats] = anova1(all_diff, group_labels, 'off');
fprintf('ANOVA p = %.4f\n', p_anova);

c = multcompare(stats, 'Display', 'off');

th_posthoc_p = nan(1, length(region_names));
for i = 1:size(c,1)
g1 = c(i,1); g2 = c(i,2); pval = c(i,6);
if g1 == th_idx
th_posthoc_p(g2) = pval;
elseif g2 == th_idx
th_posthoc_p(g1) = pval;
end
end
ks_p = nan(1, length(region_names));
method_used = 'One-way ANOVA + Tukey';
else
% Kruskal-Wallis + fallback to KS
fprintf('\nNon-normal group(s) detected. Running Kruskal-Wallis...\n');
[p_kw, ~, stats] = kruskalwallis(all_diff, group_labels, 'off');
fprintf('Kruskal-Wallis p = %.8f\n', p_kw);

th_data = diff_groups{th_idx};
th_posthoc_p = nan(1, length(region_names));
ks_p = nan(1, length(region_names));

fprintf('\nPost-hoc comparisons (vs Thalamus):\n');
for i = 1:length(diff_groups)
if i == th_idx, continue; end
[p_rs, ~] = ranksum(th_data, diff_groups{i});
th_posthoc_p(i) = p_rs;

if p_rs >= 0.05
[~, p_ks] = kstest2(th_data, diff_groups{i});
ks_p(i) = p_ks;
if p_ks < 0.05
sig_marker = '*';
else
sig_marker = '';
end
fprintf('%s: ranksum p = %.8f (NS) -> KS test p = %.3f %s\n', region_names{i}, p_rs, p_ks, sig_marker);
else
fprintf('%s: ranksum p = %.8f *\n', region_names{i}, p_rs);
end
end

method_used = 'Kruskal-Wallis + Ranksum (KS fallback)';
end

% === Step 6: Define colors for each group (corrected) ===
colors = [
0, 1, 1;      % Thalamus 
0, 0, 1; % Brainstem 
1.0000, 0.4118, 0.7059; % Midbrain
1, 0.5, 0; % Pons 
0.5804, 0.0000, 0.8275; % Medulla
1.0000, 0.0000, 0000; % Basal Ganglia 
];

% === Step 7: Manual box + median + whiskers + jittered scatter plot ===
figure('Color','w','Position',[100 100 900 500]);
hold on;

num_groups = length(region_names);
positions = 1:num_groups;
box_width = 0.5;  % Width of each box

for i = 1:num_groups
group_data = diff_groups{i};
x = positions(i);

% Calculate boxplot stats
q1 = prctile(group_data, 25);
med = median(group_data);
q3 = prctile(group_data, 75);
iqr_val = q3 - q1;
whisker_low = max(min(group_data), q1 - 1.5*iqr_val);
whisker_high = min(max(group_data), q3 + 1.5*iqr_val);

% Draw box
rectangle('Position', [x - box_width/2, q1, box_width, q3 - q1], ...
'EdgeColor', colors(i,:), 'LineWidth', 1, 'LineStyle', '-');

% Draw median line (color matched to box)
plot([x - box_width/2, x + box_width/2], [med, med], ...
'Color', colors(i,:), 'LineWidth', 1);

% Draw whiskers
plot([x, x], [whisker_low, q1], 'Color', colors(i,:), 'LineWidth', 1);
plot([x, x], [q3, whisker_high], 'Color', colors(i,:), 'LineWidth', 1);

% Draw whisker caps
cap_width = box_width * 0.25;
plot([x - cap_width/2, x + cap_width/2], [whisker_low, whisker_low], 'Color', colors(i,:), 'LineWidth', 1);
plot([x - cap_width/2, x + cap_width/2], [whisker_high, whisker_high], 'Color', colors(i,:), 'LineWidth', 1);

% Scatter individual points with jitter
jitterAmount = 0.15;
x_jitter = x + (rand(size(group_data)) - 0.5) * 2 * jitterAmount;
scatter(x_jitter, group_data, 80, 'MarkerEdgeColor', colors(i,:), ...
'MarkerFaceColor', colors(i,:), 'MarkerFaceAlpha', 0.6);
end

% Formatting
xlim([0.5, num_groups + 0.5]);
set(gca, 'XTick', positions, 'XTickLabel', region_names, 'FontSize', 12);
ylabel('Z-scored RFP - GFP Difference', 'FontSize', 14);
title(sprintf('%s across regions', method_used), 'FontSize', 16);

% Annotate significant p-values near the top of plot
yl = ylim;
text_y = yl(2) * 0.95;
th_idx = find(strcmp(region_names, 'Thalamus'));

for i = 1:length(region_names)
if i == th_idx, continue; end

if ~isnan(th_posthoc_p(i)) && th_posthoc_p(i) < 0.05
text(i, text_y, sprintf('p=%.7f', th_posthoc_p(i)), 'HorizontalAlignment', 'center', 'FontSize', 12);
text(i, text_y * 0.98, '*', 'FontSize', 18, 'FontWeight', 'bold', ...
'HorizontalAlignment', 'center', 'Color', 'r');
elseif exist('ks_p','var') && ~isnan(ks_p(i)) && ks_p(i) < 0.05
text(i, text_y, sprintf('KS p=%.7f', ks_p(i)), 'HorizontalAlignment', 'center', 'FontSize', 12);
text(i, text_y * 0.98, '*', 'FontSize', 18, 'FontWeight', 'bold', ...
'HorizontalAlignment', 'center', 'Color', 'm');
end
end

hold off;

% === Step 8: Save current figure as vector EPS (color) with painters renderer ===
set(gcf, 'Renderer', 'painters'); % ensure vector output
print(gcf, fullfile(filepath1, 'All_thalamus_compare'), '-depsc', '-painters');

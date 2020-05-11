% ========================= Project information ===========================
% Authors: Junjie Zhang, Kourosh Khoshelham
% Paper title: 3D reconstruction of internal wood decay using
% photogrammetry and sonic tomography
% =========================================================================
% ========================= Script information ============================
% This script reads the valid control points and summarise the semi-
% variograms later used for interpolation.
% =========================================================================

% Load data
load truncated_control_points.mat;

% Calculate the differences between the values of every two control points
% [index_n index_m distance delta_velocity]
differences = [];
for n = 1:size(control_points, 1)
    for m = 1:size(control_points, 1)
        if (n >= m)
            continue;
        end
        
        % Show progress
        display(['Computing for n = ', num2str(n), ' ', 'and m = ', num2str(m)])
        
        % First control point
        xn = control_points(n, 5);
        yn = control_points(n, 6);
        zn = control_points(n, 7);
        vn = control_points(n, 8);
        
        % Second control point
        xm = control_points(m, 5);
        ym = control_points(m, 6);
        zm = control_points(m, 7);
        vm = control_points(m, 8);
        
        % Distance between the two control points
        d = sqrt((xn-xm)^2 + (yn-ym)^2 + (zn-zm)^2);
        
        differences = [differences; n m d vn-vm];
    end
end

% Save the results
save('differences.mat', 'differences');

% Plot
figure
scatter(differences(:, 3), differences(:, 4).^2, 'r.')
title('Semivariogram cloud', 'FontSize', 20)
xlabel('Distance(m)', 'FontSize', 20)
ylabel('Squared velocity differences (m^{2}/s^{2})', 'FontSize', 20)

% Experimental Semivariogram
% Lag width
delta = 0.05;

% Semivariogram cloud shown as a boxplot
S_boxplot = [];
S_boxplot_g = [];

for n = 1:size(differences, 1)
    distance = differences(n, 3);
    difference = differences(n, 4);
    
    lag = round(distance / delta) * delta;
    
    S_boxplot = [S_boxplot; difference^2];
    S_boxplot_g = [S_boxplot_g; lag];
end

% Plot the squared velocity differences within lags
figure
boxplot(S_boxplot, S_boxplot_g)
title('Boxplot of the semivariances', 'FontSize', 20)
xlabel('Lag (m)', 'FontSize', 20)
ylabel('Squared velocity differences (m^{2}/s^{2})', 'FontSize', 20)

% Experimental semivariogram
% (gamma, Eq 10.13, P293, Geographic Information Analysis 2nd Edition)
semivariances = [];
for n = min(unique(S_boxplot_g)):delta:max(unique(S_boxplot_g))
    count = 0;
    sum = 0;
    
    for m = 1:length(S_boxplot_g)
        if S_boxplot_g(m, 1) == n
            count = count + 1;
            sum = sum + S_boxplot(m, 1);
        end
    end
    
    semivariances = [semivariances; n sum/count];
end

% Remove the decreasing semivariances at the end (sill & range)

% Fit a function to the semivariances. Use cftool and save as svfun.m

figure
scatter(semivariances(:, 1), semivariances(:, 2), 'r', 'DisplayName', 'Semivariances')
title('Semivariance function', 'FontSize', 20)
xlabel('Distance (m)', 'FontSize', 20)
ylabel('Semivariance / gamma (m^{2}/s^{2})', 'FontSize', 20)

% svfun.m contains the function I fit to this case
hold on
fplot(@(x) svfun(x), [semivariances(1, 1) 1.2], 'r', 'DisplayName', 'Semivariance function')
legend

% Save the semivariances for Kriging interpolation
save('semivariances.mat', 'semivariances');
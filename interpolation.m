% ========================= Project information ===========================
% Authors: Junjie Zhang, Kourosh Khoshelham
% Paper title: 3D reconstruction of internal wood decay using
% photogrammetry and sonic tomography
% =========================================================================
% ========================= Script information ============================
% This script reads the tree model and the valid control points, and uses
% the function in svfun.m which is fit to the semi-variogram to interpolate
% the internal wood health. It exports the point cloud with coordinates and
% wood health indicated by stress wave velocities in results.csv. It can be
% visualised in software such as CloudCompare.
% =========================================================================

% Load data
load trunk.mat;
% load measurements.mat;
load truncated_control_points.mat;
load Markers_transformed.mat;

% Boundaries for interpolation
% xmin = trunk.XLimits(1, 1);
% xmax = trunk.XLimits(1, 2);
% ymin = trunk.YLimits(1, 1);
% ymax = trunk.YLimits(1, 2);
% zmin = trunk.ZLimits(1, 1);
% zmax = trunk.ZLimits(1, 2);
xmin = min(Markers(:, 1));
xmax = max(Markers(:, 1));
ymin = min(Markers(:, 2));
ymax = max(Markers(:, 2));
zmin = min(Markers(:, 3));
zmax = max(Markers(:, 3));
inc = 0.01; % Increment; density of the interpolated points
[mx, my, mz] = meshgrid(xmin:inc:xmax, ymin:inc:ymax, zmin:inc:zmax);
mx = double(mx(:));
my = double(my(:));
mz = double(mz(:));

% Remove the irrelevant points, same algorithm used for control points
xt = trunk.Location(:, 1);
yt = trunk.Location(:, 2);
zt = trunk.Location(:, 3);
it = 1:100:length(xt);
xt_s = double(xt(it));
yt_s = double(yt(it));
zt_s = double(zt(it));
ri = sqrt(xt_s.^2 + yt_s.^2); % radius of tree points
ri1 = griddata(xt_s, yt_s, zt_s, ri, mx, my, mz);

ic = ~isnan(ri1); % inside the tree

% X, Y, Z of the points inside the tree trunk for interpolation
imx = mx(ic);
imy = my(ic);
imz = mz(ic);

% Ordinary Kriging
% Semivariance function can be found in svfun.m
% P304, Geographic Information Systems 2nd Edition
A = [];
for n = 1:size(control_points, 1)
    for m = 1:size(control_points, 1)
        d = sqrt((control_points(n, 5)-control_points(m, 5))^2 + (control_points(n, 6)-control_points(m, 6))^2 + (control_points(n, 7)-control_points(m, 7))^2);
        A(n, m) = svfun(d);
    end
end

A = [A ones(size(A, 1), 1)];
A = [A; ones(1, size(A, 2))];
A(end, end) = 0;

ipv = [];
for n = 1:size(imx, 1)
    px = imx(n, 1);
    py = imy(n, 1);
    pz = imz(n, 1);
    
    b = [];
    for m = 1:size(control_points, 1)
        cx = control_points(m, 5);
        cy = control_points(m, 6);
        cz = control_points(m, 7);
        
        d = sqrt((px-cx)^2 + (py-cy)^2 + (pz-cz)^2);
        
        b = [b; svfun(d)];
    end
    
    b = [b; 1];
    
    w = A \ b;
    
    % Remove negative weights and standardise the rest
    weights = w(1:end-1);
    positive_weights_sum = sum(weights(weights>0));
    for i = 1:length(w)-1
        if w(i) < 0
            w(i) = 0;
        else
            w(i) = w(i) / positive_weights_sum;
        end
    end

    v = 0;
    for q = 1:size(control_points, 1)
        v = v + w(q) * control_points(q, 8);
    end
    
    ipv = [ipv; v];
    
    display(['Progress: ', num2str(n/size(imx, 1)*100), '%'])
end

% Interpolated points
% [X Y Z Velocity/Density or health of wood]
ip = [imx imy imz ipv];

% Export the data for CloudCompare
csvwrite('results.csv', ip)

% Plot the results with the trunk
pcshow(trunk)
hold on
scatter3(ip(:, 1), ip(:, 2), ip(:, 3), 14, ip(:, 4), 'filled')

% ----------
% Execute the following separately
% ----------

% Quality check
load check_points.mat;

cpvc = []; % [Computed_velocity Interpolated_velocity RMS]
for n = 1:size(check_points, 1)
    cpx = check_points(n, 5);
    cpy = check_points(n, 6);
    cpz = check_points(n, 7);
    b = [];
    for m = 1:size(control_points, 1)
        cx = control_points(m, 5);
        cy = control_points(m, 6);
        cz = control_points(m, 7);
        
        d = sqrt((cpx-cx)^2 + (cpy-cy)^2 + (cpz-cz)^2);
        
        b = [b; svfun(d)];
    end
    
    b = [b; 1];
    
    w = A \ b;
    
    % Remove negative weights and standardise the rest
    weights = w(1:end-1);
    positive_weights_sum = sum(weights(weights>0));
    for i = 1:length(w)-1
        if w(i) < 0
            w(i) = 0;
        else
            w(i) = w(i) / positive_weights_sum;
        end
    end

    v = 0;
    for q = 1:size(control_points, 1)
        v = v + w(q) * control_points(q, 8);
    end
    
    cpvc = [cpvc; v];
    
    display(['Progress: ', num2str(n/size(check_points, 1)*100), '%'])
end

cpvc = [cpvc check_points(:, 8)];
cpvc(:, 3) = (cpvc(:, 1) - cpvc(:, 2)).^2;
rmse = sqrt(mean(cpvc(:, 3)));
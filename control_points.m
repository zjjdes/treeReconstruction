% ========================= Project information ===========================
% Authors: Junjie Zhang, Kourosh Khoshelham
% Paper title: 3D reconstruction of internal wood decay using
% photogrammetry and sonic tomography
% =========================================================================
% ========================= Script information ============================
% This script reads all the input data (TOF, sensor coordinates and 3D 
% model of the tree trunk) and computes the information of all the
% candidate control points for the interpolation.
% =========================================================================

format long g;

% Cordinates of the markers handpicked in CloudCompare
% [X Y Z]
load Markers.mat;

% Measured Time-of-Flight obtained from Arbotom
% Runtimes(n, m) --> TOF between sensor #n and sensor #m in us
load Runtimes.mat;

% Points of the point cloud of the surface of the trunk
% [X Y Z R G B nX nY nZ]
ptCloud = pcread("trunk.ply");

% Empirically determined settings
dist_lim = 0.01; % Threshold of the distance between two rays to generatea control point
weight_hori = 0.75; % The weight of the horizontal ray when computing the velocity of the control point

% Rotate the point cloud around x axis by 75 degrees
a = 75;
A = [1 0 0 0; 0 cosd(a) sind(a) 0; 0 -sind(a) cosd(a) 0; 0 0 0 1];
ptCloud_r = pctransform(ptCloud, affine3d(A));

% Rotate the markers in the same way
Markers = transformPointsForward(affine3d(A), Markers);

% Use the mean (x, y) as the origin
xt = ptCloud_r.Location(:, 1);
yt = ptCloud_r.Location(:, 2);
x0 = mean(xt);
y0 = mean(yt);
A = [1 0 0 0; 0 1 0 0; 0 0 1 0; -x0 -y0 0 1];

% Translate the point cloud and the markers
trunk = pctransform(ptCloud_r, affine3d(A));
Markers = transformPointsForward(affine3d(A), Markers);

% Save the transformed point cloud
save('trunk.mat', 'trunk');
pcwrite(trunk,'rotated_trunk.ply','Encoding','ascii');

% Save the transformed markers
save('Markers_transformed.mat', 'Markers');

% Plot the tree trunk point cloud
% Requires Computer Vision Toolbox
% pcshow(trunk)

% Coordinates of the measurement rays
% [start_# end_# x_start x_end y_start y_end z_start z_end TOF Length Velocity]
measurements = [];
for n = 1:size(Markers, 1)
    for m = 1:size(Markers, 1)
        if n >= m
            continue
        end
        
        x_start = Markers(n, 1);
        y_start = Markers(n, 2);
        z_start = Markers(n, 3);
        x_end = Markers(m, 1);
        y_end = Markers(m, 2);
        z_end = Markers(m, 3);
        
        % Due to random errors, measurements between the same two
        % sensors might have different measured runtimes. Taking the
        % average of the runtimes as the input.
        tof = mean(Runtimes(n, m), Runtimes(m, n));
        
        measurements = [measurements; n m x_start x_end y_start y_end ...
            z_start z_end tof];
    end
end

% Plot the measurement rays
% hold on
% plot3(measurements(:, 3:4)', measurements(:, 5:6)', measurements(:, 7:8)')

% Compute the length of each measurement
measurements = [measurements sqrt((measurements(:, 4) -...
    measurements(:, 3)).^2 + (measurements(:, 6) -...
    measurements(:, 5)).^2 + (measurements(:, 8) -...
    measurements(:, 7)).^2)];

% Compute the velocities of each measurement
measurements = [measurements measurements(:, 10) ./ (measurements(:, 9)...
    ./ 10^6)];

% Save the measurements
save('measurements.mat', 'measurements');

% Compute the coordinates of the control points
% [Start_sensor_first_measurement End_sensor_first_measurement Start_sensor_second_measurement End_sensor_second_measurement X Y Z Velocity]
all_control_points = [];
for n = 1:size(measurements, 1)
    for m = 1:size(measurements, 1)
        % No sensor is shared between the two measurements
        if n >= m || measurements(n, 1) == measurements(m, 1) ||...
                measurements(n, 2) == measurements(m, 2) ||...
                measurements(n, 1) == measurements(m, 2) ||...
                measurements(n, 2) == measurements(m, 1)
            continue
        end
        
        % If the two measurements are vertical, skip
        % Only needed when the sensors are places in a grid shape
        % In this specific case, there are 8 sensors on each level
        if rem(measurements(m, 1), 8) == rem(measurements(m, 2), 8)
            if rem(measurements(m, 1), 8) ==...
                    rem(measurements(n, 1), 8) ||...
                    rem(measurements(m, 1), 8) ==...
                    rem(measurements(n, 2), 8)
                continue
            end
        end
        
        if rem(measurements(n, 1), 8) == rem(measurements(n, 2), 8)
            if rem(measurements(n, 1), 8) ==...
                    rem(measurements(m, 1), 8) ||...
                    rem(measurements(n, 1), 8) ==...
                    rem(measurements(m, 2), 8)
                continue
            end
        end
        
        [xi, yi, zi] = deal(measurements(n, 3), measurements(n, 5),...
            measurements(n, 7));
        [xj, yj, zj] = deal(measurements(n, 4), measurements(n, 6),...
            measurements(n, 8));
        [xc, yc, zc] = deal(measurements(m, 3), measurements(m, 5),...
            measurements(m, 7));
        [xd, yd, zd] = deal(measurements(m, 4), measurements(m, 6),...
            measurements(m, 8));
        
        % Direction vectors of Line ij and Line cd
        vij = [xj - xi yj - yi zj - zi];
        vcd = [xd - xc yd - yc zd - zc];
        
        % Direction vector of Line l' (perpendicular to Lines ij and cd)
        v3 = cross(vij, vcd);
        
        % Solve for t1, t2 and t3
        % Equation 12: (xi, yi, zi)+t1*vij+t3*v3=(xc, yc, zc)+t2*vcd
        syms st1 st2 st3
        S = solve(xi+st1*vij(1, 1)+st3*v3(1, 1)==xc+st2*vcd(1, 1),...
            yi+st1*vij(1, 2)+st3*v3(1, 2)==yc+st2*vcd(1, 2),...
            zi+st1*vij(1, 3)+st3*v3(1, 3)==zc+st2*vcd(1, 3));

        
        % Convert syms of t1, t2, t3 to numeric double precision
        t1 = double(S.st1);
        t2 = double(S.st2);
        t3 = double(S.st3);
        
        % Compute Pij
        Pij = [xi yi zi] + t1 * vij;
        
        % Compute Pcd
        Pcd = [xc yc zc] + t2 * vcd;
        
        % Compute P (the mid point between Pij and Pcd)
        P = (Pij + Pcd) / 2;
        
        % Add to the list of valid control points if the two rays are
        % reasonably close. The velocity value of the control point is the
        % average of the two rays.
        if sqrt((Pcd(1, 1) - Pij(1, 1))^2 + (Pcd(1, 2) - Pij(1, 2))^2 +...
                (Pcd(1, 3) - Pij(1, 3))^2) < dist_lim
            
            % At least one of the measurements lies on a horizontal plane
            if floor(measurements(n, 1)/8) == floor(measurements(n, 2)/8) || floor(measurements(m, 1)/8) == floor(measurements(m, 2)/8)
                % The measurement on the horizontal planes weighs more than
                % the one which is not
                num_sen = 8; % Number of sensors on each level
                if floor(measurements(n, 1)/num_sen) == floor(measurements(n, 2)/num_sen) && floor(measurements(m, 1)/num_sen) == floor(measurements(m, 2)/num_sen)
                    V = mean([measurements(n, 11) measurements(m, 11)]);
                elseif floor(measurements(n, 1)/num_sen) == floor(measurements(n, 2)/num_sen)
                    V = measurements(n, 11) * weight_hori + measurements(m, 11) * (1 - weight_hori);
                elseif floor(measurements(m, 1)/num_sen) == floor(measurements(m, 2)/num_sen)
                    V = measurements(n, 11) * (1 - weight_hori) + measurements(m, 11) * weight_hori;
                end
                
                all_control_points = [all_control_points; measurements(n, 1) measurements(n, 2) measurements(m, 1) measurements(m, 2) P V];
            end
        end
    end
end

% hold on
% plot3(all_control_points(:, 1), all_control_points(:, 2), all_control_points(:, 3), 'og')

% Save all the control points since the computation takes a long time
save('all_control_points.mat', 'all_control_points');
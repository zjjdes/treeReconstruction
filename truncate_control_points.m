% ========================= Project information ===========================
% Authors: Junjie Zhang, Kourosh Khoshelham
% Paper title: 3D reconstruction of internal wood decay using
% photogrammetry and sonic tomography
% =========================================================================
% ========================= Script information ============================
% This script reads the rotated model of the tree trunk and the generated
% control points, and removes those outside of the tree trunk.
% =========================================================================

% Load data
load trunk.mat;
load all_control_points.mat;

% https://www.mathworks.com/matlabcentral/answers/474394-how-to-test-if-points-are-inside-a-point-cloud-model?s_tid=prof_contriblnk
% Remove the control points outside of the tree trunk
% generate new radius for control points
% i imagine it like a kind of projection on the tree
xt = trunk.Location(:, 1);
yt = trunk.Location(:, 2);
zt = trunk.Location(:, 3);
% xc = double(all_control_points(:, 1));
% yc = double(all_control_points(:, 2));
% zc = double(all_control_points(:, 3));
all_control_points = double(all_control_points);
it = 1:100:length(xt); % reduce data for faster calculation

% Keep only the data inside the trunk
xt_s = double(xt(it));
yt_s = double(yt(it));
zt_s = double(zt(it));
ri = sqrt(xt_s.^2 + yt_s.^2); % radius of tree points
ri1 = griddata(xt_s, yt_s, zt_s, ri, all_control_points(:, 5), all_control_points(:, 6), all_control_points(:, 7));

ic = ~isnan(ri1);  % inside the tree

fstart = all_control_points(:, 1);
fend = all_control_points(:, 2);
sstart = all_control_points(:, 3);
send = all_control_points(:, 4);
xc = all_control_points(:, 5);
yc = all_control_points(:, 6);
zc = all_control_points(:, 7);
vc = all_control_points(:, 8);

% [Start_sensor_first_measurement End_sensor_first_measurement Start_sensor_second_measurement End_sensor_second_measurement X Y Z Velocity]
control_points = [fstart(ic) fend(ic) sstart(ic) send(ic) xc(ic) yc(ic) zc(ic) vc(ic)];

% Remove 50 random points for validation
check_points = [];
for i=1:50
    j = randi([1 size(control_points, 1)]);
    check_points = [check_points; control_points(j, :)];
    control_points(j, :) = [];
end

% Save the selected control points
save('truncated_control_points.mat', 'control_points');
save('check_points.mat', 'check_points');

pcshow(trunk)
hold on
scatter3(control_points(:, 5), control_points(:, 6), control_points(:, 7), 14, control_points(:, 8), 'filled')
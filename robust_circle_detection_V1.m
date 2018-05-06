% Detects circle in laser scan data
% Searches full dataset to find circle
clear
load playpensample.mat
d = 0.1; % distance between lidar and neado origin
rot_angle = 0;
translate_x = 0;
translate_y = 0;

%remove zero radius data
theta_clean1 = [];
r_clean1 = [];

for i = 1:size(r)
    if r(i) > 0
        r_clean1 = [r_clean1;r(i)];
        theta_clean1 = [theta_clean1;theta(i)];
    end
end

%convert to cartesian
scan1_points = [];
for i = 1:size(r_clean1)
    theta = degtorad(theta_clean1(i));
    r = r_clean1(i);
    [x1,y1] = pol2cart(theta,r);
    scan1_points = [scan1_points, [x1;y1;1]];
end

%translate data to Neato coordinate frame
scan1_points_neato = translate(0,d)*scan1_points;

%translate Neato to room coordinate frame
scan1_points_room = translate(translate_x,translate_y)*rotate(rot_angle)*scan1_points_neato;
scan1_points = scan1_points_room(1:2,:)';


%% RANSAC FOR CIRCLES
r = 0.22; %circle radius
% 
% %TEST WITH CIRCLE DATA
% x = [0.3*cos(0.2) + 2; 0.3*cos(1) + 2; 0.3*cos(1.5) + 2; 0.3*cos(pi) + 2; 0.3*cos(2) + 2];
% y = [0.3*sin(0.2) + 2; 0.3*sin(1) + 2; 0.3*sin(1.5) + 2; 0.3*sin(pi) + 2; 0.3*sin(2) + 2];
% data1 = horzcat(x,y);
% data = data1;

%TEST WITH SCAN DATA
data = scan1_points;

% global variables to tune
num_times = 1000;
threshold = 0.01;
bound = 0.5; %select points within bound

%store result
result = struct();
result.num_inlier = 0;
result.c = [];

for k = 1:num_times

    %select a pair
    %randomly select one point
    point1 = datasample(data,1);
    % point1 = [1.2913026,0.25855201];
    xmin = point1(1) - bound;
    xmax = point1(1) + bound;
    ymin = point1(2) - bound;
    ymax = point1(2) + bound;
    [xrow,xcol] = find(data(:,1)<xmax & data(:,1)>xmin);
    [yrow,ycol] = find(data(:,2)<ymax & data(:,2)>ymin);

    row = intersect(xrow, yrow);

    neighbours = [];
    for i = 1:size(row)
        neighbour = data(row(i),:);
        neighbours = [neighbours; neighbour];
    end
    
    % avoid fitting line for noise
     if size(neighbours,1) > 1
        %select the second point
        point2 = datasample(neighbours,1);

        %compute c
        %find direction
        l = [point2(1) - point1(1), point2(2) - point1(2)];
        l_mag = norm(l);
        l_unit = 1./l_mag*l;
        n1 = [-l(2),l(1)];
        n2 = [l(2),-l(1)];
        midpoint = [(point1(1) + point2(1))/2, (point1(2) + point2(2))/2];

        %compute dot product to determine n pointing away from origin
        dotproduct1 = dot(midpoint,n1);

        if dotproduct1 > 0
            c_dir = n1;
        else
            c_dir = n2;
        end

        %compute c magnitude
        c_mag = sqrt(r.^2 - (l_mag/2).^2);
        c = c_mag*c_dir + midpoint;

        %define inlier or outlier
        inliers = [];
        outliers = [];
        for i = 1:size(data)
            P = data(i,:);
            %compute distance
            d = abs(norm(P-c) - r);
            %classify inlier and outlier
            if d <= threshold
                inliers = [inliers;P];
            else
                outliers = [outliers;P];
            end   
        end

        %calculate number of inliers
        size_inlier = size(inliers);
        num_inlier = size_inlier(1);

        %update result
        if num_inlier >= result.num_inlier
            result.num_inlier = num_inlier;
            result.inliers = inliers;
            result.c = c;
            result.midpoint = midpoint;
        end
     end
end


%define a circle
x_predict = [];
y_predict = [];
theta = 0:0.1:2*pi;
for angle = 1:size(theta,2)
    x = r*cos(theta(1,angle)) + result.c(1);
    y = r*sin(theta(1,angle)) + result.c(2);
    x_predict = [x_predict; x];
    y_predict = [y_predict; y];
end
figure
hold on
plot(x_predict(:,1), y_predict(:,1),'g-','DisplayName', 'RANSAC','LineWidth',3)
plot(data(:,1), data(:,2),'bx','DisplayName', 'Data')
plot(result.inliers(:,1),result.inliers(:,2),'ro','DisplayName','inliers')
%plot(result.midpoint(:,1),result.midpoint(:,2),'yo','DisplayName','midpoint')
%plot(result.c(:,1),result.c(:,2),'mo','DisplayName','center')

xlabel('x')
ylabel('y')
axis('equal')
title('Circle Fitting for PlayPenSample')
legend('show')
axis('equal')

function res = rotate(theta)
    res = [cos(theta), sin(theta), 0; -sin(theta), cos(theta), 0; 0, 0, 1];
end

function res = translate(a,b)
    res = [1 0 a;0 1 b;0 0 1];
end
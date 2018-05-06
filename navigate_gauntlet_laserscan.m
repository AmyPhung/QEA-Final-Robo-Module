% Generates wheel velocity commands for a given set of data

clear
load line_detection_scan1.mat

% Creating meshgrid for gauntlet
x = -2.413:0.1:.6096;
y = -.3048:0.1:2.54;
[X,Y] = meshgrid(x,y);
all_line_eq = {};
cellindex = 1;

for i = 1:size(lines,1)
    for u = 1:size(lines(i).segments, 2)
        points = lines(i).segments{u};
        
        %compute length of segment
        total_num = size(points,1);
        start_point = points(1,:);
        x1 = start_point(1);
        y1 = start_point(2);
        end_point = points(total_num,:);
        x2 = end_point(1);
        y2 = end_point(2);
        length = norm([x2-x1, y2-y1]);
        
        if length > 0.2
            slope = lines(i).m;
            line_eq = gradient_field(x1,y1,x2,slope,X,Y);
            all_line_eq{cellindex,1} = line_eq;
            cellindex = cellindex + 1;
        end
    end
end



% Create lines for circle
cx = -.3302;
cy = 1.8034;
% cx = -0.33;
% cy = 1.80;
radius = 0.27;
fun = @(theta)log((X-(radius*cos(theta)+cx)).^2+(Y-(radius*sin(theta)+cy)).^2);
c = -integral(fun,0,2*pi,'ArrayValued',true);

Z = zeros();
for line = 1:numel(all_line_eq)
    Z = Z + all_line_eq{line};
end
Z = Z + 0.7*c;
%40*b2l1+b2l2
%3*b1l1+1.7*b1l2
%30*b2l1+b2l2
%30*b2l1+b2l2+2*b1l1+1.5*b1l2)*20
%Z = (3*b1l1+1.7*b1l2+40*b2l1+b2l2)+.7*c;
% Z = (b1l1+b1l2+b2l1+b2l2)+.7*c;
figure;
contour(X,Y,Z,'ShowText','on')
hold on;
xlim([-2.413,.6096]);
ylim([.3048,2.54]);


[FX,FY] = gradient(Z);

%% CODE START %%%
%Requires: [FX, FY], Z (results from integral)

TIME_STEP = 2; % in seconds
INITIAL_LOCATION = [-0.62,0]; % in meters from the origin
INITIAL_HEADING = pi/2; % in radians from the x-axis
%INITIAL_HEADING = 0; % in radians from the x-axis
MAX_SPEED = 0.1;
STOP_THRESHOLD = 10; % threshold when NEATO will stop
MAX_RUN_TIME = 10; % maximum seconds to drive NEATO
d = 0.24; % NEATO wheelbase size

x0 = INITIAL_LOCATION(1);
y0 = INITIAL_LOCATION(2);
max_distance = MAX_SPEED*TIME_STEP;

%scan data
%need v1, v2, v3.... (results from contour, usually named q)

current_theta = INITIAL_HEADING;
current_x = x0;
current_y = y0;
commands = zeros(MAX_RUN_TIME/TIME_STEP,2);
step = 1;

%for debugging
points = zeros(MAX_RUN_TIME/TIME_STEP,2);


for test = 1:40

    points(test,:) = [current_x,current_y];

    %Compute direction of steepest descent  
    [xindex,yindex] = nearestIndex(current_x,current_y,x,y);

    fx = FX(yindex,xindex);
    fy = FY(yindex,xindex);
    g = sqrt(fx^2 + fy^2);

    z_value = Z(yindex,xindex);
    if z_value > STOP_THRESHOLD
        break
    end

    if fy >= 0
        desired_theta = acos(dot([1,0],[fx,fy])/norm([fx,fy]));
    elseif fy < 0
        desired_theta = 2*pi-acos(dot([1,0],[fx,fy])/norm([fx,fy]));
    end

    dtheta = desired_theta-current_theta;
    %Compute angular velocity for time step
    omega = dtheta/TIME_STEP; % in radians/s

    %Compute wheel commands
    V_L = double(MAX_SPEED - omega*d/2);
    V_R = double(MAX_SPEED + omega*d/2);
    commands(step,1) = V_L;
    commands(step,2) = V_R;

    %Find next point
    [dx,dy] = pol2cart(desired_theta,max_distance);

    %Update variables
    current_x = current_x + dx;
    current_y = current_y + dy;
    current_theta = current_theta + dtheta;
    step = step + 1;

end

plot(points(:,1),points(:,2),'o-');
axis equal

hold off;

function [xindex,yindex] = nearestIndex(px,py,xs,ys)
    %{
    ax = approx. of x
    yx = approx. of y
    px = x value to be approximated
    py = y value to be approximated
    xs = matrix of x values
    ys = matrix of y values
    %}
    xindex = dsearchn(xs',px);
    yindex = dsearchn(ys',py);
end

function res = gradient_field(x1,y1,x2,slope,X,Y)
    y0 = @(x0)slope*(x0-x1) + y1; 
    fun = @(x0)log((X-x0).^2+(Y-y0(x0)).^2);
    res = integral(fun,x1,x2,'ArrayValued',true);
end
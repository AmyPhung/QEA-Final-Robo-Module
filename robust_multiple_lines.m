%% Detect multiple lines in a dataset

%clear
%load pre collected data
load laserscan.mat

%collect data
% sub = rossubscriber('/stable_scan');
% scan_message = receive(sub);
% r_1 = scan_message.Ranges(1:end-1);
% theta_1 = [0:359]';

d = 0.1;% distance between lidar and neado origin
rot_angle = 0;
translate_x = 0;
translate_y = 0;

%hard origin coordinate(0.62,0)

%remove zero radius data
theta_clean1 = [];
r_clean1 = [];

for i = 1:size(r_1)
    if r_1(i) > 0
        r_clean1 = [r_clean1;r_1(i)];
        theta_clean1 = [theta_clean1;theta_1(i)];
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

%% line detection
data = scan1_points;

% global variables to tune
num_times = 1000;
threshold = 0.01;
bound = 0.3; %select neighbour points within bound
min_neighbours = 3;
gap = 0.1; %gap between line segment

%store result and all detected lines
lines = [];
num_outliers = size(data,1);
line_cell_index = 1;

while 0.04*size(scan1_points,1) < num_outliers
    result = struct();
    result.num_inliers = 1;
    result.num_outliers = 0;

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
        if size(neighbours,1) > min_neighbours
            %select the second point
            point2 = datasample(neighbours,1);

            %define a line
            m = (point2(2) - point1(2))./(point2(1) - point1(1));
            b = point1(2) - m.*point1(1);
            
            %define inlier or outlier
            inliers = [];
            outliers = [];
            
            for i = 1:size(data)
                P = data(i,:);
                A = [P(1) - point1(1), P(2) - point1(2)];
                v = [point2(1) - point1(1), point2(2) - point1(2)];
                %find unit vector on the line
                u = 1./norm(v).*v;
                %calculate projection
                projection = (A*u')*u;
                %calculate distance
                difference = A - projection;
                d = norm(difference);
                
                %classify inlier and outlier
                if d <= threshold
                    inliers = [inliers;P];
                else
                    outliers = [outliers;P];
                end 

                
            end
            
            %calculate number of inliers and outliers
            num_inlier = size(inliers,1);
            num_outlier = size(outliers,1);

            %update result
            if num_inlier >= result.num_inliers
                result.num_inliers = num_inlier;
                result.num_outliers = num_outlier;
                result.inliers = inliers;
                result.outliers = outliers;
                result.m = m;
                result.b = b;
                
                %categorize into line segments
                line_seg = inliers(1,:); %initial first segment
                cellindex=1;
                allsegments = {}; %initiate cell array
                
                for i = 1:size(inliers) - 1
                    distance = norm(inliers(i,:)- inliers(i+1,:));
                    
                    %check which line segement belongs
                    if distance < gap
                      line_seg = [line_seg;inliers(i+1,:)];
                    else
                        allsegments{cellindex}=line_seg;
                        cellindex = cellindex+1;
                        line_seg = [inliers(i+1,:)]; %create new segment
                    end
                    allsegments{cellindex} = line_seg; %save final result
                end
                result.segments = allsegments;
            end
            
        end
    end
    if numel(fieldnames(result)) > 2
       lines = [lines;result];
       num_outliers = result.num_outliers
       data = result.outliers;
    else
       num_outliers = 0;
    end
    
end

%plot result
figure
hold on
plot(scan1_points(:,1), scan1_points(:,2),'bx')
for i = 1:size(lines,1)
    for u = 1:size(lines(i).segments, 2)
        points = lines(i).segments{u};
        plot(points(:,1),points(:,2),'r-','LineWidth',2)
    end
end
hold off
xlabel('x')
ylabel('y')
title('Multiple Line Detection for Scan1')
legend('Data','Line Fitting')


function res = rotate(theta)
    res = [cos(theta), sin(theta), 0; -sin(theta), cos(theta), 0; 0, 0, 1];
end

function res = translate(a,b)
    res = [1 0 a;0 1 b;0 0 1];
end

% Attempts to fit a circle to each group of points classified as a line
% segment
% Require Running robust_multiple_lines.m beforehand


% global variables to tune
num_times = 1000;
threshold = 0.01;
r = 0.2; %circle radius

%store result
circle = struct();
circle.num_inlier = 0;
circle.c = [];


for each_line = 1:size(lines,1)
    for each_segment = 1:size(lines(each_line).segments, 2)
        segment_points = lines(each_line).segments{each_segment};
        
        %fit circle
        result = struct();
        result.num_inlier = 0;
        
        if size(segment_points,1) > 2
            for k = 1:num_times
           
                %select a pair
                pair = datasample(segment_points,2);
                point1 = pair(1,:);
                point2 = pair(2,:);

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
                for i = 1:size(segment_points)
                    P = segment_points(i,:);
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
                num_inlier = size(inliers,1);
                
                %update result
                if num_inlier >= result.num_inlier
                    result.num_inlier = num_inlier;
                    result.inliers = inliers;
                    result.c = c;
                    result.midpoint = midpoint;
                end
            
            end
        
            %update result
            if result.num_inlier >= circle.num_inlier
                circle.num_inlier = result.num_inlier;
                circle.inliers = result.inliers;
                circle.c = result.c;
                circle.midpoint = result.midpoint;
            end
        end
        
        
    end
end

%define a circle
x_predict = [];
y_predict = [];
theta = 0:0.1:2*pi;
for angle = 1:size(theta,2)
    x = r*cos(theta(1,angle)) + circle.c(1);
    y = r*sin(theta(1,angle)) + circle.c(2);
    x_predict = [x_predict; x];
    y_predict = [y_predict; y];
end
figure
hold on
plot(x_predict(:,1), y_predict(:,1),'g-','DisplayName', 'RANSAC','LineWidth',3)
plot(scan1_points(:,1), scan1_points(:,2),'bx','DisplayName', 'Data')
plot(circle.inliers(:,1),circle.inliers(:,2),'ro','DisplayName','inliers')
%plot(circle.midpoint(:,1),circle.midpoint(:,2),'yo','DisplayName','midpoint')
%plot(circle.c(:,1),circle.c(:,2),'mo','DisplayName','center')
hold off

xlabel('x')
ylabel('y')
axis('equal')
title('scan4 in cartesian')
legend('show')
axis('equal')

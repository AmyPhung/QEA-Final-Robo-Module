%select a pair
pair = [1,3;5,1];
point1 = pair(1,:);
point2 = pair(2,:);
r = 5;

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

c = c_mag*c_dir;

quiver(0,0,c(1),c(2))
hold on
plot(point1(1),point1(2),'rx')
plot(point2(1),point2(2),'bx')
hold off
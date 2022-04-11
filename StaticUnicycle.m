%Plot static unicycle

%variables
x = 0;
y = 0;
z = 0; %assumed to be constant in model, measured from bottom of wheel to ground
theta = 0; % steering angle in x-y plane in counter clockwise direction 0 at x-axis (to the left)
delta = 0; % sideways falling angle of cycle to the right, 0 when upright and rotated around contact point with ground
alpha = 0; % forward angle of seat with respect to wheel, rotating around the center of the wheel, 0 when upright

%Unicycle parameters
r = 1; %raduis of wheel
h= 2.5; %height of saddle from center of wheel

numwheel = 100; %number of points to plot wheel
centerwheel=zeros(3,1);
centerwheel(3,1)=r;
saddle = zeros(3,1);
saddle(3,1)=h;
angles = linspace(0,2*pi,numwheel);
wheelcoord = zeros(3,numwheel);
wheelcoord(1,:) = r*cos(angles);
wheelcoord(3,:)= r+r*sin(angles);
%rotations
Rotationmatrixdelta = Rotationmatrix(delta,'x');
Rotationmatrixtheta = Rotationmatrix(theta,'z');
Rotationmatrixalpha = Rotationmatrix(alpha,'y');
Totalrotationmatrix = Rotationmatrixtheta*Rotationmatrixdelta;
for i=1:numwheel
    wheelcoord(:,i)=Totalrotationmatrix*wheelcoord(:,i);
end
centerwheel = Totalrotationmatrix*centerwheel;
saddle = Rotationmatrixalpha*saddle;
saddle(3,1)=saddle(3,1)+r;
saddle=Totalrotationmatrix*saddle;

%translations
wheelcoord(1,:) = wheelcoord(1,:) +x*ones(1,numwheel);
wheelcoord(2,:) = wheelcoord(2,:) +y*ones(1,numwheel);
wheelcoord(3,:) = wheelcoord(3,:) +z*ones(1,numwheel);
centerwheel(1,1)=centerwheel(1,1)+x;
centerwheel(2,1)=centerwheel(2,1)+y;
centerwheel(3,1)=centerwheel(3,1)+z;
saddle(1,1)=saddle(1,1)+x;
saddle(2,1)=saddle(2,1)+y;
saddle(3,1)=saddle(3,1)+z;

saddleplot=[centerwheel,saddle];
%plotting
plot3(wheelcoord(1,:),wheelcoord(2,:),wheelcoord(3,:))
hold on
plot3(centerwheel(1),centerwheel(2),centerwheel(3),'.')
plot3(saddleplot(1,:),saddleplot(2,:),saddleplot(3,:),'-o')
xlim([-5,5])
ylim([-5,5])
zlim([0,10])
xlabel('x')
ylabel('y')
zlabel('z')


function Rotmat = Rotationmatrix(angle, axis)
    Rotmat = zeros(3,3);
    if strcmp(axis,'x')
        Rotmat(1,1)=1;
        Rotmat(2,2)=cos(angle);
        Rotmat(2,3)=-sin(angle);
        Rotmat(3,2)=sin(angle);
        Rotmat(3,3)=cos(angle);
    elseif strcmp(axis,'z')
        Rotmat(3,3)=1;
        Rotmat(2,2)=cos(angle);
        Rotmat(1,2)=-sin(angle);
        Rotmat(2,1)=sin(angle);
        Rotmat(1,1)=cos(angle);
    elseif strcmp(axis,'y')
        Rotmat(1,1)=cos(angle);
        Rotmat(2,2)=1;
        Rotmat(3,1)=-sin(angle);
        Rotmat(1,3)=sin(angle);
        Rotmat(3,3)=cos(angle);
    end
end


%Plot static unicycle
function vec =Unicyclemoviemaker(x,y,theta,epsilon,alpha,r,h,xlocs)
    %variables
    % xlocs = 'CM'if x and y contain the locations of the center of mass,
    % xlocs = 'bike'if x and y contain the locations of the bottom of the
    % wheel.
    %x center of mass, center of mass is plotted as saddle
    %y center of mass, center of mass is plotted as saddle
    %theta, steering angle in x-y plane in counter clockwise direction 0 at x-axis (to the left)
    %epsilon sideways falling angle of cycle to the right, 0 when upright and rotated around contact point with ground, take minus sign for answers from lagrangians those are defined with clockwise angles 
    %alpha forward angle of seat with respect to wheel, rotating around the center of the wheel, 0 when upright
    %r radius of wheel
    %h distance from center of wheel to center of mass, aka saddle
    
    rvec = zeros(3,1);
    rvec(3,:)=rvec(3,:)+r;
    numtime=length(x);
    numwheel = 100; %number of points to plot wheel
    timeddata = zeros(numtime,numwheel+2,3);
    angles = linspace(0,2*pi,numwheel);
    if strcmp(xlocs,'CM')
        for j=1:numtime
            centerwheel=zeros(3,1);
            centerwheel(3,1)=r;
            saddleinit = zeros(3,1);
            saddleinit(3,1)=h;
            saddle = zeros(3,1);
            saddle(1,1)=x(j);
            saddle(2,1)=y(j);
            wheelcoord = zeros(3,numwheel);
            wheelcoord(1,:) = r*cos(angles);
            wheelcoord(3,:)= r+r*sin(angles);
            %rotations
            Rotationmatrixepsilon = Rotationmatrix(epsilon(j),'x');
            Rotationmatrixtheta = Rotationmatrix(theta(j),'z');
            Rotationmatrixalpha = Rotationmatrix(alpha(j),'y');
            Totalrotationmatrix = Rotationmatrixtheta*Rotationmatrixepsilon;
            xback = saddle - Totalrotationmatrix*(Rotationmatrixalpha*saddleinit + rvec);
            saddle(3,1)=-xback(3,1);
            xback(3,1)=0;
            for i=1:numwheel
                wheelcoord(:,i)=Totalrotationmatrix*wheelcoord(:,i);
            end
            centerwheel = Totalrotationmatrix*centerwheel;
            
            %translations
            wheelcoord(1,:) = wheelcoord(1,:) +xback(1,1)*ones(1,numwheel);
            wheelcoord(2,:) = wheelcoord(2,:) +xback(2,1)*ones(1,numwheel);
            wheelcoord(3,:) = wheelcoord(3,:) +xback(3,1)*ones(1,numwheel);
            centerwheel(1,1)=centerwheel(1,1)+xback(1,1);
            centerwheel(2,1)=centerwheel(2,1)+xback(2,1);
            centerwheel(3,1)=centerwheel(3,1)+xback(3,1);
            %collecting terms
            timeddata(j,:,:)=transpose([wheelcoord,centerwheel,saddle]);
    
        end
    elseif strcmp(xlocs,'bike')
        z =zeros(length(x));
        for j=1:numtime
            centerwheel=zeros(3,1);
            centerwheel(3,1)=r;
            saddleinit = zeros(3,1);
            saddleinit(3,1)=h;
            wheelcoord = zeros(3,numwheel);
            wheelcoord(1,:) = r*cos(angles);
            wheelcoord(3,:)= r+r*sin(angles);
            %rotations
            Rotationmatrixepsilon = Rotationmatrix(epsilon(j),'x');
            Rotationmatrixtheta = Rotationmatrix(theta(j),'z');
            Rotationmatrixalpha = Rotationmatrix(alpha(j),'y');
            Totalrotationmatrix = Rotationmatrixtheta*Rotationmatrixepsilon;
            for i=1:numwheel
                wheelcoord(:,i)=Totalrotationmatrix*wheelcoord(:,i);
            end
            centerwheel = Totalrotationmatrix*centerwheel;
            saddle = Rotationmatrixalpha*saddleinit;
            saddle(3,1)=saddle(3,1)+r; %rotation point of saddle was around center of wheel so it now first still needs to be translated upward
            saddle=Totalrotationmatrix*saddle;
            
            %translations
            wheelcoord(1,:) = wheelcoord(1,:) +x(j)*ones(1,numwheel);
            wheelcoord(2,:) = wheelcoord(2,:) +y(j)*ones(1,numwheel);
            wheelcoord(3,:) = wheelcoord(3,:) +z(j)*ones(1,numwheel);
            centerwheel(1,1)=centerwheel(1,1)+x(j);
            centerwheel(2,1)=centerwheel(2,1)+y(j);
            centerwheel(3,1)=centerwheel(3,1)+z(j);
            saddle(1,1)=saddle(1,1)+x(j);
            saddle(2,1)=saddle(2,1)+y(j);
            saddle(3,1)=saddle(3,1)+z(j);
            %collecting terms
            timeddata(j,:,:)=transpose([wheelcoord,centerwheel,saddle]);
            %xback(j,:) = saddle - Totalrotationmatrix*(Rotationmatrixalpha*saddleinit + rvec);
        end
    end
    vec = timeddata;
    
    %plotting
    %myVideo = VideoWriter('myVideoFile'); %open video file
    %myVideo.FrameRate = 10;  %can adjust this, 5 - 10 works well for me
    %open(myVideo)
    for i=1:numtime
        centerwheel=timeddata(i,numwheel+1,:);
        saddle = timeddata(i,numwheel+2,:);
        wheelcoord = timeddata(i,1:numwheel,:);
        saddleplot=[centerwheel,saddle];
        plot3(wheelcoord(1,:,1),wheelcoord(1,:,2),wheelcoord(1,:,3))
        hold on
        plot3(centerwheel(1,1,1),centerwheel(1,1,2),centerwheel(1,1,3),'.')
        plot3(saddleplot(1,:,1),saddleplot(1,:,2),saddleplot(1,:,3),'-o')
        xlim([-5,5])
        ylim([-5,5])
        zlim([-3,3])
        xlabel('x')
        ylabel('y')
        zlabel('z')
        drawnow
        pause(0.01)
        %frame = getframe(gcf); %get frame
        %writeVideo(myVideo, frame);
        hold off
    end
    %close(myVideo)
end

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


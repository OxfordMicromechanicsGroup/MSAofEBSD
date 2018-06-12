%% Simulation of Gnomonic projection from Stereographic image -- Tomohito
function [SimPat] = zxz_tilt_sim1(angles,PC,Width,Height,MasterPattern,tilt)

% angles - Euler angles in degrees
% PC - pattern centre geometry
% Width - output simulated pattern width
% Height - output simulated pattern height
% Master Pattern - full simulated pattern on stereographic projection from Dynamics
% tilt - effective tilt degs - accounts for sample and camera tilts
%
%   AJW 8/5/18


k=length(MasterPattern);

phi1=angles(1);
PHI=angles(2);
phi2=angles(3);
Xstar=PC(1); %pattern centre X position as fraction of pattern width
Ystar=PC(2); %pattern centre Y position as fraction of pattern height
Zstar=PC(3);


%eq 8 in 'which way is up'
Rz=@(theta)[cosd(theta) sind(theta) 0;-sind(theta) cosd(theta) 0;0 0 1];

%eq 9 in 'which way is up'
Rx=@(theta)[1 0 0;0 cosd(theta) sind(theta);0 -sind(theta) cosd(theta)];

% crystal orientation matrix from Euler angles 
OriMatrix=Rz(phi2)*Rx(PHI)*Rz(phi1);  % - eq10 in 'which way is up'
crystal2detector=OriMatrix*Rx(tilt);  % - eq17 in 'which way is up'
detector2crystal=pinv(crystal2detector);

rot=detector2crystal;

Sphere_Radius = Zstar*Height;
% Sphere_Radius = Zstar*(Height-1);

[x_gnomon,y_gnomon]=meshgrid(1:Width,1:Height);
% posX = x_gnomon(:)-(Xstar*(Width-1)+1);
% posY = y_gnomon(:)-(Ystar*(Height-1)+1);
posX = x_gnomon(:)-(Xstar*Width);
posY = y_gnomon(:)-(Ystar*Height);
vec = [posX,posY,Sphere_Radius*ones(Width*Height,1)];
Vector_On_Sphere = vec./sqrt(sum(abs(vec).^2,2)); %vec(:,1:3) ./ norm(vec(:,1:3));


Vector_On_Sphere_original=Vector_On_Sphere*detector2crystal;
% Vector_On_Sphere_original=pinv(rot*eye(3))*Vector_On_Sphere';
% Vector_On_Sphere_original=Vector_On_Sphere_original';

Vector_On_Stereo=Vector_On_Sphere_original(:,1:2)./(1+abs(Vector_On_Sphere_original(:,3)));

Vector_On_Stereo_shift = (k-1)/2*(Vector_On_Stereo(:,1:2)) + (k+1)/2;

F=griddedInterpolant(MasterPattern,'cubic');

% Create intensities
StereoInt=F(Vector_On_Stereo_shift(:,2),Vector_On_Stereo_shift(:,1));
SimPat = reshape(StereoInt,[Height,Width]);


end


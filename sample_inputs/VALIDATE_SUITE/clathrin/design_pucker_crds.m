%angle1 is in rads. It is the angle between the clathrin leg and the
%z-axis, so a flat clathin is pi/2. A puckered clathrin has an angle>pi/2.
%legLen1 is for once they are assembled, so twice a single leg:
%Prints clathrin-clathrin binding sites at 1/2 of the leg.
%legLen2 determines how far the AP2 sites are from the clathrin sites
%located at leglen1. So it is designed to mimic the clathrin crystal
%structure, but wiht the cla-cla bindin sites printed at half the leg
%length so that when they bind it is the full leg len (+ sigma)

function[RigidCrds]=design_pucker_crds(angle1, legLen1, legLen2)
format long
angle2=120*pi/180;%need three-fold symmetry, use 120 degrees.
Mx=rotateEulerX(angle1);

%rotate a vector from the z-axis to angle 1 degrees from it.
v1=[0; 0; legLen1];

vRot=Mx*v1;
%length after rotation:
display('length of leg')
sqrt(vRot'*vRot)
%Now rotate vRot around the Z-axis, to create the lower leg. 

Mz=rotateEulerZ(angle2);

lowerleg=Mz*vRot;
%shift so that the end of lower leg is now the origin
origin=[0;0;0]
tmpCrds=[lowerleg, origin, vRot];
display('Current lengths')
sqrt(tmpCrds(:,1)'*tmpCrds(:,1))
d=tmpCrds(:,3);
sqrt(d'*d)

%angle between lowerleg and vrot.
l1=sqrt(lowerleg'*lowerleg)
l2=sqrt(vRot'*vRot)
display('Angle between lowerleg and upper');
theta=acos(lowerleg'*vRot/(l1*l2))

%angle between lowerleg and Z.
zAxis=[0;0;1];
display('Angle between lowerleg and Z axis');
theta=acos(lowerleg'*zAxis/(l1))

%angle between lowerleg and Z.
zAxis=[0;0;1];
display('Angle between upper leg and Z-axis');
theta=acos(vRot'*zAxis/(l2))


%Normal to upper and lower leg (should be z-axis!)
n1=cross(vRot, lowerleg);
nlen=sqrt(n1'*n1);
n1=n1/nlen

figure(1)
hold on;
p1=[origin';vRot'];
%plot3(p1(:,1),p1(:,2),p1(:,3),'ro-','LineWidth',3);

p1=[origin';lowerleg'];
%plot3(p1(:,1),p1(:,2),p1(:,3),'bo-','LineWidth',3);

%subtract off the lowerleg crd to move that to origin
off=lowerleg.*ones(3,3);
legCrds=tmpCrds-off;

%After translation, length of legs.
display('Current lengths after translation')
sqrt(legCrds(:,2)'*legCrds(:,2))
d=legCrds(:,3)-legCrds(:,2);
sqrt(d'*d)
%angle between distal and the z-axis. 

display('After Translation, angle between upper leg and Z-axis');
theta=acos(legCrds(:,2)'*zAxis/sqrt(legCrds(:,2)'*legCrds(:,2)))

%normal to this vector and z-axis.
n1=cross(zAxis, legCrds(:,2));

n1=n1/sqrt(n1'*n1);
rotAngle=angle1-theta;
display('angle to rotate')
rotAngle
n1
M=rotation_axis_angle(rotAngle, n1);
%Now apply M to both points in the legCrds.
j1=M*legCrds(:,2);
j2=M*legCrds(:,3);
newCrds=[legCrds(:,1),j1, j2]
display('Current lengths after Rotation')
sqrt(newCrds(:,2)'*newCrds(:,2))
d=newCrds(:,3)-newCrds(:,2);
sqrt(d'*d)

p1=[origin';zAxis'];
plot3(p1(:,1),p1(:,2),p1(:,3),'ko-','LineWidth',3);

p1=[origin';newCrds(:,2)'];
plot3(p1(:,1),p1(:,2),p1(:,3),'mo-','LineWidth',3);



%Now rotate this entire leg (two points) around the z axis by angle2 each
%time.

Mz=rotateEulerZ(angle2);
leg2=zeros(3,2);
leg3=zeros(3,2);
leg2(:,1)=Mz*newCrds(:,2);
leg2(:,2)=Mz*newCrds(:,3);
p1=[origin';leg2(:,1)'];
plot3(p1(:,1),p1(:,2),p1(:,3),'mo-','LineWidth',3);



%Now rotate leg2 by the same angle
leg3(:,1)=Mz*leg2(:,1);
leg3(:,2)=Mz*leg2(:,2);
p1=[origin';leg3(:,1)'];
plot3(p1(:,1),p1(:,2),p1(:,3),'mo-','LineWidth',3);




TotCrds=[newCrds, leg2, leg3];
%Arrange the coordinates for the input files.
%Put clathrin-clathrin binding sites at 1/2 of the leg.
display('clath leg length')
sqrt(newCrds(:,2)'*newCrds(:,2))
sqrt(leg2(:,1)'*leg2(:,1))
cc1=newCrds(:,2)*0.5;
cc2=leg2(:,1)*0.5;
cc3=leg3(:,1)*0.5;

%PLOT clath binding sites
plot3(cc1(1),cc1(2),cc1(3),'bo','MarkerSize',10);
plot3(cc2(1),cc2(2),cc2(3),'bo','MarkerSize',10);
plot3(cc3(1),cc3(2),cc3(3),'bo','MarkerSize',10);
%put AP2 binding site at leg2 distance from the knee. Currently it is at
%leg1 distance. 
ratio=legLen2/legLen1;
aptmp=(newCrds(:,3)-newCrds(:,2))*ratio;
ac1=newCrds(:,2)+aptmp;
aptmp=(leg2(:,2)-leg2(:,1))*ratio;
ac2=leg2(:,1)+aptmp;
aptmp=(leg3(:,2)-leg3(:,1))*ratio;
ac3=leg3(:,1)+aptmp;


%plot AP crds;

p1=[newCrds(:,2)';ac1'];
plot3(p1(:,1),p1(:,2),p1(:,3),'co-','LineWidth',3);

p1=[leg2(:,1)';ac2'];
plot3(p1(:,1),p1(:,2),p1(:,3),'co-','LineWidth',3);

p1=[leg3(:,1)';ac3'];
plot3(p1(:,1),p1(:,2),p1(:,3),'co-','LineWidth',3);

%Draw lines from CC to AP.

p1=[cc1';ac1'];
plot3(p1(:,1),p1(:,2),p1(:,3),'yo-','LineWidth',3);

p1=[cc2';ac2'];
plot3(p1(:,1),p1(:,2),p1(:,3),'yo-','LineWidth',3);

p1=[cc3';ac3'];
plot3(p1(:,1),p1(:,2),p1(:,3),'yo-','LineWidth',3);

RigidCrds=[newCrds(:,1)';cc1';cc2';cc3';ac1';ac2';ac3'];

%Angle of COM-AP binding site with the x-y plane:
v=ac1-newCrds(:,1);
lenV=sqrt(v'*v);
v2=v;
v2(3)=0;%v in the x-y plane. 
lenV2=sqrt(v2'*v2);
theta=acos(v'*v2/lenV/lenV2)


%%create rotation matrix for angle (in rads) around the X axis. 
function[M]=rotateEulerX(angle)
format long
s=sin(angle);
c=cos(angle);

M=[1 0 0 ;
    0 c -s;
    0 s c]


%create rotation matrix for angle (in rads) around the z axis. 
function[M]=rotateEulerZ(angle)
format long
s=sin(angle);
c=cos(angle);

M=[c -s 0 ;
    s c 0;
    0 0 1]


%perform rotation of theta around an axis u
function[M]=rotation_axis_angle(theta, u)
format long
ux=u(1);
uy=u(2);
uz=u(3);
ct=cos(theta);
st=sin(theta);
m=1-ct;

M(1,1)=ct+ux^2*m;
M(1,2)=ux*uy*m-uz*st;
M(1,3)=ux*uz*m+uy*st;
M(2,1)=uy*ux*m+uz*st;
M(2,2)=ct+uy^2*m;
M(2,3)=uy*uz*m-ux*st;
M(3,1)=uz*ux*m-uy*st;
M(3,2)=uz*uy*m+ux*st;
M(3,3)=ct+uz^2*m;


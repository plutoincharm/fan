%%% assuming ankle is at post, foot end is meta
clc;
clear;
clearvars;
%%%
global L rc m MI g
%%%
L1 = 0.4418;
L2 = 0.4687;
L3 = 0.2035;
r1 = 0.1872;
r2 = 0.1738;
r3 = 0.09;
m1 = 7.4;
m2 = 3.411;
m3 = 1.073;
Ixx1 = 0.1038;
Iyy1 = 0.0;
Izz1 = 0.1038;
Ixx2 = 0.0;
Iyy2 = 0.05916;
Izz2 = 0.05916;
Ixx3 = 0.0;
Iyy3 = 0.01;
Izz3 = 0.01;
L = [L1 L2 L3];
rc = [r1 r2 r3];
m = [m1 m2 m3];
MI = [Ixx1 Iyy1 Izz1 Ixx2 Iyy2 Izz2 Ixx3 Iyy3 Izz3];
g = 9.81;
%%%
coord = readmatrix("JointCoords3D_edit.xlsx");
time = coord(:,1);
% obtaining GRF and GRT data for the four stages of the gait cycle
GRdata = readmatrix("GRdata.xlsx");
t1 = time(24:76);
f1 = GRdata(24:76,:);
t2 = time(77:91);
f2 = GRdata(77:91,:);
t3 = time(92:144);
f3 = GRdata(92:144,:);
t4 = time(145:159);
f41 = GRdata(77:91,2:7);
f42 = GRdata(77:91,8:13);
f4 = [t4 f42 f41];
time_plot = [t1;t2;t3;t4];
GRdata_new = [f1;f2;f3;f4];
%%%
gtc = 24:159; % gait cycle time is from 0.23 to 1.58 seconds
tfit = time(gtc);
dt = 0.01; % dt should be set to 0.01 only
tt = 0.23:dt:1.58;
GRFx = -GRdata_new(:,2);
GRFy = GRdata_new(:,3);
GRFz = GRdata_new(:,4);
%%% right side
xhip = 2-0.001*coord(gtc,20);
yhip = 0.001*coord(gtc,21);
zhip = 0.001*coord(gtc,22);
xknee = 2-0.001*coord(gtc,23);
yknee = 0.001*coord(gtc,24);
zknee = 0.001*coord(gtc,25);
xankleold = 2-0.001*coord(gtc,26);
yankleold = 0.001*coord(gtc,27);
zankleold = 0.001*coord(gtc,28);
xpost = 2-0.001*coord(gtc,29);
ypost = 0.001*coord(gtc,30);
zpost = 0.001*coord(gtc,31);
xmeta = 2-0.001*coord(gtc,32);
ymeta = 0.001*coord(gtc,33);
zmeta = 0.001*coord(gtc,34);
xcop = 2-0.001*coord(gtc,35);
ycop = 0.001*coord(gtc,36)+0.02;
zcop = 0.001*coord(gtc,37);
dcop = sqrt((xcop-xpost).^2 + (ycop-ypost).^2 + (zcop-zpost).^2);
rGRF = r3 - dcop; % location of cop wrt post
%%%
xankle = xpost;
yankle = ypost;
zankle = zpost;
xfoot = xmeta;
yfoot = ymeta;
zfoot = zmeta;
%%%
Len_thigh = mean(sqrt((xhip-xknee).^2 + (yhip-yknee).^2 + (zhip-zknee).^2));
Len_shank = mean(sqrt((xknee-xankle).^2 + (yknee-yankle).^2 + (zknee-zankle).^2));
Len_foot = mean(sqrt((xankle-xfoot).^2 + (yankle-yfoot).^2 + (zankle-zfoot).^2));
%%% Inverse kinematics
tht2 = acos((zhip-zknee)/L1);
tht1 = atan2((yhip-yknee)./(L1*sin(tht2)),(xhip-xknee)./(L1*sin(tht2)));
s4 = L2.^(-1).*(zankle.*cos(tht2)+(-1).*zknee.*cos(tht2)+xankle.*cos( ...
  tht1).*sin(tht2)+(-1).*xknee.*cos(tht1).*sin(tht2)+yankle.*sin( ...
  tht1).*sin(tht2)+(-1).*yknee.*sin(tht1).*sin(tht2));
tht4 = asin(s4);
s3c4 = L2.^(-1).*(((-1).*yankle+yknee).*cos(tht1)+xankle.*sin(tht1)+(-1) ...
  .*xknee.*sin(tht1));
c3c4 = L2.^(-1).*(xankle.*cos(tht1).*cos(tht2)+(-1).*xknee.*cos(tht1).* ...
  cos(tht2)+yankle.*cos(tht2).*sin(tht1)+(-1).*yknee.*cos(tht2).* ...
  sin(tht1)+(-1).*zankle.*sin(tht2)+zknee.*sin(tht2));
tht3 = atan2(s3c4./cos(tht4),c3c4./cos(tht4));
s6 = L3.^(-1).*((-1).*yankle.*cos(tht1).*cos(tht3)+yfoot.*cos(tht1).* ...
  cos(tht3)+(-1).*yankle.*cos(tht2).*sin(tht1).*sin(tht3)+yfoot.* ...
  cos(tht2).*sin(tht1).*sin(tht3)+zankle.*sin(tht2).*sin(tht3)+(-1) ...
  .*zfoot.*sin(tht2).*sin(tht3)+xankle.*(cos(tht3).*sin(tht1)+(-1).* ...
  cos(tht1).*cos(tht2).*sin(tht3))+xfoot.*((-1).*cos(tht3).*sin( ...
  tht1)+cos(tht1).*cos(tht2).*sin(tht3)));
tht6 = asin(s6);
s5c6 = L3.^(-1).*((-1).*zankle.*cos(tht2).*cos(tht4)+zfoot.*cos(tht2).* ...
  cos(tht4)+(-1).*yankle.*cos(tht4).*sin(tht1).*sin(tht2)+yfoot.* ...
  cos(tht4).*sin(tht1).*sin(tht2)+yankle.*cos(tht2).*cos(tht3).*sin( ...
  tht1).*sin(tht4)+(-1).*yfoot.*cos(tht2).*cos(tht3).*sin(tht1).* ...
  sin(tht4)+(-1).*zankle.*cos(tht3).*sin(tht2).*sin(tht4)+zfoot.* ...
  cos(tht3).*sin(tht2).*sin(tht4)+(-1).*yankle.*cos(tht1).*sin(tht3) ...
  .*sin(tht4)+yfoot.*cos(tht1).*sin(tht3).*sin(tht4)+xfoot.*((-1).* ...
  sin(tht1).*sin(tht3).*sin(tht4)+cos(tht1).*(cos(tht4).*sin(tht2)+( ...
  -1).*cos(tht2).*cos(tht3).*sin(tht4)))+xankle.*(sin(tht1).*sin( ...
  tht3).*sin(tht4)+cos(tht1).*((-1).*cos(tht4).*sin(tht2)+cos(tht2) ...
  .*cos(tht3).*sin(tht4))));
c5c6 = L3.^(-1).*((-1).*yankle.*cos(tht2).*cos(tht3).*cos(tht4).*sin( ...
  tht1)+yfoot.*cos(tht2).*cos(tht3).*cos(tht4).*sin(tht1)+zankle.* ...
  cos(tht3).*cos(tht4).*sin(tht2)+(-1).*zfoot.*cos(tht3).*cos(tht4) ...
  .*sin(tht2)+yankle.*cos(tht1).*cos(tht4).*sin(tht3)+(-1).*yfoot.* ...
  cos(tht1).*cos(tht4).*sin(tht3)+(-1).*zankle.*cos(tht2).*sin(tht4) ...
  +zfoot.*cos(tht2).*sin(tht4)+(-1).*yankle.*sin(tht1).*sin(tht2).* ...
  sin(tht4)+yfoot.*sin(tht1).*sin(tht2).*sin(tht4)+(-1).*xankle.*( ...
  cos(tht4).*sin(tht1).*sin(tht3)+cos(tht1).*(cos(tht2).*cos(tht3).* ...
  cos(tht4)+sin(tht2).*sin(tht4)))+xfoot.*(cos(tht4).*sin(tht1).* ...
  sin(tht3)+cos(tht1).*(cos(tht2).*cos(tht3).*cos(tht4)+sin(tht2).* ...
  sin(tht4))));
tht5 = atan2(s5c6./cos(tht6),c5c6./cos(tht6));
%%% fitting fourier series for theta
f1 = fit(tfit,tht1,'fourier8');
tht1 = f1(tt);
f2 = fit(tfit,tht2,'fourier8');
tht2 = f2(tt);
f3 = fit(tfit,tht3,'fourier8');
tht3 = f3(tt);
f4 = fit(tfit,tht4,'fourier8');
tht4 = f4(tt);
f5 = fit(tfit,tht5,'fourier8');
tht5 = f5(tt);
f6 = fit(tfit,tht6,'fourier8');
tht6 = f6(tt);
fxhip = fit(tfit,xhip,'fourier8');
xhip = fxhip(tt);
fyhip = fit(tfit,yhip,'fourier8');
yhip = fyhip(tt);
fzhip = fit(tfit,zhip,'fourier8');
zhip = fzhip(tt);
%%% acceleration of hip, and time derivatives of thetas
[eps1,eta1] = differentiate(f1,tt);
[eps2,eta2] = differentiate(f2,tt);
[eps3,eta3] = differentiate(f3,tt);
[eps4,eta4] = differentiate(f4,tt);
[eps5,eta5] = differentiate(f5,tt);
[eps6,eta6] = differentiate(f6,tt);
[vxhip,axhip] = differentiate(fxhip,tt);
[vyhip,ayhip] = differentiate(fyhip,tt);
[vzhip,azhip] = differentiate(fzhip,tt);
%%% estimating joint torques
c = 1;
sol = zeros(length(tt),6); %% preallocating
for k = 1:length(tt)
    q = [xhip(k) yhip(k) zhip(k) tht1(k) tht2(k) tht3(k) tht4(k) tht5(k) tht6(k)];
    qdot = [vxhip(k) vyhip(k) vzhip(k) eps1(k) eps2(k) eps3(k) eps4(k) eps5(k) eps6(k)];
    qddot = [axhip(k) ayhip(k) azhip(k) eta1(k) eta2(k) eta3(k) eta4(k) eta5(k) eta6(k)];
    GRF = [GRFx(k) GRFy(k) GRFz(k)];
    r = rGRF(k);
    sol(c,:) = torques3D_II(q,qdot,qddot,GRF,r);
    c = c+1;
end
T4 = sol(:,1);
T5 = sol(:,2);
T6 = sol(:,3);
T7 = sol(:,4);
T8 = sol(:,5);
T9 = sol(:,6);
%%% plotting
figure;
plot(tt,tht1,'b-','LineWidth',1);
grid on;
title('\theta_1 vs time');
xlabel('time (sec) \rightarrow');
ylabel('\theta_1 (rad) \rightarrow');
figure;
plot(tt,tht2,'b-','LineWidth',1);
grid on;
title('\theta_2 vs time');
xlabel('time (sec) \rightarrow');
ylabel('\theta_2 (rad) \rightarrow');
figure;
plot(tt,tht3,'b-','LineWidth',1);
grid on;
title('\theta_3 vs time');
xlabel('time (sec) \rightarrow');
ylabel('\theta_3 (rad) \rightarrow');
figure;
plot(tt,tht4,'b-','LineWidth',1);
grid on;
title('\theta_4 vs time');
xlabel('time (sec) \rightarrow');
ylabel('\theta_4 (rad) \rightarrow');
figure;
plot(tt,tht5,'b-','LineWidth',1);
grid on;
title('\theta_5 vs time');
xlabel('time (sec) \rightarrow');
ylabel('\theta_5 (rad) \rightarrow');
figure;
plot(tt,tht6,'b-','LineWidth',1);
grid on;
title('\theta_6 vs time');
xlabel('time (sec) \rightarrow');
ylabel('\theta_6 (rad) \rightarrow');
%
figure;
plot(tt,T4,'b-','LineWidth',1);
grid on;
title('Hip flexion/extension moment vs time');
xlabel('time (sec) \rightarrow');
ylabel('T_4 (N-m) \rightarrow');
figure;
plot(tt,T5,'b-','LineWidth',1);
grid on;
title('Hip abduction/adduction moment vs time');
xlabel('time (sec) \rightarrow');
ylabel('T_5 (N-m) \rightarrow');
figure;
plot(tt,T6,'b-','LineWidth',1);
grid on;
title('Hip internal/external rotation moment vs time');
xlabel('time (sec) \rightarrow');
ylabel('T_6 (N-m) \rightarrow');
figure;
plot(tt,T7,'b-','LineWidth',1);
grid on;
title('Knee flexion/extension moment vs time');
xlabel('time (sec) \rightarrow');
ylabel('T_7 (N-m) \rightarrow');
figure;
plot(tt,T8,'b-','LineWidth',1);
grid on;
title('Ankle dorsiflexion/plantarflexion moment vs time');
xlabel('time (sec) \rightarrow');
ylabel('T_8 (N-m) \rightarrow');
figure;
plot(tt,T9,'b-','LineWidth',1);
grid on;
title('Ankle inversion/eversion moment vs time');
xlabel('time (sec) \rightarrow');
ylabel('T_9 (N-m) \rightarrow');
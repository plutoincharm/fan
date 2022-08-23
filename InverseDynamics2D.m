clc;
clear;
clearvars;
%%%
global L r m MI g
%%%
L1 = 0.4418;
L2 = 0.4033;
r1 = 0.1872;
r2 = 0.1738;
r3 = 0.044;
m1 = 7.4;
m2 = 3.411;
m3 = 1.073;
MI1 = 0.1038;
MI2 = 0.05916;
MI3 = 0.01;
L = [L1 L2];
r = [r1 r2 r3];
m = [m1 m2 m3];
MI = [MI1 MI2 MI3];
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
%%% right side
xhip = 2-0.001*coord(gtc,20);
yhip = 0.001*coord(gtc,21);
xknee = 2-0.001*coord(gtc,23);
yknee = 0.001*coord(gtc,24);
xankle = 2-0.001*coord(gtc,26);
yankle = 0.001*coord(gtc,27);
xpost = 2-0.001*coord(gtc,29);
ypost = 0.001*coord(gtc,30);
xmeta = 2-0.001*coord(gtc,32);
ymeta = 0.001*coord(gtc,33);
xcop = 2-0.001*coord(gtc,35);
ycop = 0.001*coord(gtc,36)+0.02;
x3 = (xankle+xpost+xmeta)/3;
y3 = (yankle+ypost+ymeta)/3;
%%% Inverse kinematics
tht1 = pi/2 + atan2(yhip-yknee,xhip-xknee);
tht2 = pi/2 + atan2(yknee-yankle,xknee-xankle);
tht3 = pi/2 + atan2(yankle-y3,xankle-x3);
%%% fitting fourier series for theta
f1 = fit(tfit,tht1,'fourier8');
tht1 = f1(tt);
f2 = fit(tfit,tht2,'fourier8');
tht2 = f2(tt);
f3 = fit(tfit,tht3,'fourier8');
tht3 = f3(tt);
fxhip = fit(tfit,xhip,'fourier8');
xhip = fxhip(tt);
fyhip = fit(tfit,yhip,'fourier8');
yhip = fyhip(tt);
%%% acceleration of hip, and time derivatives of thetas
[omg1,alp1] = differentiate(f1,tt);
[omg2,alp2] = differentiate(f2,tt);
[omg3,alp3] = differentiate(f3,tt);
[vxhip,axhip] = differentiate(fxhip,tt);
[vyhip,ayhip] = differentiate(fyhip,tt);
%%%
Len_thigh = mean(sqrt((xhip-xknee).^2 + (yhip-yknee).^2));
Len_shank = mean(sqrt((xknee-xankle).^2 + (yknee-yankle).^2));
%%% estimating joint torques
c = 1;
sol = zeros(length(tt),3); %% preallocating
for k = 1:length(tt)
    q = [xhip(k) yhip(k) tht1(k) tht2(k) tht3(k)];
    qdot = [vxhip(k) vyhip(k) omg1(k) omg2(k) omg3(k)];
    qddot = [axhip(k) ayhip(k) alp1(k) alp2(k) alp3(k)];
    GRF = [GRFx(k) GRFy(k)];
    cop = [xcop(k) ycop(k)];
    sol(c,:) = torques2D(q,qdot,qddot,GRF,cop);
    c = c+1;
end
T1 = sol(:,1);
T2 = sol(:,2);
T3 = sol(:,3);
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
plot(tt,T1,'b-','LineWidth',1);
grid on;
title('Hip torque vs time');
xlabel('time (sec) \rightarrow');
ylabel('T_1 (N-m) \rightarrow');
figure;
plot(tt,T2,'b-','LineWidth',1);
grid on;
title('Knee torque vs time');
xlabel('time (sec) \rightarrow');
ylabel('T_2 (N-m) \rightarrow');
figure;
plot(tt,T3,'b-','LineWidth',1);
grid on;
title('Ankle torque vs time');
xlabel('time (sec) \rightarrow');
ylabel('T_3 (N-m) \rightarrow');
clc;
clearvars -except X1 Y1x Y1y;
%%%
coord = readmatrix("gait_data_2.xlsx");
time = coord(:,1);
tfit = time;
tt = time;
%%%
xtrunk = 0.001*coord(:,2);
ytrunk = 0.001*coord(:,3);
%%% left side
lxhip = 0.001*coord(:,4);
lyhip = 0.001*coord(:,5);
lxknee = 0.5*(0.001*coord(:,6)+0.001*coord(:,8));
lyknee = 0.5*(0.001*coord(:,7)+0.001*coord(:,9));
lxankle = 0.5*(0.001*coord(:,10)+0.001*coord(:,12));
lyankle = 0.5*(0.001*coord(:,11)+0.001*coord(:,13));
lxpost = 0.001*coord(:,14);
lypost = 0.001*coord(:,15);
lxmeta = 0.001*coord(:,16);
lymeta = 0.001*coord(:,17);
%%% right side
rxhip = 0.001*coord(:,18);
ryhip = 0.001*coord(:,19);
rxknee = 0.5*(0.001*coord(:,20)+0.001*coord(:,22));
ryknee = 0.5*(0.001*coord(:,21)+0.001*coord(:,23));
rxankle = 0.5*(0.001*coord(:,24)+0.001*coord(:,26));
ryankle = 0.5*(0.001*coord(:,25)+0.001*coord(:,27));
rxpost = 0.001*coord(:,28);
rypost = 0.001*coord(:,29);
rxmeta = 0.001*coord(:,30);
rymeta = 0.001*coord(:,31);
%%% mid-point of pelvis
x_p = (lxhip+rxhip)/2;
y_p = (lyhip+ryhip)/2;
rxfoot = 1/3*(rxankle + rxpost + rxmeta);
ryfoot = 1/3*(ryankle + rypost + rymeta);
lxfoot = 1/3*(lxankle + lxpost + lxmeta);
lyfoot = 1/3*(lyankle + lypost + lymeta);
%%% angles for inverse dynamics code
x1 = xtrunk;
y1 = ytrunk;
tht1 = pi + atan2(ytrunk-y_p,xtrunk-x_p);
tht2 = pi + atan2(ryhip-ryknee,rxhip-rxknee);
tht3 = pi + atan2(lyhip-lyknee,lxhip-lxknee);
tht4 = pi + atan2(ryknee-ryankle,rxknee-rxankle);
tht5 = pi + atan2(lyknee-lyankle,lxknee-lxankle);
tht6 = pi + atan2(ryankle-ryfoot,rxankle-rxfoot);
tht7 = pi + atan2(lyankle-lyfoot,lxankle-lxfoot);
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
f7 = fit(tfit,tht7,'fourier8');
tht7 = f7(tt);
fx1 = fit(tfit,x1,'fourier8');
x1 = fx1(tt);
fy1 = fit(tfit,y1,'fourier8');
y1 = fy1(tt);
%%% accelerations
[vx1,ax1] = differentiate(fx1,tt);
[vy1,ay1] = differentiate(fy1,tt);
[omg1,alp1] = differentiate(f1,tt);
[omg2,alp2] = differentiate(f2,tt);
[omg3,alp3] = differentiate(f3,tt);
[omg4,alp4] = differentiate(f4,tt);
[omg5,alp5] = differentiate(f5,tt);
[omg6,alp6] = differentiate(f6,tt);
[omg7,alp7] = differentiate(f7,tt);
%%% GRFs
f = readmatrix("grf_data_2.xlsx");
tf = f(:,1);
gt = 1:15:2311;
GRFx = f(gt,2);
GRFy = f(gt,3);
mark = 6:102; % full gait cycle
%mark = 6:67; % non zero GRF
%mark = 6:16; % first DS
%mark = 53:67; % second DS
%%%
X = [vx1 vy1 omg1 omg2 omg3 omg4 omg5...
    omg6 omg7 ax1 ay1 alp1 alp2 alp3 alp4 alp5 alp6 alp7];
Yx = GRFx/67;
Yy = GRFy/67;
X2 = X(mark,:);
Y2x = Yx(mark,:);
Y2y = Yy(mark,:);
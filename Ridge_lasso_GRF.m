clc;
clearvars -except X1 Y1x Y1y X2 Y2x Y2y X3 Y3x Y3y X4 Y4x Y4y X5 Y5x Y5y X6 Y6x Y6y betax coefx coef0x betay coefy coef0y;
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
%
xtrunk = 2-0.001*coord(gtc,2);
ytrunk = 0.001*coord(gtc,3);
%%% left side
lxhip = 2-0.001*coord(gtc,5);
lyhip = 0.001*coord(gtc,6);
lxknee = 2-0.001*coord(gtc,8);
lyknee = 0.001*coord(gtc,9);
lxankle = 2-0.001*coord(gtc,11);
lyankle = 0.001*coord(gtc,12);
lxpost = 2-0.001*coord(gtc,14);
lypost = 0.001*coord(gtc,15);
lxmeta = 2-0.001*coord(gtc,17);
lymeta = 0.001*coord(gtc,18);
%%% right side
rxhip = 2-0.001*coord(gtc,20);
ryhip = 0.001*coord(gtc,21);
rxknee = 2-0.001*coord(gtc,23);
ryknee = 0.001*coord(gtc,24);
rxankle = 2-0.001*coord(gtc,26);
ryankle = 0.001*coord(gtc,27);
rxpost = 2-0.001*coord(gtc,29);
rypost = 0.001*coord(gtc,30);
rxmeta = 2-0.001*coord(gtc,32);
rymeta = 0.001*coord(gtc,33);
%%% mid-point of pelvis
x_p = (lxhip+rxhip)/2;
y_p = (lyhip+ryhip)/2;
rxfoot = 1/3*(rxankle + rxpost + rxmeta);
ryfoot = 1/3*(ryankle + rypost + rymeta);
lxfoot = 1/3*(lxankle + lxpost + lxmeta);
lyfoot = 1/3*(lyankle + lypost + lymeta);
%
R2 = mean(sqrt((rxhip-xtrunk).^2 + (ryhip-ytrunk).^2));
R3 = mean(sqrt((lxhip-xtrunk).^2 + (lyhip-ytrunk).^2));
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
%%% misc
del2 = mean(tht1 - atan2(y1-ryhip,x1-rxhip) - pi);
del3 = mean(pi - tht1 + atan2(y1-lyhip,x1-lxhip));
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
%%% forces
X = [vx1 vy1 omg1 omg2 omg3 omg4 omg5...
    omg6 omg7 ax1 ay1 alp1 alp2 alp3 alp4 alp5 alp6 alp7];
Yx_ridge = betax(1) + X*betax(2:end); % ridge
Yx_lasso = X*coefx + coef0x;  % lasso
GRFx_ridge = Yx_ridge*73.94;
GRFx_lasso = Yx_lasso*73.94;
%
Yy_ridge = betay(1) + X*betay(2:end); % ridge
Yy_lasso = X*coefy + coef0y;  % lasso
GRFy_ridge = Yy_ridge*73.94;
GRFy_lasso = Yy_lasso*73.94;

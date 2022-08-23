%%% assuming ankle is at post, foot end is meta
clc;
clearvars -except GRFxr GRFyr GRFxl GRFyl;
%%%
global L rc m g misc MI
%%%
L1 = 0.4418;
L2 = 0.4687;
L3 = 0.2035;
r1 = 0.1872;
r2 = 0.1738;
r3 = 0.09;
m1 = 50.172;
m2 = 7.4;
m3 = 3.411;
m4 = 1.073;
m5 = 7.4;
m6 = 3.411;
m7 = 1.073;
Ixx2 = 0.1038;
Iyy2 = 0.0;
Izz2 = 0.1038;
Ixx3 = 0.0;
Iyy3 = 0.05916;
Izz3 = 0.05916;
Ixx4 = 0.0;
Iyy4 = 0.01;
Izz4 = 0.01;
L = [L1 L2 L3];
rc = [r1 r2 r3];
m = [m1 m2 m3 m4 m5 m6 m7];
MI = [Ixx2 Iyy2 Izz2 Ixx3 Iyy3 Izz3 Ixx4 Iyy4 Izz4];
g = 9.81;
%%%
coord = readmatrix("JointCoords3D_edit.xlsx");
time = coord(:,1);
gtc = 24:159; % gait cycle time is from 0.23 to 1.58 seconds
tfit = time(gtc);
dt = 0.01; % dt should be set to 0.01 only
tt = 0.23:dt:1.58;
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
GRFx = -GRdata_new(:,2);
GRFy = GRdata_new(:,3);
GRFz = GRdata_new(:,4);
xtrunk = 2-0.001*coord(gtc,2);
ytrunk = 0.001*coord(gtc,3);
ztrunk = 0.001*coord(gtc,4);
%%% left side
lxhip = 2-0.001*coord(gtc,5);
lyhip = 0.001*coord(gtc,6);
lzhip = 0.001*coord(gtc,7);
lxknee = 2-0.001*coord(gtc,8);
lyknee = 0.001*coord(gtc,9);
lzknee = 0.001*coord(gtc,10);
lxankleold = 2-0.001*coord(gtc,11);
lyankleold = 0.001*coord(gtc,12);
lzankleold = 0.001*coord(gtc,13);
lxpost = 2-0.001*coord(gtc,14);
lypost = 0.001*coord(gtc,15);
lzpost = 0.001*coord(gtc,16);
lxmeta = 2-0.001*coord(gtc,17);
lymeta = 0.001*coord(gtc,18);
lzmeta = 0.001*coord(gtc,19);
lxcop = 2-0.001*coord(gtc,38);
lycop = 0.001*coord(gtc,39)+0.02;
lzcop = 0.001*coord(gtc,40);
dlcop = sqrt((lxcop-lxpost).^2 + (lycop-lypost).^2 + (lzcop-lzpost).^2);
lxankle = lxpost;
lyankle = lypost;
lzankle = lzpost;
lxfoot = lxmeta;
lyfoot = lymeta;
lzfoot = lzmeta;
%%% right side
rxhip = 2-0.001*coord(gtc,20);
ryhip = 0.001*coord(gtc,21);
rzhip = 0.001*coord(gtc,22);
rxknee = 2-0.001*coord(gtc,23);
ryknee = 0.001*coord(gtc,24);
rzknee = 0.001*coord(gtc,25);
rxankleold = 2-0.001*coord(gtc,26);
ryankleold = 0.001*coord(gtc,27);
rzankleold = 0.001*coord(gtc,28);
rxpost = 2-0.001*coord(gtc,29);
rypost = 0.001*coord(gtc,30);
rzpost = 0.001*coord(gtc,31);
rxmeta = 2-0.001*coord(gtc,32);
rymeta = 0.001*coord(gtc,33);
rzmeta = 0.001*coord(gtc,34);
rxcop = 2-0.001*coord(gtc,35);
rycop = 0.001*coord(gtc,36)+0.02;
rzcop = 0.001*coord(gtc,37);
drcop = sqrt((rxcop-rxpost).^2 + (rycop-rypost).^2 + (rzcop-rzpost).^2);
rxankle = rxpost;
ryankle = rypost;
rzankle = rzpost;
rxfoot = rxmeta;
ryfoot = rymeta;
rzfoot = rzmeta;
%%% pelvis
xp = 0.5*(rxhip+lxhip);
yp = 0.5*(ryhip+lyhip);
zp = 0.5*(rzhip+lzhip);
%%%
Len_thigh = mean(sqrt((rxhip-rxknee).^2 + (ryhip-ryknee).^2 + (rzhip-rzknee).^2));
Len_shank = mean(sqrt((rxknee-rxankle).^2 + (ryknee-ryankle).^2 + (rzknee-rzankle).^2));
Len_foot = mean(sqrt((rxmeta-rxpost).^2 + (rymeta-rypost).^2 + (rzmeta-rzpost).^2));
a = mean(sqrt((xtrunk-xp).^2 + (ytrunk-yp).^2 + (ztrunk-zp).^2));
b = mean(sqrt((rxhip-xp).^2 + (ryhip-yp).^2 + (rzhip-zp).^2));
misc = [a b];
%%% Inverse kinematics
tht1 = atan2((ryhip-ytrunk)/a,(rxhip-xtrunk)/a);
%%% right side
c3 = -(rzknee-rzhip)/L1;
tht3 = acos(c3);
s12 = -(ryknee-ryhip)./(L1*sin(tht3));
c12 = -(rxknee-rxhip)./(L1*sin(tht3));
tht2 = atan2(s12,c12) - tht1;
s5=L2.^(-1).*((rzankle+(-1).*rzknee).*cos(tht1+tht2).^2.*cos(tht3)+( ...
  rxankle+(-1).*rxknee).*cos(tht1+tht2).*sin(tht3)+sin(tht1+tht2).*( ...
  (rzankle+(-1).*rzknee).*cos(tht3).*sin(tht1+tht2)+(ryankle+(-1).* ...
  ryknee).*sin(tht3)));
tht5 = asin(s5);
s4c5=L2.^(-1).*(((-1).*ryankle+ryknee).*cos(tht1+tht2)+(rxankle+(-1).* ...
  rxknee).*sin(tht1+tht2));
c4c5=L2.^(-1).*((rxankle+(-1).*rxknee).*cos(tht1+tht2).*cos(tht3)+(-1) ...
  .*(rzankle+(-1).*rzknee).*cos(tht1+tht2).^2.*sin(tht3)+sin(tht1+ ...
  tht2).*((ryankle+(-1).*ryknee).*cos(tht3)+((-1).*rzankle+rzknee).* ...
  sin(tht1+tht2).*sin(tht3)));
tht4 = atan2(s4c5./cos(tht5),c4c5./cos(tht5));
s7=L3.^(-1).*(ryankle.*cos(tht4).*sin(tht1).*sin(tht2)+(-1).*ryfoot.* ...
  cos(tht4).*sin(tht1).*sin(tht2)+rxankle.*cos(tht3).*sin(tht1).* ...
  sin(tht2).*sin(tht4)+(-1).*rxfoot.*cos(tht3).*sin(tht1).*sin(tht2) ...
  .*sin(tht4)+rzankle.*sin(tht3).*sin(tht4)+(-1).*rzfoot.*sin(tht3) ...
  .*sin(tht4)+cos(tht2).*sin(tht1).*((rxankle+(-1).*rxfoot).*cos( ...
  tht4)+((-1).*ryankle+ryfoot).*cos(tht3).*sin(tht4))+cos(tht1).*( ...
  cos(tht2).*(((-1).*ryankle+ryfoot).*cos(tht4)+((-1).*rxankle+ ...
  rxfoot).*cos(tht3).*sin(tht4))+sin(tht2).*((rxankle+(-1).*rxfoot) ...
  .*cos(tht4)+((-1).*ryankle+ryfoot).*cos(tht3).*sin(tht4))));
tht7 = asin(s7);
s56c7=L3.^(-1).*((-1).*(rzankle+(-1).*rzfoot).*cos(tht1+tht2).^2.*cos( ...
  tht3)+(-1).*(rxankle+(-1).*rxfoot).*cos(tht1+tht2).*sin(tht3)+sin( ...
  tht1+tht2).*((-1).*(rzankle+(-1).*rzfoot).*cos(tht3).*sin(tht1+ ...
  tht2)+((-1).*ryankle+ryfoot).*sin(tht3)));
c56c7=L3.^(-1).*(rxankle.*cos(tht3).*cos(tht4).*sin(tht1).*sin(tht2)+( ...
  -1).*rxfoot.*cos(tht3).*cos(tht4).*sin(tht1).*sin(tht2)+rzankle.* ...
  cos(tht4).*sin(tht3)+(-1).*rzfoot.*cos(tht4).*sin(tht3)+(-1).* ...
  ryankle.*sin(tht1).*sin(tht2).*sin(tht4)+ryfoot.*sin(tht1).*sin( ...
  tht2).*sin(tht4)+cos(tht2).*sin(tht1).*((-1).*(ryankle+(-1).* ...
  ryfoot).*cos(tht3).*cos(tht4)+((-1).*rxankle+rxfoot).*sin(tht4))+ ...
  cos(tht1).*(sin(tht2).*((-1).*(ryankle+(-1).*ryfoot).*cos(tht3).* ...
  cos(tht4)+((-1).*rxankle+rxfoot).*sin(tht4))+cos(tht2).*((-1).*( ...
  rxankle+(-1).*rxfoot).*cos(tht3).*cos(tht4)+(ryankle+(-1).*ryfoot) ...
  .*sin(tht4))));
tht56 = atan2(s56c7./cos(tht7),c56c7./cos(tht7));
for ii = 1:length(tht56)
    if tht56(ii) > 0
        tht56(ii) = tht56(ii) - 2*pi;
    end
end
tht6 = tht56 - tht5;
%%% left side
c9 = -(lzknee-lzhip)/L1;
tht9 = acos(c9);
s18 = -(lyknee-lyhip)./(L1*sin(tht9));
c18 = -(lxknee-lxhip)./(L1*sin(tht9));
tht8 = atan2(s18,c18) - tht1;
s11=L2.^(-1).*((lzankle+(-1).*lzknee).*cos(tht1+tht8).^2.*cos(tht9)+( ...
  lxankle+(-1).*lxknee).*cos(tht1+tht8).*sin(tht9)+sin(tht1+tht8).*( ...
  (lzankle+(-1).*lzknee).*cos(tht9).*sin(tht1+tht8)+(lyankle+(-1).* ...
  lyknee).*sin(tht9)));
tht11 = asin(s11);
s10c11=L2.^(-1).*(((-1).*lyankle+lyknee).*cos(tht1+tht8)+(lxankle+(-1).* ...
  lxknee).*sin(tht1+tht8));
c10c11=L2.^(-1).*((lxankle+(-1).*lxknee).*cos(tht1+tht8).*cos(tht9)+(-1) ...
  .*(lzankle+(-1).*lzknee).*cos(tht1+tht8).^2.*sin(tht9)+sin(tht1+ ...
  tht8).*((lyankle+(-1).*lyknee).*cos(tht9)+((-1).*lzankle+lzknee).* ...
  sin(tht1+tht8).*sin(tht9)));
tht10 = atan2(s10c11./cos(tht11),c10c11./cos(tht11));
s13=L3.^(-1).*(cos(tht8).*sin(tht1).*((lxankle+(-1).*lxfoot).*cos( ...
  tht10)+((-1).*lyankle+lyfoot).*cos(tht9).*sin(tht10))+lyankle.* ...
  cos(tht10).*sin(tht1).*sin(tht8)+(-1).*lyfoot.*cos(tht10).*sin( ...
  tht1).*sin(tht8)+lxankle.*cos(tht9).*sin(tht1).*sin(tht10).*sin( ...
  tht8)+(-1).*lxfoot.*cos(tht9).*sin(tht1).*sin(tht10).*sin(tht8)+ ...
  cos(tht1).*(cos(tht8).*(((-1).*lyankle+lyfoot).*cos(tht10)+((-1).* ...
  lxankle+lxfoot).*cos(tht9).*sin(tht10))+((lxankle+(-1).*lxfoot).* ...
  cos(tht10)+((-1).*lyankle+lyfoot).*cos(tht9).*sin(tht10)).*sin( ...
  tht8))+lzankle.*sin(tht10).*sin(tht9)+(-1).*lzfoot.*sin(tht10).* ...
  sin(tht9));
tht13 = asin(s13);
s1112c13=L3.^(-1).*((-1).*(lzankle+(-1).*lzfoot).*cos(tht1+tht8).^2.*cos( ...
  tht9)+(-1).*(lxankle+(-1).*lxfoot).*cos(tht1+tht8).*sin(tht9)+sin( ...
  tht1+tht8).*((-1).*(lzankle+(-1).*lzfoot).*cos(tht9).*sin(tht1+ ...
  tht8)+((-1).*lyankle+lyfoot).*sin(tht9)));
c1112c13=L3.^(-1).*(cos(tht8).*sin(tht1).*((-1).*(lyankle+(-1).*lyfoot).* ...
  cos(tht10).*cos(tht9)+((-1).*lxankle+lxfoot).*sin(tht10))+ ...
  lxankle.*cos(tht10).*cos(tht9).*sin(tht1).*sin(tht8)+(-1).* ...
  lxfoot.*cos(tht10).*cos(tht9).*sin(tht1).*sin(tht8)+(-1).* ...
  lyankle.*sin(tht1).*sin(tht10).*sin(tht8)+lyfoot.*sin(tht1).*sin( ...
  tht10).*sin(tht8)+cos(tht1).*(cos(tht8).*((-1).*(lxankle+(-1).* ...
  lxfoot).*cos(tht10).*cos(tht9)+(lyankle+(-1).*lyfoot).*sin(tht10)) ...
  +((-1).*(lyankle+(-1).*lyfoot).*cos(tht10).*cos(tht9)+((-1).* ...
  lxankle+lxfoot).*sin(tht10)).*sin(tht8))+lzankle.*cos(tht10).*sin( ...
  tht9)+(-1).*lzfoot.*cos(tht10).*sin(tht9));
tht1112 = atan2(s1112c13./cos(tht13),c1112c13./cos(tht13));
for ii = 1:length(tht1112)
    if tht1112(ii) > 0
        tht1112(ii) = tht1112(ii) - 2*pi;
    end
end
tht12 = tht1112 - tht11;
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
f8 = fit(tfit,tht8,'fourier8');
tht8 = f8(tt);
f9 = fit(tfit,tht9,'fourier8');
tht9 = f9(tt);
f10 = fit(tfit,tht10,'fourier8');
tht10 = f10(tt);
f11 = fit(tfit,tht11,'fourier8');
tht11 = f11(tt);
f12 = fit(tfit,tht12,'fourier8');
tht12 = f12(tt);
f13 = fit(tfit,tht13,'fourier8');
tht13 = f13(tt);
fxtrunk = fit(tfit,xtrunk,'fourier8');
xtrunk = fxtrunk(tt);
fytrunk = fit(tfit,ytrunk,'fourier8');
ytrunk = fytrunk(tt);
fztrunk = fit(tfit,ztrunk,'fourier8');
ztrunk = fztrunk(tt);
%%% acceleration of hip, and time derivatives of thetas
[vxtrunk,axtrunk] = differentiate(fxtrunk,tt);
[vytrunk,aytrunk] = differentiate(fytrunk,tt);
[vztrunk,aztrunk] = differentiate(fztrunk,tt);
[eps1,eta1] = differentiate(f1,tt);
[eps2,eta2] = differentiate(f2,tt);
[eps3,eta3] = differentiate(f3,tt);
[eps4,eta4] = differentiate(f4,tt);
[eps5,eta5] = differentiate(f5,tt);
[eps6,eta6] = differentiate(f6,tt);
[eps7,eta7] = differentiate(f7,tt);
[eps8,eta8] = differentiate(f8,tt);
[eps9,eta9] = differentiate(f9,tt);
[eps10,eta10] = differentiate(f10,tt);
[eps11,eta11] = differentiate(f11,tt);
[eps12,eta12] = differentiate(f12,tt);
[eps13,eta13] = differentiate(f13,tt);
%%% estimating GRFs
c = 1;
solGRF = zeros(length(tt),2); %% preallocating
%for k = 1:length(tt)
for k = 1:130
    q = [xtrunk(k) ytrunk(k) ztrunk(k) tht1(k) tht2(k) tht3(k) tht4(k) tht5(k) tht6(k) tht7(k) tht8(k) tht9(k) tht10(k) tht11(k) tht12(k) tht13(k)];
    qdot = [vxtrunk(k) vytrunk(k) vztrunk(k) eps1(k) eps2(k) eps3(k) eps4(k) eps5(k) eps6(k) eps7(k) eps8(k) eps9(k) eps10(k) eps11(k) eps12(k) eps13(k)];
    qddot = [axtrunk(k) aytrunk(k) aztrunk(k) eta1(k) eta2(k) eta3(k) eta4(k) eta5(k) eta6(k) eta7(k) eta8(k) eta9(k) eta10(k) eta11(k) eta12(k) eta13(k)];
    if (k>=1 && k<=53)
        solGRF(c,1) = 0;
        solGRF(c,2) = 0;
    elseif (k>=54 && k<=68)
        grf = double_stance_GRF(q,qdot,qddot);
        solGRF(c,1) = grf(1);
        solGRF(c,2) = grf(2);
    elseif (k>=69 && k<=121)
        grf = single_stance_GRF(q,qdot,qddot);
        solGRF(c,1) = grf(1);
        solGRF(c,2) = grf(2);
    elseif (k>=122 && k<=130)
        grf = double_stance_GRF(q,qdot,qddot);
        solGRF(c,1) = grf(1);
        solGRF(c,2) = grf(2);
    end
    c = c+1;
end
GRFx_new = solGRF(:,1);
GRFy_new = solGRF(:,2);
GRFz_new = zeros(length(tt),1);
%%%
fgrfx = fit(tfit,GRFx_new,'fourier8');
GRFx_new = fgrfx(tt);
GRFx_new(1:53) = 0;
fgrfy = fit(tfit,GRFy_new,'fourier8');
GRFy_new = fgrfy(tt);
GRFy_new(1:53) = 0;
%%%
figure;
plot(tt,GRFx_new,'b-','LineWidth',1);
grid on;
hold on;
plot(tt,GRFx,'r-','LineWidth',1);
xlabel('time (sec) \rightarrow');
ylabel('GRF_x (N) \rightarrow');
legend('3D 16 dof humanoid','Experimental');
%
figure;
plot(tt,GRFy_new,'b-','LineWidth',1);
grid on;
hold on;
plot(tt,GRFy,'r-','LineWidth',1);
xlabel('time (sec) \rightarrow');
ylabel('GRF_y (N) \rightarrow');
legend('3D 16 dof humanoid','Experimental');
%
rGRF = r3 - drcop;
%%% estimating joint torques
c = 1;
sola = zeros(length(tt),6); %% preallocating
solb = zeros(length(tt),6); %% preallocating
for k = 1:length(tt)
    q = [xtrunk(k) ytrunk(k) ztrunk(k) tht1(k) tht2(k) tht3(k) tht4(k) tht5(k) tht6(k) tht7(k) tht8(k) tht9(k) tht10(k) tht11(k) tht12(k) tht13(k)];
    qdot = [vxtrunk(k) vytrunk(k) vztrunk(k) eps1(k) eps2(k) eps3(k) eps4(k) eps5(k) eps6(k) eps7(k) eps8(k) eps9(k) eps10(k) eps11(k) eps12(k) eps13(k)];
    qddot = [axtrunk(k) aytrunk(k) aztrunk(k) eta1(k) eta2(k) eta3(k) eta4(k) eta5(k) eta6(k) eta7(k) eta8(k) eta9(k) eta10(k) eta11(k) eta12(k) eta13(k)];
    %GRF = [GRFx(k) GRFy(k) GRFz(k)];
    GRFa = [GRFx_new(k) GRFy_new(k) GRFz_new(k)];
    GRFb = [GRFxr(k) GRFyr(k) 0];
    r = rGRF(k);
%     sola(c,:) = torques3D_II(q,qdot,qddot,GRFa,r);
%     solb(c,:) = torques3D_II(q,qdot,qddot,GRFb,r);
    sola(c,:) = torquesHumanoid3D(q,qdot,qddot,GRFa,r);
    solb(c,:) = torquesHumanoid3D(q,qdot,qddot,GRFb,r);
    c = c+1;
end
T5a = sola(:,1);
T6a = sola(:,2);
T7a = sola(:,3);
T8a = sola(:,4);
T9a = sola(:,5);
T10a = sola(:,6);
T5b = solb(:,1);
T6b = solb(:,2);
T7b = solb(:,3);
T8b = solb(:,4);
T9b = solb(:,5);
T10b = solb(:,6);
%%% plotting
figure;
plot(tt,T5a,'b-','LineWidth',1);
grid on;
hold on;
plot(tt,T5b,'r-','LineWidth',1);
title('Hip flexion/extension moment');
legend('3D newton','2D lagrangian');
xlabel('time (sec) \rightarrow');
ylabel('joint moment (N-m) \rightarrow');
figure;
plot(tt,T6a,'b-','LineWidth',1);
grid on;
hold on;
plot(tt,T6b,'r-','LineWidth',1);
title('Hip abduction/adduction moment');
legend('3D newton','2D lagrangian');
xlabel('time (sec) \rightarrow');
ylabel('joint moment (N-m) \rightarrow');
figure;
plot(tt,T7a,'b-','LineWidth',1);
grid on;
hold on;
plot(tt,T7b,'r-','LineWidth',1);
title('Hip internal/external rotation moment');
legend('3D newton','2D lagrangian');
xlabel('time (sec) \rightarrow');
ylabel('joint moment (N-m) \rightarrow');
figure;
plot(tt,T8a,'b-','LineWidth',1);
grid on;
hold on;
plot(tt,T8b,'r-','LineWidth',1);
title('Knee flexion/extension moment');
legend('3D newton','2D lagrangian');
xlabel('time (sec) \rightarrow');
ylabel('joint moment (N-m) \rightarrow');
figure;
plot(tt,T9a,'b-','LineWidth',1);
grid on;
hold on;
plot(tt,T9b,'r-','LineWidth',1);
title('Ankle dorsiflexion/plantarflexion moment');
legend('3D newton','2D lagrangian');
xlabel('time (sec) \rightarrow');
ylabel('joint moment (N-m) \rightarrow');
figure;
plot(tt,T10a,'b-','LineWidth',1);
grid on;
hold on;
plot(tt,T10b,'r-','LineWidth',1);
title('Ankle inversion/eversion moment');
legend('3D newton','2D lagrangian');
xlabel('time (sec) \rightarrow');
ylabel('joint moment (N-m) \rightarrow');
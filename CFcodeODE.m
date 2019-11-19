%BMED 6310 Group Project Code 11-6
%Cystic Fibrosis Conditions
function CFcodeODE
CF = 1; %1 means the patient has CF, 0 means no CF
if CF == 1
    y0 = [1278.037, 64.7532, 15.04461, 12.75542, 1713.089, 100];
    % initial conditions: B = y0(1), P = y0(2), A = y0(3), M = y0(4), H = y0(5),
    % D = y0(6), LF = y0(7), and PH = y0(8)
elseif CF == 0
    y0 = [1.383617,6.091615, 2.393149, 1.075577, 2431.41, 100];
end
%AB = 0; %No antibiotic treatment so AB not needed

t0 = 0;
tf = 400;
hr = 1; %hr = step size
t = t0:hr:19;
[~,yode] = ode45(@diffeqs,t,y0);

t2 = 20:hr:159; %reset bacteria (B) at t = 20, 160, and 300
y02 = [3000; yode(end,2); yode(end,3); yode(end,4); yode(end,5); yode(end,6)];
[~,yode2] = ode45(@diffeqs,t2,y02);

t3 = 160:hr:299;
y03 = [5000; yode2(end,2); yode2(end,3); yode2(end,4); yode2(end,5); yode2(end,6)];
[~,yode3] = ode45(@diffeqs,t3,y03);

t4 = 300:hr:tf;
y04 = [10000; yode3(end,2); yode3(end,3); yode3(end,4); yode3(end,5); yode3(end,6)];
[~,yode4] = ode45(@diffeqs,t4,y04);

if CF == 1
    H0 = 1713.089; %initial H and D (doesn't change throughout)
    D0 = 100;
elseif CF == 0
    H0 = 2431.41;
    D0 = 100;
end
H = [yode(:,5); yode2(:,5); yode3(:,5); yode4(:,5)];
D = [yode(:,6); yode2(:,6); yode3(:,6); yode4(:,6)];
LF = 100.*((H./D)./(H0./D0));
PH = 100.*H./H0;

tt = t0:hr:tf;
%B = y(1); P = y(2); A = y(3); M = y(4); H = y(5); D = y(6); LF = y(7); PH = y(8);
figure;
subplot(2,2,1);
plot(tt, LF, 'g') %LF
hold on;
plot(tt, PH, 'r') %PH
plot(tt, [yode(:,6); yode2(:,6); yode3(:,6); yode4(:,6)], 'b') %D
hold off;
xlabel('time')
legend('LF','PH','D');
ylim([60 120]);
subplot(2,2,2);
plot(tt, [yode(:,1); yode2(:,1); yode3(:,1); yode4(:,1)], 'r') %B
xlabel('time')
legend('B');
ylim([0 10000]);
subplot(2,2,3);
plot(tt, [yode(:,3); yode2(:,3); yode3(:,3); yode4(:,3)],'color',[0 .5 0]) %A
hold on;
plot(tt, [yode(:,4); yode2(:,4); yode3(:,4); yode4(:,4)],'color',[1 .5 0]) %M
hold off;
xlabel('time')
legend('A','M');
ylim([0 30]);
subplot(2,2,4);
plot(tt, [yode(:,2); yode2(:,2); yode3(:,2); yode4(:,2)],'color',[0 0 .7]) %P
xlabel('time')
legend('P');
ylim([60 140]);
end

function f = diffeqs(~,y)
CF = 1; %change whether CF is present

B = y(1); P = y(2); A = y(3); M = y(4); H = y(5); D = y(6);
Bp = (((B^0.8)*(M^1.2)) - ((0.8*B^1.2)*(P^0.1)))*(M^-2.4);
Pp = (20*(B^0.2)*(D^0.4)*(A^-0.1))-50*(P^0.5);
Ap = (15*(B^0.1)*(P^0.1))-12*(A^0.5);
Mp = (0.16*(B^0.25)*(2^CF))-0.15*(M^(2-CF));
Hp = 500-40*(H^0.3)*(P^0.2)*(A^(-0.2));
Dp = (40*(H^0.3)*(P^0.2)*(A^-0.2))-5*D;
f = [Bp;Pp;Ap;Mp;Hp;Dp];
end

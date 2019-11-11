%BMED 6310 Group Project Code 11-6

%Cystic Fibrosis Conditions
CF = 1;
B = 1278.037;
P = 64.7532;
A = 15.04461;
H0 = 1713.089;
D0 = 100;
M = 12.75542;
 
H = H0;
D = D0;
LF = 100*((H/D)/(H0/D0));
PH = 100*H/H0;

AB = 0;

t0 = 0;
tf = 400;
tvec = [t0:1:tf];
 
for i = 1:10
    i
    Bp = (B^0.8 * M^1.2 - 0.8 * B^1.2 * P^.1 * (AB+1)) * M^-2.4
    Pp = 20 * B^0.2 * D^0.4 * A^-0.1 - 50 * P^0.5
    Ap = 15 * B^0.1 * P^0.1 - 12 * A^0.5
    Hp = 500 - 40 * H^0.3 * P^0.2 * A^0.5
    Dp = 40 * H^0.3 * P^0.2 * A^-0.2 - 5 * D
    Mp = 0.16 * B^0.25 * 2^CF - 0.15 * M^(2-1*CF)
    LFp = 100*((Hp/Dp)/(H0/D0))
    PHp = 100*Hp/H0
    if i == 20
        B = 3000;
    elseif i == 22
        AB = 5;
    elseif i == 32
        AB = 0;
    elseif i == 160
        B = 5000;
    elseif i == 162
        AB = 5;
    elseif i == 172
        AB = 0;
    elseif i == 300
        B = 10000;
    elseif i == 302
        AB = 5;
    elseif i == 312
        AB = 0;
    end
    plot(tvec, LFp, 'g', tvec, PHp, 'r', tvec, Dp, 'b')
    hold on
end 


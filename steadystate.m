function F = steadystate(y)
CF=0.1;
% where x is B,P,A,M,H,D
% B = y(1); P = y(2); A = y(3); M = y(4); H = y(5); D = y(6);
F(1) = (((y(1)^0.8)*(y(4)^1.2)) - ((0.8*y(1)^1.2)*(y(2)^0.1)))*(y(4)^-2.4);
F(2) = (20*(y(1)^0.2)*(y(6)^0.4)*(y(3)^-0.1))-50*(y(2)^CF);
F(3) = (15*(y(1)^0.1)*(y(2)^0.1))-12*(y(3)^CF);
F(4) = (0.16*(y(1)^0.25)*(2^CF))-0.15*(y(4)^(2-CF));
F(5) = 500-40*(y(5)^0.3)*(y(2)^0.2)*(y(3)^(-0.2));
F(6) = (40*(y(5)^0.3)*(y(2)^0.2)*(y(3)^-0.2))-5*y(6);

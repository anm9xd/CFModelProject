% Solve Steady state values for CF= 0.5
fun=@steadystate;
y0 = [1000,25,10,10,1000,100];
[y,fval] = fsolve(fun,y0)

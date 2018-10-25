function C = example2()

syms t;

C = (0.25 + 0.01*sin(3*t) + 0.02*cos(5*t) + 0.01*cos(6*t) + 0.01*cos(8*t) + 0.01*cos(10*t)).*exp(1i*t);



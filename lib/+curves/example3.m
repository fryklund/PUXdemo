function C = example3()

syms t;

C = (1+0.2*(cos(-5*t) + sin(-1*t))).*exp(1i*t) + 0.17*1i;



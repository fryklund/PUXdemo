function C = example3Exterior()

syms t;

C = 0.3*(1+0.1*(sin(-6*t) + sin(-3*t))).*exp(-1i*t)+ 0.17*1i;
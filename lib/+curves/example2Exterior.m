function C = example2Exterior()
syms t;
trans = pi/4;
C = (0.05 + 0.005*cos(-2*(t+trans)) + 0.005*sin(-3*(t+trans)) + 0.005*cos(-5*(t+trans)) + 0.005*cos(-7*(t+trans))).*exp(-1i*(t+trans));
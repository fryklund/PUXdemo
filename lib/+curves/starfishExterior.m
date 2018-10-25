function C = starfish_exterior(amplitude, n_arms)
% C = starfish(amplitude, n_arms)

syms t;
C = 0.4*(1 + amplitude*cos(-n_arms*t)).*exp(-1i*t);
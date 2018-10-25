function C = circle(center, radius)
% C = starfish(amplitude, n_arms)

syms t;
C = center + radius * exp(1i*t);


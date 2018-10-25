function C = circle_exterior(center, radius)
% C = starfish(amplitude, n_arms)

syms t;
C = center + radius * exp(-1i*t);


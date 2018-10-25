function curve = diffCurve(parametrisation)
% INTPUTS:
% Parametrisation of curve s in paramater t, such that s = s(t) for t in
% [0,2*pi).

% OUTPUTS:
%       curves.tau:
%           The parametrization of the curve
%       curves.dtau:
%           The derivative, with respect to t, of the parametrization
%       curves.d2tau:
%           The second derivative, with respect to t, of the parametrization
%       curves.normal:
%           The inward normal at point curves.tau(t) to the curve given by tau.
%       curve.tangent:
%           The tangent at point curves.tau(t) to the curve given by tau.
%       curve.kappa: 
%           The curvature at point curves.tau(t) to the curve given by tau.


d = diff(parametrisation);
d2 = diff(d);

curve.sym = parametrisation;
curve.tau = matlabFunction(parametrisation);
curve.dtau = matlabFunction(d);
curve.d2tau = matlabFunction(d2);
curve.normal = matlabFunction(-1i*d / abs(d));
curve.tangent = matlabFunction(d / abs(d));
curve.kappa = matlabFunction(-imag(d*conj(d2))/abs(d)^3);

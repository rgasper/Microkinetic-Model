function dydt = odes(t,y,flags,params)
 % params = [alpha k1p km1 k2 km2p];
 alpha = params(1);
 k1p = params(2);
 km1 = params(3);
 k2 = params(4);
 km2 = params(5);
 k3 = params(6);
 km3p = params(7);
 
 % y(1) = psi-A, y(2) = psi-C, y(3) = theta-A, y(4) = theta-B
 dy1dt = alpha*( -k1p*y(1)*(1-y(3)-y(4)) + km1*y(3) );
 dy2dt = alpha*( k3*y(4) - km3p*y(2)*(1-y(3)-y(4)) );
 dy3dt = k1p*y(1)*(1-y(3)-y(4)) - km1*y(3) - k2*y(3) + km2*y(4);
 dy4dt = k2*y(3) - km2*y(4) - k3*y(4) + km3p*y(2)*(1-y(3)-y(4));
 
 dydt = [dy1dt ; dy2dt ; dy3dt; dy4dt]; % must return as column vector

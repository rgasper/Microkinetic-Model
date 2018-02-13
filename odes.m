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
 r1 = k1p*y(1)*(1-y(3)-y(4));
 rm1 = km1*y(3);
 r2 = k2*y(3);
 rm2 = km2*y(4);
 r3 = k3*y(4);
 rm3 = km3p*y(2)*(1-y(3)-y(4));
 

 dy1dt = alpha*( -r1 + rm1 );
 dy2dt = alpha*(  r3 - rm3 );
 dy3dt = r1 - rm1 - r2 + rm2;
 dy4dt = r2 - rm2 - r3 + rm3;
 
 dydt = [dy1dt ; dy2dt ; dy3dt; dy4dt]; % must return as column vector

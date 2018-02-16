function dydt = odes(t,y,flags,params)
% params = [alpha k1p km1 k2p km2p k3 km3p k4p km4p];
 alpha = params(1);
 k1p = params(2);
 km1 = params(3);
 k2p = params(4);
 km2p = params(5);
 k3 = params(6);
 km3p = params(7);
 k4p = params(8);
 km4p = params(9);
 
%  y(1) = psi-A, y(2) = psi-B, y(3) =psi-C2
%  y(4) = theta-A, y(5) = theta-B, y(6) = theta-C
 theta_s = 1-y(4)-y(5)-y(6);
 
 %defining reaction rates
 r1 = k1p*y(1)*(theta_s);
 rm1 = km1*y(4);
 r2 = k2p*y(4)*theta_s^2;
 rm2 = km2p*y(5)*y(6)^2;
 r3 = k3*y(5);
 rm3 = km3p*y(2)*(theta_s);
 r4 = k4p*y(6)^2;
 rm4 = km4p*y(3)*theta_s^2;
 
 %odes
 dy1dt = alpha*( -r1 + rm1 );
 dy2dt = alpha*(  r3 - rm3 );
 dy3dt = alpha*(  r4 - rm4 );
 dy4dt = r1 - rm1 - r2 + rm2;
 dy5dt = r2 - rm2 - r3 + rm3;
 dy6dt = r2 - rm2 - r4 + rm4;

 % must return as column vector
 dydt = [dy1dt ; dy2dt ; dy3dt; dy4dt; dy5dt; dy6dt];

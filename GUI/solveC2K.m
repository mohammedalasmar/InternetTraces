function y=solveC2K(C2)
global ave VAR epsi T kk

% epsi=10^(-2) ,VAR=(4.9e6)^2 , T=1 , ave=18.9e6;
rr=(C2*T-ave(kk)*T).^2;
y=rr+VAR(kk)*log(2*pi*rr/VAR(kk))+2*log(epsi(kk))*VAR(kk) ;

 
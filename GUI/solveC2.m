function y=solveC2(C2)
global ave VAR epsi T  

% epsi=10^(-2) ,VAR=(4.9e6)^2 , T=1 , ave=18.9e6;
rr=(C2*T-ave *T).^2;
y=rr+VAR *log(2*pi*rr/VAR )+2*log(epsi )*VAR  ;

 
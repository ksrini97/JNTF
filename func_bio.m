function [f]=func_bio(x,Z_ftir,Z_hnmr,P_ftir,P_hnmr)
global F;global H;global FH; 
global alpha; global beta; global gamma; global lambda;global R_f;global R_h;


p=3*R_f;q=p+9*R_f;r=q+1765*R_f;s=r+2048*R_f;
A=reshape(x(1:p),[3,R_f]);
B=reshape(x(p+1:q),[9,R_f]);
H1=reshape(x(q+1:r),[1765,R_f]);
H2=reshape(x(r+1:s),[2048,R_f]);

f=norm(P_ftir.*(Z_ftir-ktensor({A,B,H1})))+norm(P_hnmr.*(Z_hnmr-ktensor({A,B,H2})))+alpha*norm(H1'*FH*H2)+beta*norm(H1*F*H1')+beta*norm(H2*H*H2')+gamma*norm(ktensor({A,B}))+lambda*(norm(H1)+norm(H2));

end

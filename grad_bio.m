function [g]=grad_bio(x,Z_ftir,Z_hnmr,P_ftir,P_hnmr)
 
global F;global H;global FH; 
global alpha; global beta; global gamma; global lambda;global R_f; global R_h;


p=3*R_f;q=p+9*R_f;r=q+1765*R_f;s=r+2048*R_f;
A=reshape(x(1:p),[3,R_f]);
B=reshape(x(p+1:q),[9,R_f]);
H1=reshape(x(q+1:r),[1765,R_f]);
H2=reshape(x(r+1:s),[2048,R_f]);



 gradA=@(A,B,H1,H2) double(2*(tenmat((P_ftir.*ktensor({A,B,H1})),1)-tenmat((P_ftir.*Z_ftir),1))*khatrirao(H1,B))+double(2*(tenmat((P_hnmr.*ktensor({A,B,H2})),1)-tenmat((P_hnmr.*Z_hnmr),1))*khatrirao(H2,B))+gamma*double(ktensor({A,B}))*B;
 gradB=@(A,B,H1,H2) double(2*(tenmat((P_ftir.*ktensor({A,B,H1})),2)-tenmat((P_ftir.*Z_ftir),2))*khatrirao(H1,A))+double(2*(tenmat((P_hnmr.*ktensor({A,B,H2})),2)-tenmat((P_hnmr.*Z_hnmr),2))*khatrirao(H2,A))+gamma*double((ktensor({A,B})))'*A;
 gradH1=@(A,B,H1,H2)double(2*(tenmat((P_ftir.*ktensor({A,B,H1})),3)-tenmat((P_ftir.*Z_ftir),3))*khatrirao(B,A))+2*alpha*((H1'*FH*H2)*(FH*H2)')'+2*beta*(H1*F')+2*lambda*H1;
 gradH2=@(A,B,H1,H2)double(2*(tenmat((P_hnmr.*ktensor({A,B,H2})),3)-tenmat((P_hnmr.*Z_hnmr),3))*khatrirao(B,A))+2*alpha*((H1'*FH*H2)*(H1'*FH))'+2*beta*H2*H'+2*lambda*H2;


grA=gradA(A,B,H1,H2);
grB=gradB(A,B,H1,H2);
grH1=gradH1(A,B,H1,H2);
grH2=gradH2(A,B,H1,H2);
g = [grA(:);grB(:);grH1(:);grH2(:)];
end

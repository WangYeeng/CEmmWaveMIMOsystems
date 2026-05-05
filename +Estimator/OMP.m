function[xhat]=OMP(y,A,t)
%y:measurement
%A:sensing matrix
%t:number of iteration
%k:sparse level

[M,N]=size(A);
xhat=zeros(N,1);
At=zeros(M,t);
pos_xhat=zeros(1,t);
r=y;  %Initialization of residual
for i=1:t
	product=A'*r;  %Inner product of each column in sensing matrix and residual
	[~,pos]=max(abs(product)); %Find the column most relevant to residual
    pos=pos(end); 
	At(:,i)=A(:,pos); 
	pos_xhat(i)=pos; %save the column index
	A(:,pos)=zeros(M,1); 
	xhat_ls=(At(:,1:i)'*At(:,1:i))^(-1)*At(:,1:i)'*y;  %Least Square 
	r=y-At(:,1:i)*xhat_ls;   %update residual
end
xhat(pos_xhat)=xhat_ls;


function[hhat]=MiAMP(y,A,alf)

% y: measurement
% A: sensing matrix

[M,N]=size(A);
S=size(y,2);  
eta_l1=@(r,lam) exp(1i*angle(r)).*max(bsxfun(@minus,abs(r),lam),0); %soft-thresholding shrinkage function
Bmf=A'; % matched filter
hhat=zeros(N,S); % Initialization of signal estimation
v=y; % Initialization of residual
i = 0;
while i < 40
    rhat=hhat+Bmf*v;
    rvar=sum(abs(v).^2,1)/M;
    lam=alf*sqrt(rvar);
    hhat=eta_l1(rhat,lam);
	g=(1/M)*sum((hhat~=0).*(1-lam./(2*abs(rhat)+eps)),1);
    c=(1/M)*sum((hhat~=0).*(0.5*lam.*sqrt(rhat).*(1./((conj(rhat)).^(3/2)+eps))),1);
	v=y-A*hhat+bsxfun(@times,v,g)+bsxfun(@times,conj(v),c);
    i = i+1;
end
end

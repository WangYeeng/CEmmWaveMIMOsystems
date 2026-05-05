function X = SSD(Y,Phi,L,N,M,fs,fc,I)
%%% initialization
R = Y;
set = (-(N-1)/2:1:(N-1)/2)/N;
%%% JOMP
for i = 1 : L
    P = sum(abs(Phi'* R).^2,2);
    for n = 1 : N
        theta1 = set(n);
        delta1 = abs(round(theta1 * fs/(M*fc)*((M-1)/2) * N));
        power_window = mod((n - delta1 : n + delta1) - 1,N) + 1;
        P1(n) = sum(P(power_window))/(2*delta1 + 1);
    end
    [value(i),order] = max(P1);
    theta = set(order);
    for m = 1 : M
        delta = round(theta * fs/(M*fc)*(m-1-(M-1)/2) * N);
        path_support(:,i,m) = mod((order + delta - 4 : order + delta + 4) - 1,N) + 1;
        Phi_s = Phi(:,path_support(:,i,m));
        Temp = zeros(N,1);
        Temp(path_support(:,i,m),:) = inv(Phi_s'*Phi_s)*Phi_s'*R(:,m);
        
        for ii = 1 : I
            [~,order1] = max(abs(Temp));
            path_support(:,i,m) = mod((order1 - 4 : order1 + 4) - 1,N) + 1;
            Phi_s = Phi(:,path_support(:,i,m));
            Temp = zeros(N,1);
            Temp(path_support(:,i,m),:) = inv(Phi_s'*Phi_s)*Phi_s'*R(:,m);
        end

        R(:,m) = R(:,m) - Phi*Temp;
    end
end
for m = 1 : M
    support = unique(reshape(path_support(:,:,m),[],1));
%     X(:,m) = WOMP(Y(:,m),Phi,support,15);
    Phi_final = Phi(:,support);
    X(:,m) = zeros(N,1);
    X(support,m) = inv(Phi_final'*Phi_final)*Phi_final'*Y(:,m);
end
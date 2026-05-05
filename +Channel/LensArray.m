function [UN] = LensArray(N_ver,N_hor,d,type)

M = N_ver*N_hor;
switch type
      case 'ULA'
         UN=(1/sqrt(M))*exp(1i*2*pi*(-(M-1):2:M-1)'/2*d*(-(M-1):2:M-1)*(1/M));
      case 'UPA'
         UN1=(1/sqrt(N_ver))*exp(1i*2*pi*(-(N_ver-1):2:N_ver-1)'/2*d*(-(N_ver-1):2:N_ver-1)*(1/N_ver));
         UN2=(1/sqrt(N_hor))*exp(1i*2*pi*(-(N_hor-1):2:N_hor-1)'/2*d*(-(N_hor-1):2:N_hor-1)*(1/N_hor));
         UN=kron(UN1,UN2);
end
end


function [Fk,Gk] = theodorsen(k)

%T = zeros(1,2);
C = 1-0.165/(1-1i*0.0455/k)-0.335/(1-1i*0.33/k);
Fk = real(C);
Gk = imag(C);
return
end
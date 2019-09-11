function [x,C,amp] = model_varpro(y,fre, damp, fun_type,t)

K = length(fre);
C = zeros(length(t),K);
count = 1;
for j = 1:K
    switch fun_type(j)
        case 1
            d_l = damp(count);
            C(:,j) = exp((-d_l + 1i*2*pi*fre(j)).*t);
            count = count + 1;
        case 2
%             C(:,j) = exp((-damp(count)*t + 1i*2*pi*fre(j)).*t);
            d_g = damp(count)^2;
            C(:,j) = exp((-d_g*t + 1i*2*pi*fre(j)).*t);
            count = count+1;
        case 3
            d_l = damp(count);
            count = count+1;
%             d_g = damp(count);
            d_g = damp(count)^2;
            count = count+1;
            C(:,j) = ((alpha(j)*exp(-d_l*t) + (1-alpha(j))*exp(-d_g*t.*t)).*exp(1i*2*pi*fre(j).*t));
        case 4
            d_l = damp(count);
            count = count+1;
%             d_g = damp(count);
            d_g = damp(count)^2;
            count = count+1;
            C(:,j) = exp((-d_l -d_g*t + 1i*2*pi*fre(j)).*t);                
    end
end
C1 = [real(C);imag(C)];%abs(C);%[real(C);imag(C)];
ap = C1\[real(y);imag(y)];%C1\[abs(y)];%C1\[real(y);imag(y)];
% ll = find(ap<0);
% ap(ll) = 0;
x = C*ap;
amp = ap';

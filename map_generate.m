function [qq]=map_generate(flip, FOV, a0,res)

%psi = -ax*(pi/2)^2/(a0+4*ax);
theta = flip/(pi/2)^2;
k0 = 1/FOV/2;
%theta*a0 = -(5+theta)*ax;
ax=-theta*a0/(5+theta);


ay =ax;

for j=1:res
    for k=1:res
        qq(j,k) = a0 +2*ax*cos(2*pi*k0*(j-5)) + 2*ay*cos(2*pi*k0*(k-5));
        q(j,k) = a0 +2*ax*(2*pi*(1/FOV/2))^2*(j)^2; 
                +2*ay - ay*(2*pi*(1/FOV/2))^2*(k)^2;
    end
end
%mesh(q)
% figure;

posx = [-res./2:((res./2)-1)]; 
posy = [-res./2:((res./2)-1)]; 
% mesh(posx,posy,qq);
B0_map_spoke=qq;
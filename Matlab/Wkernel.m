function [ W ] = Wkernel(r,h,n)
%WDEF Summary of this function goes here
%   Detailed explanation goes here

W = 0;

if(n == 1)
    %Wdef
    W = (315/(64*pi*h^9))*(h^2 -(r*r'))^3;
end

if(n == 2)
    %gradWpresure 
    W = -(45/(pi*h^6))*(r/(sqrt(r*r')))*(h-sqrt((r*r')))^2;
end

if(n == 3)
    %laplacianWviscosity 
    W = (45/(pi*h^6))*(h-(sqrt((r*r'))));  
end

if(n == 4)
    %gradWdef 
    W = -(945/(32*pi*h^9))*r*(h^2-((r*r')))^2;
end

if(n == 5)
    %laplacianWdef 
    W = -(945/(32*pi*h^9))*((h^2)-(r*r'))*(3*h^2 - 7*(r*r'));
end

end


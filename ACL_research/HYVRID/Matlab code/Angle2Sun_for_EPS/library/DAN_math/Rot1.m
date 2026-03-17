function [ y ] = Rot1( x )
%#codegen
y=[1 0 0;
   0 cos(x) sin(x);
   0 -sin(x) cos(x)];
end
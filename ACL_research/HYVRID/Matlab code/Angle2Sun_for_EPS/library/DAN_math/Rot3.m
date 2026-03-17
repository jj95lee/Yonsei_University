function [ y ] = Rot3( x )
%#codegen
y=[cos(x) sin(x) 0;
   -sin(x) cos(x) 0;
   0 0 1];
end
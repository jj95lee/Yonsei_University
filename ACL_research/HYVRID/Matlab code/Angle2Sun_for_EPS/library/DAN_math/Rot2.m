function [ y ] = Rot2( x )
%#codegen
y=[cos(x)  0  -sin(x);
   0       1  0;
   sin(x) 0  cos(x)];
end
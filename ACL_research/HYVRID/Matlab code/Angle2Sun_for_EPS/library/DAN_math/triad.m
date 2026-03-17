function R   = triad( u1,v1,u2,v2 )
%#codegen
% Transformation from vector to direction cosine matrix
% u1, v1 are the vectors before transformation
% u2, v2 are the vectors after transformation

q1 = u1;
r1 = cross(u1,v1)/norm(cross(u1,v1));
s1 = cross(q1,r1);

q2 = u2;
r2 = cross(u2,v2)/norm(cross(u2,v2));
s2 = cross(q2,r2);

M1=[q1 r1 s1];
M2=[q2 r2 s2];

R = M2*M1';

end


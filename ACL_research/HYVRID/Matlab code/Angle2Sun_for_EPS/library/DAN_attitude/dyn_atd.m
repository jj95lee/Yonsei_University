function dX = dyn_atd(t, X, u, MOI)

w		= [X(1); X(2); X(3)];
q		= [X(4); X(5); X(6); X(7)];
u		= [u(1); u(2); u(3)];

dw		= -( MOI^-1 )*cross( w, MOI*w ) + ( MOI^-1 )*u;

dq		= [	-0.5*q(2:4,1)'*w;
			0.5*( q(1)*w + cross(q(2:4,1), w) )];

dX		= [	dw;
			dq];
		
end
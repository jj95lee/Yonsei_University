function dX = dyn_atd_rwa(t, X, options, u, MOI, MOI_w)

q		= [X(1); X(2); X(3); X(4)];
w		= [X(5); X(6); X(7)];
w_w		= [X(8); X(9); X(10)];
u_t		= [u(1); u(2); u(3)];
u_w		= [u(4); u(5); u(6)];

dw		= -( MOI^-1 )*cross( w, (MOI*w + MOI_w*w_w) ) + ( MOI^-1 )*u_t;

dq		= [	-0.5*q(2:4,1)'*w;
			0.5*( q(1)*w + cross(q(2:4,1), w) )];

dw_w	= MOI_w\u_w;
		
dX		= [	dq;
			dw;
			dw_w];
		
end
function xout = sat(xin, lim)
[a, b] = size(xin);
xout = zeros(a, b);

if a == 1
	for i = 1:b
		if abs(xin(1, i)) > lim
			xout(1, i) = sign(xin(1, i))*lim;
		else
			xout(1, i) = xin(1, i);
		end
	end
	
elseif b == 1
	for i = 1:a
		if abs(xin(i, 1)) > lim
			xout(i, 1) = sign(xin(i, 1))*lim;
		else
			xout(i, 1) = xin(i, 1);
		end
	end
	
end
function[t y]= rk4c(ODE,ti,tf,h,yi) 

t=ti:h:tf;
n=length(t);
l=length(yi);


y=zeros(l,n);
y(:,1)=yi;


for i=1:n-1
    
    tm=t(i)+h/2;
    
    
    K1=h.*feval(ODE,t(i),y(:,i));
    y2=y(:,i)+(K1./2);
    
    K2=h.*feval(ODE,tm,y2);
    y3=y(:,i)+(K2./2);
    
    K3=h.*feval(ODE,tm,y3);
    
    y4=y(:,i)+K3;
    K4=h.*feval(ODE,t(i)*h,y4);
    
    y(:,i+1)=y(:,i)+1/6*(K1+2*K2+2*K3+K4);
end
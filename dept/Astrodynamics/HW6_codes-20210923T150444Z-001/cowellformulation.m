function Xdot = cowellformulation(t,X,opt)

opt = 3;

mu = 398600.4418;
Re = 6378.137;
we = [0;0;7.2921*10^-5];

r=sqrt(X(1)^2+X(2)^2+X(3)^2);

%% Case
if opt==1       % Two body case
    a = zeros(3,1);
elseif opt==2   % J2 case
    %% Zonal harmonics
    j2 = 1082.645*10^-6;
    
    a_j2= [-3*j2*mu*Re^2*X(1)/(2*r^5) * (1 - 5*X(3)^2/r^2) ;
        -3*j2*mu*Re^2*X(2)/(2*r^5) * (1 - 5*X(3)^2/r^2) ;
        -3*j2*mu*Re^2*X(3)/(2*r^5) * (3 - 5*X(3)^2/r^2) ];
    %%
    a = a_j2;
elseif opt==3   % J3, sectorial, tesseral, air drag, solar and moon gravity, solar radiation
    %% Zonal harmonics
    j3 = -2.546*10^-6;
    
    a_j3=[-5*j3*mu*Re^3*X(1)/(2*r^7) * (3*X(3) - 7*X(3)^3/r^2) ;
        -5*j3*mu*Re^3*X(2)/(2*r^7) * (3*X(3) - 7*X(3)^3/r^2) ;
        -5*j3*mu*Re^3/(2*r^7) * (6*X(3)^2 - 7*X(3)^4/r^2 - 3*r^2/5) ];
    
    %% Sectorial (l=2, m=2)
    %a_s = 0
    
    %% Tesseral (l=3, m=1)
    %a_t = 0
    
    %% Air drag
    ApM = 0.01; % 0.01 m^2/kg
    Cd = 2.2;  % Drag coeff
    den = 1.585*10^-12*exp(-(r-Re-1000)/60.828);   % Exponential model at 450-500km
    v_rel = X(4:6) - cross(we,X(1:3));
    a_air = -1/2*Cd*ApM*den*norm(v_rel)*v_rel;

    %% 3rd body (Sun)
    mu_sun = 132712440018;
    jd = jday(2014, 12, 2, 13, 20, 10);
    r_sun = (sun(jd))'*149597870.700;
    a_sun = -mu_sun/norm(r_sun)^3*(X(1:3) - 3*(X(1:3))'*r_sun/norm(r_sun)^2*r_sun);

    %% 3rd body (Moon)
    mu_moon = 4902.8;
    r_moon = (moon(jd))'*Re;
    a_moon = -mu_moon/norm(r_moon)^3*(X(1:3) - 3*(X(1:3))'*r_moon/norm(r_moon)^2*r_moon - 15/2*((X(1:3))'*r_moon/norm(r_moon)^2)^2*r_moon);

    %% Solar radiation
    psr = 4.57*10^-6;   % kg*km/sec^2*m^2
    Cr = 1.5;    % reflectivity
    a_sr = psr*Cr*ApM*(X(1:3)-r_sun)/norm(X(1:3)-r_sun);

    %%
    a = a_j3 + a_air + a_sun + a_moon + a_sr;
    
end

%% EOM
Xdot = [X(4);
        X(5);
        X(6);
        -mu/r^3*X(1) + a(1);
        -mu/r^3*X(2) + a(2);
        -mu/r^3*X(3) + a(3)];
    
    
end


constFuselage = calcFuselageDerivedInputs();
for i=1:3
S_wet_body(i) =2*pi*constFuselage.radius*constFuselage.length;
cf_body(i) = .455/((log10(constFuselage.Re))^2.58*(1+(.144*constFuselage.Mach)^2)^.65);
FF_body(i) = .9+5/(constFuselage.f)^1.5+constFuselage.f/400;
cf_wing(i) = 0;
FF_wing(i) = (1+(.6/constFuselage.Xc)*(constFuselage.t/constFuselage.c)+100*(constFuselage.t/constFuselage.c)^4)*((1.35*constFuselage.Mach)^.18*(cos(constFuselage.lambda))^.28);
S_wetwing(i) = 0;
Dq(i) = (.139+.419*(constFuselage.Mach-.161)^2)*constFuselage.A_base;
Cd0_misc(i) = Dq/constFuselage.S_ref;
Cd0_tot(i) = ((cf_body*FF_body*S_wet_body)+(FF_wing*constFuselage.Q_wing*cf_wing*S_wetwing))/constFuselage.S_ref+Cd0_misc;
end
function constFuselage = calcFuselageDerivedInputs()% this function calculates all the derived inputs needed for the fuselage
%% Defining Given Constants
% Standard Atmosphere (assume sea level)
constFuselage.alt = 0; % altitude
[T, constFuselage.a, P, constFuselage.rho] = atmoscoesa(constFuselage.alt) % temp, speed of sound, pressure, density


constFuselage.vel = [,,]; % velocity
constFuselage.length = 0; % length
constFuselage.mu = 0; % viscosity
constFuselage.k = 0; % surface whatever
constFuselage.Amax = 0; % max area
constFuselage.radius=0;
constFuselage.Q_body=1;
constFuselage.X = 0;
constFuselage.c = 0;
constFuselage.Xc=constFuselage.X/constFuselage.c;
constFuselage.t = 0;
constFuselage.lambda=  0;%sweep angle
constFuselage.Q_wing = 0;
constFuselage.S_ref = 0;
constFuselage.A_base=0;
constFuselage.Wingspan=0;
%% Calculating reynolds number

Re_1 = (constFuselage.rho * constFuselage.vel * constFuselage.length) / constFuselage.mu;

Re_cutoff = 38.21 * (constFuselage.length / constFuselage.k)^1.053; 

if Re_1 < Re_cutoff
    constFuselage.Re = Re_1;
else 
    constFuselage.Re = Re_cutoff;
end

%% Calculating the fineness ratio

constFuselage.f = constFuselage.length / sqrt((4 * constFuselage.Amax) / pi);

%% Calculating Mach

constFuselage.Mach = constFuselage.vel / constFuselage.a;

%% Calculating S_wet

 % wetted area estimate, m^2


end
clc; close all; clear all;

constFuselage = calcFuselageDerivedInputs();
S_wet_body = zeros(3,1);
cf_body = zeros(3,1);
FF_body = zeros(3,1);
cf_wing = zeros(3,1);
FF_wing = zeros(3,1);
S_wetwing = zeros(3,1);
Dq = zeros(3,1);
Cd0_misc = zeros(3,1);
Cd0_tot = zeros(3,1);

for i=1:3
S_wet_body(i) =2.*pi.*constFuselage.radius(i).*constFuselage.length(i); % wetted area of body (m^2)
cf_body(i) = .455./((log10(constFuselage.Re(i))).^2.58 .* (1+(.144.*constFuselage.Mach(i)).^2).^.65); % skin friction over fuselage (assume all turbulent)
FF_body(i) = .9+5./(constFuselage.f(i)).^1.5+constFuselage.f(i)./400; % fuselage form factor
cf_wing(i) = .074 ./ ((constFuselage.ReWing(i)).^.2); % skin friction over wing (assume all turbulent)
FF_wing(i) = (1+(.6./constFuselage.Xc(i)).*(constFuselage.t(i)./constFuselage.c(i))+100.*(constFuselage.t(i)./constFuselage.c(i)).^4).*((1.35*constFuselage.Mach(i)).^.18.*(cos(constFuselage.lambda(i))).^.28);
S_wetwing(i) = 2 .* constFuselage.wingPlanformArea(i); % wetted area of wing (m^2)
% Dq(i) = (.139+.419.*(constFuselage.Mach(i)-.161).^2).*constFuselage.A_base(i);
% Cd0_misc(i) = Dq./constFuselage.S_ref(i);
Cd0_tot(i) = ((cf_body(i)*FF_body(i)*S_wet_body(i))+(FF_wing(i)*constFuselage.Q_wing(i)*cf_wing(i)*S_wetwing(i)))/constFuselage.S_ref(i);
end

%% Function
function constFuselage = calcFuselageDerivedInputs()% this function calculates all the derived inputs needed for the fuselage
%% Defining Given Constants
% Standard Atmosphere (assume sea level)
constFuselage.alt = [0, 4000, 13000]; % altitude(m)
for i=1:3
[T(i), constFuselage.a(i), P(i), constFuselage.rho(i)] = atmoscoesa(constFuselage.alt(i)); % temp, speed of sound, pressure, density
end

constFuselage.vel = [21,63.9,260.8]; % velocity(m/s)
constFuselage.length = [1.56,8.33,70.7]; % length(m)
constFuselage.mu = [1.7894*(10^(-5)),1.6612*(10^(-5)),1.4216*(10^(-5))]; % viscosity(kg/m*s)
constFuselage.k = [.052 * 10^-5, .634*10^-5, .634*10^-5]; % surface k value
constFuselage.Amax = [0.02,1.13,33.18]; % max area(m^2)
constFuselage.radius=[0.08,0.6,3.25]; % radius(m)
constFuselage.Q_body=[1,1,1];
constFuselage.c = [.23,1.472,9.27]; % chord length (m) (cessna and 747 average chord)
constFuselage.Xc= [.302,.3,.4]; % normalized chord position of max thickness
constFuselage.t = [constFuselage.c(1) * .087, constFuselage.c(2) * .12,constFuselage.c(3) * .101]; % airfoil max thickness
constFuselage.lambda= [0, 0, 37.5]; % sweep angle (degrees)
constFuselage.Q_wing = [1.4, 1.4, 1.4]; % Q value for wing, assumed 1.4 for all since values unknown, this is upper limit so conservative estimate
constFuselage.Width=[.16, 1.2, 6.5]; % width of fuselage
constFuselage.bodyPlanformArea= constFuselage.Width .* constFuselage.length;
constFuselage.wingPlanformArea = [.63, 16.17, 511]; %(m^2)
constFuselage.S_ref = constFuselage.bodyPlanformArea + constFuselage.wingPlanformArea; % area of .Wiplane's shadow (m^2)
% constFuselage.A_base=0; none of the planes have an aggressively sloped
% aft end

constFuselage.Wingspan=[3.22,11,59.6]; % m
%% Calculating fuselage reynolds number

for i = 1:3
Re_1 = (constFuselage.rho(i) .* constFuselage.vel(i) .* constFuselage.length(i)) ./ constFuselage.mu(i);

Re_cutoff = 38.21 * (constFuselage.length(i) / constFuselage.k(i))^1.053; 

if Re_1 < Re_cutoff
    constFuselage.Re(i) = Re_1;
else 
    constFuselage.Re(i) = Re_cutoff;
end

%% Calculating wing reynolds number

constFuselage.ReWing(i) = (constFuselage.rho(i) .* constFuselage.vel(i) .* constFuselage.c(i)) ./ constFuselage.mu(i);

%% Calculating the fineness ratio

constFuselage.f(i) = constFuselage.length(i) / sqrt((4 * constFuselage.Amax(i)) / pi);

%% Calculating Mach

constFuselage.Mach(i) = constFuselage.vel(i) / constFuselage.a(i);


end
end

function constFuselage = calcFuselageDerivedInputs()
% this function calculates all the derived inputs needed for the fuselage
%% Defining Given Constants
% Standard Atmosphere (assume sea level)
constFuselage.alt = 0; % altitude
[T, constFuselage.a, P, constFuselage.rho] = atmoscoesa(constFuselage.alt) % temp, speed of sound, pressure, density


constFuselage.vel = 0; % velocity
constFuselage.length = 0; % length
constFuselage.mu = 0; % viscosity
constFuselage.k = 0; % surface whatever
constFuselage.Amax = 0; % max area


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

S_wet = 2.26; % wetted area estimate, m^2


end
%% U.S. Standard Atmosphere of 1976 v1.2
%  Written by J.A. Ferrand B.Sc (ID: 2431646)
%  Embry-Riddle Aeronautical University - Daytona Beach
%  College of Engineering (COE)
%  For use in AE 308, AE 403, AE 435, AE 440, and other Aerospace
%  Engineering (AE) coursework.
%% Description
% Not to be confused with the International Standard Atmosphere (ISA), U.S.
% Standard Atmosphere of 1976 was the joint work of NOAA, NASA, and the
% USAF to provide a reference atmospheric model for civilian and military
% applications. The 1976 model is derived from the ideal gas law and the
% hydrostatic equation. The atmospheric model, comprising of values for
% temperature, pressure, and density consists of piecewise functions which
% curve fit atmospheric data. The temperature variations are assumed to be
% linear, and the fluid is taken to be air. This function includes modified
% versions of the original 1976 model that account for atmospheric
% fluctuations in temperature. Sub models for "hot," "cold," and "tropical"
% days developed in MIL-STD-210C are included. The function performs all
% calculations using SI units, optional conversion to BE units can be
% specified.
%% Syntax
% * INPUT(*h*): Height vector. (Stations at which to evaluate atmosphere)
% * INPUT(*units*): String denoting the unit system used. Either "BE" for
% British Engineering or "SI" for systemme Internationale. (Default is SI)
% * INPUT(*day*): String denoting the type of nonstandard atmospheric
% model. Available options are "Hot," "Cold," and "Tropical."
% * OUTPUT(*T*): Temperature values.
% * OUTPUT(*P*): Pressure values.
% * OUTPUT(*r*): Density values.
%% Formulae
% * $$g_{0} = 9.80665 \frac{m}{s^{2}} \quad r_{0} = 6356577 m \quad R_{u} =
% 8.31432 \frac{J}{mol-K} \quad MW_{air} = 28.9644 \frac{g}{mol-K}$
% * $$T = T_{i} + L_{i}(h- h_{i})$
% * $$z = \frac{(h)(r_{0})}{h + r_{0}}$
% * $$P = \begin{array}{c} P_{i}\left(\frac{T_{i}}{T}\right)^{\frac{g_{0}
% W_{0}}{R L_{i}}}\\ P_{i}e^{\frac{-g_{0}W_{0}(z-z_{i})}{RT_{i}}}
% \end{array} \quad   \begin{array}{c} L_{i} = 0 \\ L_{i} \ne 0
% \end{array}$
% * $$\rho = \frac{MW_{air}P}{R_{u}T}$
%% Data Visualization
% <<atmos1976.PNG>>
%% Changelog
%  v1.2,(09/25/2021): Incorporated the non-standard temperature models from
%  appendix B of AED 2nd ed. Cleaned up the script.
%  v1.1,(05/08/2021): Added unit conversion capability from SI to BE.
%  Hardcoded piecewise values of the atmospheric regions. Replaced "if"
%  statements and "for" loops with boolean and linear index vectorizations.
%  v1.0,(04/08/2021): Initial release. All four (4) atmospheric models from
%  AED 2nd ed. are programmed. Perfectly replicates ATMOS in the standard
%  day, notable discrepancies in the Hot, Cold, and Tropical days.
%% Function Definition
function [T,P,rho] = atmos(h,units,day)
%Error checking on altitude input.
if sum(h(h < 0)) > 0 || isnumeric(h) ~= 1
    error('Invalid altitude input entered!')
else
    h = h(:); %Flatten the altitude input.
end
if h > 85722.64420076866
    warning('Input altitude is outside standard atmospheric model''s range [h < 85,723 [m](281,243 [ft])]. Some outputs are extrapolated.')
end

%Program defaults
if nargin < 2
    units = 'SI'; %Default unit system is Systeme International (SI).
end
if nargin < 3
    day = 'Standard'; %Default atmospheric model is the standard day.
end

%Unit conversions
if strcmpi(units,'BE') == 1
    h = h*0.3048; % Input altitude in BE units converted to SI units. [m]
    Tf = 1.8; %K to R
    Pf = 0.020885434273039; %Pa to psf
    rf = 0.001940333398433; %kg/m3 to slug/ft3
elseif strcmpi(units,'SI') == 0
    error('Invalid system of units! (BE or SI only)')
end

%Memory preallocation
[N,c] = size(h); %Obtain the number of altitude queries.
P = zeros(N,c); %Preallocate memory for the Pressure output vector.
region = zeros(N,c); %Preallocate memory for a vectorization.
I = (1:N)'; %Pre-compute the linear indices of the output arrays.

%Physical constants
g0W0oR = 0.034163194736310; %"g0 times MW divided by Ru"
W0oR = 0.003483676355974;%"MW divided by Ru"[kg-K/J]
r0 = 6356577; % Earth's radius. [m]

%Isothermal/lapsing region vectorizations:
z = h*r0./(r0 + h); % Geopotential altitude [m]
region(z >= 0.000 & z < 11000) = 1; %Region 1
region(z >= 11000 & z < 20000) = 2; %Region 2
region(z >= 20000 & z < 32000) = 3; %Region 3
region(z >= 32000 & z < 47000) = 4; %Region 4
region(z >= 47000 & z < 51000) = 5; %Region 5
region(z >= 51000 & z < 71000) = 6; %Region 6
region(z >= 71000 & z < 84582) = 7; %Region 7
iso = region == 2 | region == 5; %Boolean array for isothermal regions.
lap = ~iso; %Boolean array for regions with temperature gradients.
iso_reg = region(iso); %Isothermal region indices.
lap_reg = region(lap); %Lapsing region indices.

%Standard day stations of...
Ti = [288.15;216.65;216.65;228.65;270.65;270.65;214.65]; %Temperature [K]
Li = [-0.0065;0;+0.001;+0.0028;0;-0.0028;-0.002]; %Lapse rate [K/m]
Pi = [101325;22632.06397346292;5474.888669677776;868.0186847552283;...
    110.9063055549660;66.938873118687326;3.956420428040731]; %Pressure [Pa]
zi = [0;11000;20000;32000;47000;51000;71000]; %Geopotential [m]

%Compute standard temperature and pressure:
T = (Ti(region) + Li(region).*(z - zi(region)));%Temperature @ input altitude. [K]
P(I(iso)) = Pi(iso_reg).*exp(-g0W0oR*(z(I(iso)) - zi(iso_reg))./Ti(iso_reg));
P(I(lap)) = Pi(lap_reg).*(Ti(lap_reg)./T(I(lap))).^(g0W0oR./Li(lap_reg));

%Nonstandard atmospheric temperature models.
switch day
    case 'Standard'
        %No temperature correction is applied
    case {'Hot','Cold','Tropical'}
        if h > 35000
            warning('Input altitude is outside nonstandard atmospheric model''s range [h < 35,000 [m](114,829 [ft])]. Some outputs are extrapolated.')
        end
        switch day
            case 'Hot' %day stations of...
                Ti = [312.6;228.6;235.4]; %Temperature [K]
                Li = [-0.007;+0.008;+0.014]; %Lapse rate [K/m]
                hi = [0;12000;20500]; %Altitude [m]
            case 'Cold' %day stations of...
                Ti = [222.1;247.1;247.1;208.1;208.1;185.9;185.9;204.3]; %Temperature [K]
                Li = [+0.025;0;-0.006;0;-0.00888;0;+0.0046;-0.000775]; %Lapse rate [K/m]
                hi = [0;1000;3000;9500;13000;15500;18500;22500]; %Altitude [m]
            case 'Tropical' %day stations of...
                Ti = [305.27;193.27;212.27]; %Temperature [K]
                Li = [-0.007;+0.0038;+0.00248]; %Lapse rate [K/m]
                hi = [0;16000;21000]; %Altitude [m]
        end
        %Overwrite standard day temperature with nonstandard one.
        T = (Ti(region) + Li(region).*(h - hi(region))); %Temperature @ input altitude. [K]
    otherwise
        error('Invalid atmospheric model input.')
end
rho = P*W0oR./T; %Static density @ input altitude (subject to nonstandard temperature). [kg/m^3]
if strcmpi(units,'BE') == 1
    T = T*Tf; P = P*Pf; rho = rho*rf; %Apply unit conversion factors if BE units.
end
end
%% Sources
% * "Aircraft Engine Design," 2nd ed. (2002) by Jack D. Mattingly, William
% H. Heiser, and David T. Pratt
% * "U.S. Standard Atmosphere, 1976," (1976) U.S. Government Printing Office,
% Washington, DC
% * "Climactic Information to Determine Design and Test Requirements for
% Military Equipment," (1997) MIL-STD-210C, Rev C, Washington, DC.
%% IG-CSH Theta-Beta-Mach Relation v1.3
%  Written by J.A. Ferrand B.Sc (ID: 2431646)
%  Embry-Riddle Aeronautical University - Daytona Beach
%  College of Engineering (COE)
%  For use in AE 308, AE 403, AE 435, AE 440, and other Aerospace
%  Engineering (AE) coursework.
%% Description
% A function that solves the the $\theta - \beta -M$ relation when given
% two of the constituent inputs. This relationship is the basis for the
% solution to oblique shock problems. This function is a building block for
% higher-level solver functions for inlets, wedges, and sliplines.
%% Required Plugins
% * NONE.
%% Formulae
% $\cot{\theta} = \tan{\beta}[\frac{(\gamma+1)M^{2}}{2(M^{2}\sin^{2}{\beta}-1)}-1]$
%% Syntax
% * INPUT(*Variable Argument(s)*): This function operates with "key-value"
% pairs. A "key" is a string that identifies the input. The "value" is the
% numeric array that corresponds the key. Valid keys: "mach" for Mach
% Number (M), "theta" for Wedge angle ($\theta$), "beta" for shockwave
% angle ($\beta$), "units" to specify angular inputs as either "deg" for
% degrees or "rad" for radians, and "gam" for ratio of specific heats
% ($\gamma$).
% * NOTE: "mach," "theta," and "beta" represent physical inputs. Only two
% of these must be specified.
% * NOTE: The default value for units is "deg."
% * OUTPUT(*out*): Of the three physical inputs (Mach Number, wedge angle,
% and shockwave angle), *out* contains the numeric values of the one that
% was not specified.
%% Changelog
%  v1.1,(09/05/2021): Vectorized inputs and improved default values.
%  v1.2,(09/07/2021): Enhanced variable argument recognition.
%  v1.3,(09/08/2021): Implemented error checking on all cases, reduced
%  bloat and redundant variables. When M and B are specified, compativility
%  is established through the "sin2Bmax" relation instead of computing the
%  discriminant.
%% Function Definition
function out = tbm(varargin)
%Variable argument logic
if mod(nargin,2) == 1 %Incomplete input pair.
    error('Uneven number of inputs! Must supply complete pairs.')
elseif nargin > 12 %Absolutely too many inputs.
    error('Excessive number of input pairs! This function deals with only five(5).')
elseif nargin < 4 %Absolutely too few inputs
    error('Insufficient input pairs! (minimum of two(2) physical inputs!)')
else %Viable number of input pairs given.
    Mspec = 0; %Declare register variable for Mach number.
    Bspec = 0; %Declare register variable for shock angle.
    THspec = 0; %Declare register variable for ramp angle.
    yspec = 0; %Declare register variable for Cp/Cv ratio.
    unispec = 0; %Declare register variable for angular units.
    solspec = 0; %Declare register variable for shock solution.
    for i = 2:2:nargin %For each argument pair.
        switch varargin{i-1} %Attempt to recognize valid inputs.
            case {'mach','M'} %If input recognized as Upstream Mach numer.
                if isnumeric(varargin{i}) == 0  %1 < M < infty
                    error('Invalid Mach number. Must be numeric.')
                else
                    M = varargin{i}(:);%Accept input and flatten array.
                    val = M >=1;%Preliminary validation condition on M.
                    M(~val) = NaN;
                    Mspec = 1;%Register that M was correctly input.
                end
            case {'beta','B'} %If input recognized as shockwave angle.
                if isnumeric(varargin{i}) == 0
                    error('Invalid input for Shockwave angle')
                else
                    B = varargin{i}(:); %Accept input and flatten array.
                    val = B > 0; %Preliminary validation condition on B.
                    B(~val) = NaN;
                    Bspec = 1; %Register that B was correctly input.
                end
            case {'theta','TH'} %If input recognized as the wedge angle.
                if isnumeric(varargin{i}) == 0
                    error('Invalid input for Wedge Angle')
                else
                    TH = varargin{i}(:); %Accept input and flatten array.
                    val  = TH > 0; %Preliminary validation condition on TH.
                    TH(~val) = NaN;
                    THspec = 1; %Register that TH was correctly input.
                end
            case {'gam','y'} %If input recognized as the ratio of specific heats.
                if isnumeric(varargin{i}) == 0
                    error('Invalid input for gamma.')
                else
                    gam = varargin{i};%Accept input.
                    yspec = 0;
                end
            case 'units' %If input recognized as sytem of angular units.
                if strcmpi(varargin{i},'deg') == 0 && strcmpi(varargin{i},'rad') == 0
                    error('Invalid input for units')
                else
                    units = varargin{i};
                    unispec = 1;
                end
            case 'sol' %If input recognized as choice between shock solutions.
                if strcmpi(varargin{i},'strong') == 0 && strcmpi(varargin{i},'weak') == 0
                    error('Enter either "weak" or "strong" for shock solution.')
                else
                    sol = varargin{i};
                    solspec = 1;
                end
            otherwise %Input cannot be recognized.
                error('Unrecognizable input detected.')
        end
    end
    if Mspec + Bspec + THspec == 3 %If all three(3) variables specified.
        error('Three(3) physical variables input. System is overconstrained!')
    elseif Mspec + Bspec + THspec == 1 %If only one(1) variable specified.
        error('Only one(1) physical variable input. System is underconstrained!')
    end
    if solspec == 0 %If no shock solution specified.
        sol = 'weak'; %Default to the weak solution.
    end
    if  unispec == 0 %If no angular units specifed.
        units = 'deg'; %Default to degrees.
    end
    if yspec == 0 %If no ratio of specific heats was specified.
        gam = 1.4; %Default to the one for air.
    end
    gp1 = gam + 1;
    if strcmpi(units,'deg') == 1 %If working with degrees
        sine = @sind; cotan = @cotd; arctan = @atand;
        arccos = @acosd; arcsin = @asind; tangent = @tand;
        nineT = 90; %orthogonal angle [deg]
    else %If working with radians.
        sine = @sin; cotan = @cot; arctan = @atan;
        arccos = @acos; arcsin = @asin; tangent = @tan;
        nineT = pi/2; %orthogonal angle [rad]
    end
    if THspec == 1 %If one of the physical inputs is the wedge angle.
        val = TH <= arctan(sine(arccos(-1/gam))/(gam-1/gam));
        TH(~val) = NaN; %Additional preliminary validation condition on TH.
    end
end
%Procedural Solution:
if THspec == 0 %(Mach number and shockwave angle specified)
    val = arcsin(1./M) <= B & B <= nineT; %Validity condition
    if sum(~val) > 0%If non-physical inputs detected.
        warning('Non-physical shockwave angles detected and ignored.')
        out(~val) = NaN;
    end
    M2 = M(val).^2;
    sin2B = (sine(B(val))).^2;
    out(val) = arctan(2*cotan(B(val)).*(M2.*sin2B - 1)./(M2.*(gp1 - 2*sin2B) + 2));
elseif Bspec == 0 %(Mach number and wedge angle specified)
    M2 = M.^2; %Mach Number squared.
    M4 = M2.*M2; %Fourth power of Mach Number.
    s2B = (M2*gp1 - 4 + sqrt(gp1*(M4*gp1+8*M2*(gam-1)+16)))./(4*M2*gam);
    val = TH < arctan(2*sqrt(1./s2B-1).*(M2.*s2B-1)./(2+M2.*(gp1-2*s2B)));
    sin2TH = (sine(TH(val))).^2;
    b = -((M2(val) + 2)./M2(val) + gam*sin2TH);
    c = (2*M2(val) + 1)./M4(val) + (gp1^2/4 + (gam-1)./M2(val)).*sin2TH;
    d = (sin2TH - 1)./M4(val);
    Q = (3.*c - b.^2)/9;
    R = (9.*b.*c - 27*d - 2*(b.^3))/54;
    %D = Q.^3 + R.^2; %Discriminant of the cubic in sin2(B).
    %valid = D < 0;
    if sum(~val) > 0
        out(~val) = NaN;
        warning('Non-physical mach-wedge pairs detected and ignored!')
    end
    if strcmpi(sol,'weak') == 1 %Compute the weak shockwave angle.
        out(val) = 2*sqrt(-Q).*cos(acos(R./(sqrt(-Q.^3)))/3 + 4*pi/3) - b/3;
    else %Compute the strong shockwave angle.
        out(val) = 2*sqrt(-Q).*cos(acos(R./(sqrt(-Q.^3)))/3) - b/3;
    end
elseif Mspec == 0 %(Shockwave angle and wedge angle specified)
    tanTH = tangent(TH);
    sin2Bt2 = 2*(sine(B).^2);
    out = sqrt(-2*(tanTH + cot(B))./(tanTH*(gp1-sin2Bt2)-cotan(B)*sin2Bt2));
end
end
%% Sources
% * "Equations, Tables, and Charts for Compressible Flow," (1953), National
% Advisory Commitee for Aeronautics, Technical Report #1135, Ames Research
% Center.
% *  Thompson, M.J., "A Note on the Calculation of Oblique Shock- Wave
% Characteristics",Journalof theAeronauticalSciences,November 1950
% * "Fundamentals of Aerodynamics," 6th ed. (2017), by John D. Anderson
% * "Aircraft Engine Design," 2nd ed. (2002) by Jack D. Mattingly, William
% H. Heiser, and David T. Pratt
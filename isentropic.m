%% IG-CSH Isentropic Relations v1.4
%  Written by J.A. Ferrand B.Sc (ID: 2431646)
%  Embry-Riddle Aeronautical University - Daytona Beach
%  College of Engineering (COE)
%  For use in AE 308, AE 403, AE 435, AE 440, and other Aerospace
%  Engineering (AE) coursework.
%% Description
% This is a utility that generates data from the isentropic relations that
% result from the application of the Energy Conservation equation to 1-D,
% Inviscid, Compressible Flow in the absence of body forces of an Ideal Gas
% (IG) with constant specific Heats (CSH). The utility analytically
% resolves the thermodynamic state (defined as the collection of Mach
% Number ($M$), Isentropic Pressure Ratio ($P_{t}/P$), Isentropic 
% Temperature Ratio ($T_{t}/T$), Isentropic Density ratio ($\rho_{t}/\rho$)
% , Area Ratio ($A/A^{*}$), and Mass Flow Parameter ($MFP$) when 
% information about any of the first four(4) properties is input. If
% the latter two(2) are input, the utility needs to call a Newton-Raphson 
% solver to resolve the thermodynamic state. The utility can optionally 
% compute the Prandtl-Meyer angle (in radians) and the specific Impulse 
% ratio ($I/I^{*}$) associated with the Mach Number.
%% Required Plugins
% * newton.m (Newton-Raphson solver needed for numeric schemes)
%% Formulae
% * $$\frac{T_{t}}{T} = 1 + \frac{\gamma-1}{2}M^{2}$
% * $$\frac{P_{t}}{P} = \left(\frac{T_{t}}{T}\right)^{\frac{\gamma}
% {\gamma-1}}$
% * $$\frac{\rho_{t}}{\rho} = \left(\frac{T_{t}}{T}\right)^{\frac{1}
% {\gamma-1}}$
% * $$\frac{A}{A^{*}} = \frac{1}{M}\sqrt{\left(\frac{2}{\gamma + 1}
% \frac{T_{t}}{T}\right)^{\frac{\gamma + 1}{\gamma-1}}}$
% * $$\nu = \arctan{\left(\sqrt{M^{2}-1}\sqrt{\frac{\gamma-1}{\gamma+1}}
% \right)} \sqrt{\frac{\gamma+1}{\gamma-1}} -\arctan{\left(\sqrt{M^{2}-1}
% \right)}$
% * $$MFP \times \sqrt{R} = M \sqrt{\gamma} \left(\frac{Tt}{T}\right)
% ^{\frac{-(\gamma+1)}{2(\gamma-1)}}$
% * $$\frac{I}{I^{*}} = (1+\gamma M^{2})\left(\frac{A}{A^{*}}\right)
% \left(\frac{P}{P_{t}}\right)\left(\frac{P_{t}}{P^{*}}\right)$
%% Data Visualization
% <<isengraph.PNG>>
%% Changelog
%  v1.3,(01/01/2021): Numeric solution to MFP formula enabled.
%  v1.2,(12/30/2020): Numeric solution to Prandtl-Meyer formula enabled.
%  v1.1,(12/24/2020): Numeric solution to Stodola's Area ratio enabled.
%  v1.0,(12/02/2020): Initial Release (evaluation of analytical isentropic states only).
%% Syntax
% * INPUT(*IN*): Array of numeric values of any of the physical outputs 
% shown below.
% * INPUT(*known*): String denoting which physical property does the input
% "IN" represent. OPTIONS: {'M','TtT','PtP','AAsub','AAsup','MFPsub',
% 'MFPsup','nu'}. Default: 'M'
% * INPUT(*gam*): Ratio of specific Heats. Default: 1.4
% * OUTPUT(*TtT*): Total-to-Static isentropic Temperature ratio.
% * OUTPUT(*PtP*): Total-to-Static isentropic Pressure ratio.
% * OUTPUT(*rtr*): Total-to-Static isentropic Density ratio.
% * OUTPUT(*AAs*): Stodola's A/A* Area ratio.
% * OUTPUT(*MFP*): Mass Flow Parameter (non-dimensional mass flow rate).
% * OUTPUT(*nu*): Prandtl-Meyer Angle [rad]
% * OUTPUT(*IIs*): I/I* Specific Impulse ratio.
%% Function Definition
function [M,TtT,rtr,PtP,AAs,nu,MFP,IIs] = isentropic(IN,known,gam)
[r,c] = size(IN);
N = r*c; %Number of queries.
if c > 1
    IN = IN(:);
elseif N == 0
    error('Empty data array input!')
elseif IN < 0
    error('All normal shock relations are positive! (Negative input detected)')
end
if nargin < 2
    known = 'M'; %Default "known" is the Mach Number.
end
if nargin < 3
    gam = 1.4; %Conventional Ratio of Specific Heats of Air.
end
gm1 = gam - 1; %Constant found everywhere.
gm1o2 = gm1/2; %Constant found in TtT, PtP, and nu's "FoD" formula.
gp1 = gam + 1; %Constant found in A/A* formula.
del = gam/gm1; %Constant for PtP's exponent.
lam2 = gm1/gp1; %Constant found in A/A* formula.
lam = sqrt(lam2); %Constant found in formulae associated with PM angle.
TtT = []; rtr = []; PtP = []; AAs = []; nu = [];  MFP = [];
switch known
    case 'M' %Mach Number known (ANALYTIC)
        M = IN;
    case 'TtT' %Isentropic Temperature Ratio known (ANALYTIC)
        TtT = IN;
        M = sqrt((IN - 1)/gm1o2);
    case 'rtr' %Isentropic Density Ratio known (ANALYTIC)
        rtr = IN;
        M = sqrt((IN.^gm1 - 1)/gm1o2); 
    case 'PtP' %Isentropic Pressure Ratio known (ANALYTIC)
        PtP = IN;
        M = sqrt((IN.^(1/del) - 1)/gm1o2);
    case {'AAsup','AAsub','MFPsub','MFPsup'} %(NUMERIC)
        if strcmpi(known(end-2,end),'sub') == 1
            X0 = ones(1,N)*0.2; %Initial guess to converge-upon subsonic solution.
        else
            X0 = ones(1,N)*2; %Initial guess to converge upon supersonic solution.
        end
        M = ones(1,N)*NaN; %Preallocate memory for output Mach Numbers.
        if strcmpi(known(1:2),'AA') == 1 %Stodola's A/A* ratio known 
            AAs = IN;
            valid = AAs >= 1;
            [M(valid),~,~] = newton(IN(valid),@AM_NR,'FoD',X0,...
                struct('lam2',lam2,'gp1',gp1,'gm1',gm1));
        else %Mass Flow Parameter known 
            MFP = IN;
            valid = MFP > 0 & MFP < sqrt(gam)*((2/gp1)^(0.5/lam2));
            [M(valid),~,~] = newton(IN(valid),@MFP_NR,'FoD',X0,...
                struct('lam2',lam2,'gp1',gp1,'gam',gam));
        end
    case 'nu' %PM angle (in radians) known. (NUMERIC)
        nu = IN;
        nu_max = (1/lam - 1)*pi/2; %Compute Theoretical Maximum PM angle (safety).
        valid = (IN < nu_max & IN > 0); %PM angle is valid only if less than theoretical maximum.
        M = NaN(1,N); %Preallocate memory for output Mach Numbers.
        if sum(valid) < length(IN)
            warning('One or more input Prandtl-Meyer angles are non-physical.')
        end
        X0 = Hall_approx(IN(valid),nu_max); %Use Hall's approximation only on valid PM angles.
        [M(valid),~,~] = newton(IN(valid),@PM_NR,'FoD',X0,struct('lam',lam,'gm1o2',gm1o2));
    otherwise
        error('Invalid "known" input. Valid options: M, PtP, TtT, AAsub, AAsup, MFPsub, MFPsup, and nu')
end
M2 = M.^2;
if isempty(TtT) == 1 && nargout > 1
    TtT = 1 + M2*gm1o2; %Compute Isentropic Temperature Ratios (if needed).
end
if isempty(rtr) == 1 && nargout > 2
    rtr = TtT.^(1/gm1); %Compute Isentropic Temperature Ratios (if needed).
end
if isempty(PtP) == 1 && nargout > 3
    PtP = TtT.^del; %Compute Isentropic Pressure Ratios (if needed).
end
if isempty(AAs) == 1 && nargout > 4
    AAs = ((2*TtT./gp1).^(0.5/lam2))./M; %Compute A/A* Ratios (if needed).
end
if isempty(nu) == 1 && nargout > 5
    supersonic = M2 >= 1; %PM angle only exists for supersonic flows.
    B = sqrt(M2(supersonic) - 1); %Precompute "Beta"
    nu = ones(r,1)*NaN; %Preallocate memory for PM-angles.(NaN for all subsonic inputs)
    nu(supersonic) = atan(B*lam)/lam - atan(B);
end
if isempty(MFP) == 1 && nargout > 6 %Compute MFP (if needed).
    MFP = sqrt(gam)*((2/gp1)^(0.5/lam2))./AAs;
end
if nargout > 7 %Compute Specific Impulse (if needed).
    IIs = (1 + gam*M2).*AAs./(0.528281787717174*gp1.*PtP);
end
% Newton-Raphson "FoD" for Numeric solution to the Area-Mach Input.
    function AAsFoD = AM_NR(AAs,M,P)
        tau = 2/P.gp1 + (M.^2)*P.lam2; %Equivalent to TtT*2/(gam + 1)
        f = (tau.^(0.5/P.lam2))./M; %Stodola's Area Ratio.
        AAsFoD = (f - AAs)./(tau.^(1/P.gm1 - 0.5) - f./M);
    end
% Newton-Raphson "FoD" for Numeric solution to the PM angle Input.
    function nuFoD = PM_NR(nu,M,P)
        Msq = M.^2;
        BETA = sqrt(Msq - 1);
        nuFoD = (atan(BETA*P.lam)/P.lam - atan(BETA) - nu)...
            .*(M.*(1 + Msq*P.gm1o2))./BETA;
    end
% Newton-Raphson "FoD" for Numeric solution to the MFP Input.
    function MFPFoD = MFP_NR(MFP,M,P)
        Msq = M.^2;
        tau = 2/P.gp1 + Msq*P.lam2; %Equivalent to TtT*2/(gam + 1)
        f = (tau.^(0.5/P.lam2))./M; %Stodola's Area Ratio.
        MFPFoD = (M - MFP*f*sqrt(P.gam))/(2 - Msq./tau);
    end
% I.A Hall's Approximate Inverse Prandtl-Meyer Function.
    function M_aprx = Hall_approx(nu,nu_max)
        y = (nu./nu_max).^(2/3);
        M_aprx = (1 + y.*(1.3604 + y.*(0.0962 - 0.5127*y)))...
            ./(1 - y.*(0.6722 + 0.3278*y));
    end
end
%% Sources
% * "Fundamentals of Aerodynamics," 6th ed. (2017), by John D. Anderson
% * "Elements of Gas Turbine Propulsion," 2nd ed. (1996) by Jack D. Mattingly
% * "Inversion of the Prandtl-Meyer Relation," (1975), Aeronautical Journal Volume 79 Issue 777, pp417-418 by I. M. Hall
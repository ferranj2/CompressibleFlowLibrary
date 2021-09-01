%% IG-CSH Normal Shock (NS) Relations v1.1
%  Written by J.A. Ferrand B.Sc (ID: 2431646)
%  Embry-Riddle Aeronautical University - Daytona Beach
%  College of Engineering (COE)
%  For use in AE 308, AE 403, AE 435, AE 440, and other Aerospace
%  Engineering (AE) coursework.
%% Description
% This is a utility that generates discrete data from the Normal Shock (NS)
% Relations that result from enforcing 1-D conservation laws (mass,
% momentum, and energy) to a control volume where the fluid is inviscid,
% adiabatic, independent of body forces, comprised of a single species, 
% constant specific heats (CSH), and the equation of state is the ideal 
% gas (IG) law. Of interest in the NS problem are six(6) quantities: The 
% Upstream Mach Number ($M_{1}$), the Downstream Mach Number ($M_{2}$), the
% Static Pressure and Temperature jumps ($T_{2}/T_{1}$ and $P_{2}/P_{1}$), 
% and the Stagnation Pressure Loss ($P_{t2}/P_{t1}$). The NS problem can be 
% analytically solved for a given ratio of specific heats ($\gamma$) when 
% either $M_{1}$, $M_{2}$, or $P_{2}/P_{1}$ are known. If not, then the NS 
% problem must be solved numerically. 
%% Syntax
% * INPUT(*IN*): Numeric value of any of the physical outputs shown
% below.
% * INPUT(*known*): String denoting which physical property does the input
% "IN" represent.
% * INPUT(*gam*): Ratio of specific Heats.
% * OUTPUT(*M1*): Upstream Mach Number.
% * OUTPUT(*M2*): Downstream Mach Number.
% * OUTPUT(*p2_p1*): Static Pressure ratio across N.S.
% * OUTPUT(*T2_T1*): Static Temperature ratio across N.S.
% * OUTPUT(*r2_r1*): Static Density ratio across N.S.
% * OUTPUT(*pt2_pt1*): Stagnation Pressure ratio across N.S.
%% Required Plugins
% * newton.m (To enable numeric solutions)
% * custompolyval.m (To evaluate initial guess schema for newton.m)
%% Formulae
% * $$M_{2} = \sqrt{\frac{M_{1}{^2}\, \left(\gamma - 1\right) + 2}
% {2\gamma M_{1}{^2} - \gamma + 1}}$
% * $$\frac{P_{2}}{P_{1}} =  1 + 2\gamma\frac{M_{1}^{2}-1}{\gamma + 1}$
% * $$\frac{\rho_{2}}{\rho_{1}} = \frac{(\gamma+1)M_{1}^{2}}{(\gamma-1)
% M_{1}^{2} + 2}$
% * $$\frac{T_{2}}{T_{1}} =
% \left(\frac{P_{2}}{P_{1}}\right)\left(\frac{\rho_{1}}{\rho_{2}}\right)$
% * $$\frac{P_{t2}}{P_{t1}} =
% \left[\frac{\rho_{2}}{\rho_{1}}\right]^{\frac{\gamma}{\gamma-1}}
% \left[\frac{P_{1}}{P_{2}}\right]^{\frac{1}{\gamma-1}}$
% * $$\frac{P_{t2}}{P_{1}} = \left[\frac{(\gamma + 1)^{2} M_{1}^{2}}{4 \gamma
% M^{2}_{1} - 2(\gamma-1)}\right]^{\frac{\gamma}{\gamma-1}} \frac{1 - \gamma + 
% 2\gamma M_{1}^{2}}{\gamma + 1}$
%% Changelog
%  v1.1,(03/20/2021): Support for numeric solutions to T2_T1 input enabled.
%  v1.0,(12/10/2020): Initial Release (evaluation of analytical N.S states only).
%% Function Definition
function [M1,M2,p2p1,r2r1,T2T1,pt2pt1] = nshock(IN,known,gam)
[r,c] = size(IN);
N = r*c; %Number of queries.
if c > 1
    IN = IN(:);
elseif N == 0
    error('Empty data array input!')
elseif IN < 0
    error('All normal shock relations are positive! (Negative input detected)')
end
gm1 = gam - 1; % "gamma minus one"
gp1 = gam + 1; % "gamma plus one"
gt2 = gam*2; % "twice gamma"
M1 = NaN(1,N); %Preallocate memory for output Mach Numbers.
M2 = []; p2p1 = []; T2T1 = []; r2r1 = []; pt2pt1 = [];
switch known
    case 'M1' %Upstream Mach Number (M1) Known (ANALYTIC)
        M1 = IN;
    case 'M2' %Downstream Mach Number (M2) Known (ANALYTIC)
        M2 = IN;
        lim = sqrt(gm1/gt2); %M2 converges to a lower limit.
        valid = (M2 > lim & M2 < 1); %Bounds of M2.
        if sum(valid) < N
            warning(['Non-physical input (M2 > 1) or (M2 < ',num2str(lim),') detected!'])
        end
        M22 = M2(valid).^2;
        M1(valid) = sqrt((M22*gm1 + 2)./(M22*gt2 - gm1));
    case 'p2_p1' %Static Pressure Ratio Known (ANALYTIC)
        p2p1 = IN;
        valid = p2p1 >= 1; %Static Pressure ALWAYS rises across N.S.
        if sum(valid) < N
            warning('Non-physical input (p2_p1 < 1) detected!')
        end
        M1(valid) = sqrt(1 + (p2p1(valid) - 1)*gp1/gt2); %Evaluate M1.
    case 'r2_r1' %Static Density Ratio Known (ANALYTIC)
        r2r1 = IN;
        lim = gp1/gm1; %Static density ratio converges to a finite value.
        valid = r2r1 < lim & r2r1 >= 1; %Bounds of r2_r1.
        if sum(valid) < N
            warning(['Non-physical input (r2_r1 < 1) or (r2_r1 > ',num2str(lim),') detected!'])
        end
        M1(valid) = sqrt(2*r2r1(valid)./(gp1 - r2r1(valid)*gm1));
    case 'T2_T1' %Static Temperature Ratio Known (NUMERIC)
        T2T1 = IN;
        valid = T2T1 >= 1; %Static Temperature ALWAYS rises across N.S.
        if sum(valid) < N
            warning('Non-physical input (T2_T1 < 1) detected!')
        end
        %Coefficients of Polynomial/Laurent curve fit:
        C = [-1.569663858090678,4.797476021925396,-6.332860908258899,...
            3.692103419528282,0.417578279341740,-0.004632954445839];  
        %Evaluate initial guess:
        X0 = custompolyval(T2T1(valid)',C,2);
        %Callout to NR solver:
        [M1(valid),~,~] = newton(T2T1(valid),@T2oT1NR,'FoD',X0,...
            struct('gt2',gt2,'gp1',gp1,'gm1',gm1,'gam',gam));
    case 'pt2_pt1' %Stagnation Pressure Ratio Known
        pt2pt1 = IN;
    case 'pt2_p1' %Pitot Pressure ratio Known
        pt2p1 = IN;
    otherwise
        error('Invalid "known" input. Valid options: M1, M2, p2_p1, r2_r1, T2_T1, and pt2_pt1')
end
M12 = M1.^2; % Square of the Upstream Mach Number vector.
if isempty(M2) == 1 && nargout > 1
    M2 = sqrt((1 + gm1.*M12./2)./((gam.*M12) - gm1./2));
end
if isempty(p2p1) == 1 && nargout > 2
    p2p1 = 1 + gt2.*(M12 - 1)./gp1; 
end
if isempty(r2r1) == 1 && nargout > 3
    r2r1 = gp1*M12./(gm1.*M12 + 2);
end
if isempty(T2T1) == 1 && nargout > 4
    T2T1 = p2p1./r2r1;
end
if isempty(pt2pt1) == 1 && nargout > 5
    pt2pt1 = (r2r1.^(gam/gm1))./(p2p1.^(1/gm1));
end
% Newton-Raphson "FoD" for Numeric solution to the T2_T1 Input.
    function update = T2oT1NR(T2_T1,M,P)
        M1s = M.^2;
        PI = 1 + P.gt2*(M1s - 1)/P.gp1; %Equivalent to p2_p1
        SIG = 2./(P.gp1*M1s) + P.gm1/P.gp1; %Equivalent to 1/r2_r1
        update = (PI.*SIG - T2_T1)*P.gp1./(2*P.gt2*M.*(SIG - PI./(P.gam*(M1s.^2))));
    end
end
function [info,pars] = BCDSNL_main(matfile, pars)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file is a component of the software package BCDSNL 
% Copyright (C) 2022 Mitsuhiro Nishijima and Kazuhide Nakata
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Starting Prepare -->
load(matfile);
a = xMatrix0(:,noOfSensors+1:noOfSensors+noOfAnchors);
% <-- End Prepare

% Starting Main -->
eTimeMain = 0;
startingTime = tic;

% Initialization of U and V
[Unew, Vnew] = calc_InitialPoint(a, saDistMat);
Uold = Unew; Vold = Vnew;
clear saDistMat

% Calculation of initial gamma
coefmax2 = max(4*sum(S(1:noOfSensors,1:noOfssEdges)~=0,2) + sum(S(1:noOfSensors,noOfssEdges+1:noOfssEdges+noOfsaEdges),2));
coefmax = sqrt(full(coefmax2));
fvalnew = calc_fUV(Unew,Vnew,a,S,ssaDistance);
fvalold = fvalnew;
gamma = sqrt(2*fvalnew)*coefmax/2 * 5e-3;
gammaold = gamma;
clear coefmax2

p = 0; % the number of iteration
is_enter_step2 = 0;
is_step1_gammafix = 0;

startingTime = tic;
while 1
    p = p + 1;
    if is_step1_gammafix == 0
        gamma2old = gammaold; gammaold = gamma; 
        if p > 2
            if (fvalold-fvalnew)/fvalold >= (fval2old-fvalold)/fval2old
                gamma = gammaold * (gammaold/gamma2old);
            else
                gamma = gamma2old;
            end
        elseif p == 2
            gamma = gammaold/2;
        end
        % Optimization with respect to U
        for i = 1:noOfSensors
            idx = IncidenceIdx{i};
            [A,b] = calc_Ab(Unew, Vnew, a, S, ssaDistance, gamma, i, idx, sDim);
            Unew(:,i) = solve_Ab(A,b,sDim);
        end
        % Optimization with respect to V
        for i = 1:noOfSensors
            idx = IncidenceIdx{i};
            [A, b] = calc_Ab(Vnew, Unew, a, S, ssaDistance, gamma, i, idx, sDim);
            Vnew(:,i) = solve_Ab(A,b,sDim);
        end
        % Calculation of fvalとFval
        fval2old = fvalold; fvalold = fvalnew;
        fvalnew = calc_fUV(Unew,Vnew,a,S,ssaDistance);

        if  abs(fvalold-fvalnew)/fvalold < 1e-2
            is_enter_step2 = 1;
        end
    end
    if is_step1_gammafix == 1
        % Optimization with respect to U
        for i = 1:noOfSensors
            idx = IncidenceIdx{i};
            [A,b] = calc_Ab(Unew, Vnew, a, S, ssaDistance, gamma, i, idx, sDim);
            Unew(:,i) = solve_Ab(A,b,sDim);
        end
        % Optimization with respect to V
        for i = 1:noOfSensors
            idx = IncidenceIdx{i};
            [A, b] = calc_Ab(Vnew, Unew, a, S, ssaDistance, gamma, i, idx, sDim);
            Vnew(:,i) = solve_Ab(A,b,sDim);
        end
    end
    
    eTimeMain = toc(startingTime);
    UVdiff = norm(Unew-Vnew, 'fro'); Udiff = norm(Uold-Unew, 'fro'); Vdiff = norm(Vold-Vnew, 'fro');
    UVdiffrel = 2*UVdiff/(norm(Unew,'fro')+norm(Vnew,'fro')); Udiffrel = Udiff/norm(Uold,'fro'); Vdiffrel = Vdiff/norm(Vold,'fro');
    if is_step1_gammafix == 0
        if (UVdiffrel < pars.eps && Udiffrel < pars.eps && Vdiffrel < pars.eps) || eTimeMain > 3600 || is_enter_step2 == 1
            is_step1_gammafix = 1;
            UVmean = (Unew + Vnew) / 2;
            Unew = UVmean; Vnew = UVmean;
            Uold = UVmean; Vold = UVmean;
            clear UVmean
            fvalmean = calc_fUV(Unew,Vnew,a,S,ssaDistance);
            gamma = sqrt(2*fvalmean)*coefmax/2;
        else
            Uold = Unew; Vold = Vnew;
        end
    else 
        % Confirm whether a termination condition is satisfied
        if (UVdiffrel < pars.eps && Udiffrel < pars.eps && Vdiffrel < pars.eps) || eTimeMain > 3600
            break;
        else
            Uold = Unew; Vold = Vnew;
        end
    end     
end

info.eTimeMain = eTimeMain;
info.iter = p;
info.Udiff = Udiff;
info.Vdiff = Vdiff;
% <-- end Main

% Calculation of RMSD and \|U-V\|_F
RMSD = sum(vecnorm(Vnew(1:sDim,:) - xMatrix0(:,1:noOfSensors)).*vecnorm(Vnew(1:sDim,:) - xMatrix0(:,1:noOfSensors)));
info.RMSD = sqrt(RMSD/noOfSensors); 
info.UVdiff = norm(Unew-Vnew, 'fro');



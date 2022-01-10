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

pars.eps = 1e-5;
matfile = "/data/d2_s1000_a100_n0e+00_r1e-1.mat";
[info,pars] = BCDSNL_main(matfile, pars);                 

% Output
fprintf('|U-V|_F: %6.2e\n', info.UVdiff);
fprintf('RMSD: %6.2e\n', info.RMSD);
fprintf('The number of iteration: %d\n',info.iter);
fprintf('CPU time: %8.2f\n', info.eTimeMain);
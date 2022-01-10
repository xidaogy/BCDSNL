function u = solve_Ab(A,b,sDim)
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

if sDim == 2
    a11 = A(1,1); a12 = A(1,2); a21 = A(2,1); a22 = A(2,2);
    invA = [a22, -a12; -a21, a11]/(a22*a11-a12*a21);
    u = invA*b;
elseif sDim == 3
        a11 = A(1,1); a12 = A(1,2); a13 = A(1,3);
        a21 = A(2,1); a22 = A(2,2); a23 = A(2,3);
        a31 = A(3,1); a32 = A(3,2); a33 = A(3,3);
        invA = [a22*a33-a23*a32, -(a12*a33-a13*a32), a12*a23-a13*a22;...
            -(a21*a33-a23*a31), a11*a33-a13*a31, -(a11*a23-a13*a21);...
            a21*a32-a22*a31, -(a11*a32-a12*a31), a11*a22-a12*a21]/...
            (a11*a22*a33+a12*a23*a31+a13*a21*a32-a13*a22*a31-a12*a21*a33-a11*a23*a32);
        u = invA*b;
else
    u = linsolve(A,b);
end
  
    

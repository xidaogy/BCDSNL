function [U,V] = calc_InitialPoint(a, saDistMat)
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

sDim = size(a, 1);
m = size(saDistMat, 1);
aMin = 0.5*min(a,[],2);
aMax = 0.5*max(a,[],2);
V = zeros(sDim,m);
for i = 1:m
    idxSet = find(saDistMat(i,:)~=0); % anchors connected to sensor i
    if ~isempty(idxSet)
        [~, minIdx] = min(saDistMat(i,idxSet));
        V(:,i) = a(:,idxSet(minIdx));
    else
        V(:,i) = aMin + aMax;
    end
end
U = V;
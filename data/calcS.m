function [S,ssaDistance] = calcS(distanceMatrix, noOfSensors, noOfAnchors)
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

saDistMat = distanceMatrix(:,noOfSensors+1:noOfSensors+noOfAnchors);
[SensorIndexrow, SensorIndexcolumn, ssDistance] = find(triu(distanceMatrix(:,1:noOfSensors)));
noOfssEdges = length(SensorIndexrow);
S1 = sparse(SensorIndexrow, 1:noOfssEdges, 1, noOfSensors+noOfAnchors, noOfssEdges) ...
    - sparse(SensorIndexcolumn, 1:noOfssEdges, 1, noOfSensors+noOfAnchors, noOfssEdges);
[SensorIndex, AnchorIndex, saDistance] = find(saDistMat);
noOfsaEdges = length(SensorIndex);
S2 = sparse(SensorIndex,1:noOfsaEdges, 1, noOfSensors+noOfAnchors, noOfsaEdges) ...
    - sparse(noOfSensors+AnchorIndex, 1:noOfsaEdges, 1, noOfSensors+noOfAnchors, noOfsaEdges);
S = [S1, S2];
ssaDistance = [ssDistance; saDistance];


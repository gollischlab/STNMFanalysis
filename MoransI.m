% This function is part of the STNMFanalysis package and computes the Moran's I for a matrix.

%%%%
% Copyright (C) 2017 Tim Gollisch
%
% This file is part of STNMFanalysis.
% 
% STNMFanalysis is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% STNMFanalysis is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with STNMFanalysis.  If not, see <http://www.gnu.org/licenses/>.
%%%%

function [ output ] = MoransI( A, Nx, Ny )
% Function to compute Moran's I (measure of spatial autocorrelation) for a real-valued 2D-matrix A with dimensions Nx x Ny
% For the formula, see, e.g., https://en.wikipedia.org/wiki/Moran's_I
% The applied weight matrix gives weight 1 to pairs of horizontally or vertically neighboring pixels and weight 0 to all other pixel pairs.
    B = A - mean(mean(A));  % subtracts the total mean over the entire input matrix
    sumNeighborProducts = 0;
    sumWeights = 0;
    sumSquaredDeviations = 0;
    for x=1:Nx
        for y=1:Ny
            if x>1
                sumNeighborProducts = sumNeighborProducts + B(x,y)*B(x-1,y);
                sumWeights = sumWeights + 1;
            end;
            if x<Nx
                sumNeighborProducts = sumNeighborProducts + B(x,y)*B(x+1,y);
                sumWeights = sumWeights + 1;
            end;
            if y>1
                sumNeighborProducts = sumNeighborProducts + B(x,y)*B(x,y-1);
                sumWeights = sumWeights + 1;
            end;
            if y<Ny
                sumNeighborProducts = sumNeighborProducts + B(x,y)*B(x,y+1);
                sumWeights = sumWeights + 1;
            end;
            sumSquaredDeviations = sumSquaredDeviations + B(x,y)*B(x,y);
        end;
    end;
    output = Nx*Ny*sumNeighborProducts/sumWeights/sumSquaredDeviations;
end


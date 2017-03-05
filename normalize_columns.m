% This function is part of the STNMFanalysis package and returns a matrix with normalized columns.

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

function B = normalize_columns( A )
% function to normalize every column of the matrix A to unit Euclidean norm and return as a matrix B
[Nrows,Ncolumns] = size(A);
C=1./sqrt(sum(A.*A,1));
B = A .* repmat( C, Nrows, 1 );

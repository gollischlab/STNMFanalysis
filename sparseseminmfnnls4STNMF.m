% This function is part of the STNMFanalysis package and performs the sparse-semi-NMF.
% It is modified from the sparseseminmfnnls.m function by Yifeng Li, which
% is available under the GNU General Public License (see below).
% Modifications by Tim Gollisch and Jian K. Liu

function [A,Y,finalResidual,numIter,tElapsed]=sparseseminmfnnls4STNMF(X,Ystart,k,Nx,Ny,option)
% Modified from the original code in the NMF MATLAB Toolbox
% Mainly, we provide the initialization of the modules through an argument
% that we pass to the function (Ystart) instead of drawing random numbers
% for initialization of the modules here.
% Also, we plot intermediate results and pass the spatial dimensions Nx and
% Ny as arguments to this function.
% And then, we simplify the code slightly in ways that do not impede our application of the code.
%
% Below is the description of the original function:
%
% ===============================================
% Sparse Semi-NMF based on NNLS: X=AY, s.t. Y>0.
% Definition:
%     [A,Y,numIter,tElapsed,finalResidual]=seminmfnnls(X,k)
%     [A,Y,numIter,tElapsed,finalResidual]=seminmfnnls(X,k,option)
% X: matrix of mixed signs, dataset to factorize, each column is a sample, and each row is a feature.
% k: scalar, number of clusters.
% option: struct:
% option.iter: max number of interations. The default is 1000.
% option.dis: boolen scalar, It could be 
%     false: not display information,
%     true: display (default).
% option.residual: the threshold of the fitting residual to terminate. 
%     If the ||X-XfitThis||<=option.residual, then halt. The default is 1e-4.
% option.tof: if ||XfitPrevious-XfitThis||<=option.tof, then halt. The default is 1e-4.
% A: matrix, the basis matrix.
% Y: matrix, the coefficient matrix.
% numIter: scalar, the number of iterations.
% tElapsed: scalar, the computing time used.
% finalResidual: scalar, the fitting residual.
% References:
% ...

%%%%
% Copyright (C) <2012>  <Yifeng Li>
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
% 
% Contact Information:
% Yifeng Li
% University of Windsor
% li11112c@uwindsor.ca; yifeng.li.cn@gmail.com
% May 01, 2011
%%%%

global residualArray;
global iterationsArray;
global iterationsCounter;
global hResidualsFig;
global hModuleSearchFig;
global transformation_name;

tStart=tic;
[r,c]=size(X); % c is # of samples, r is # of features
Y=Ystart; %Y=rand(k,c);
XfitPrevious=Inf;%(size(X));
for i=1:option.iter
    iterationsCounter = iterationsCounter + 1;
    
    A=X*pseudoinverse(Y,2^32);  % Finds the optimal set of weights for the current modules by multiplication with the pseudoinverse of the module matrix
    A=normalize_columns(A);     % Normalize each column to unit Euclidean norm
    Ae=[A;sqrt(option.beta)*ones(1,k)]; % Use an extended matrix of weights to integrate the sparsity constraint
    X0=[X;zeros(1,c)];                  % Also use extended data matrix to incorporate the sparsity constraint
    Y=kfcnnls4STNMF(Ae,X0);             % Call the function for solving the non-negative least-squares problem for optimizing the modules for fixed weights
    
    % Every once in a while, evaluate and plot current solution:
    if( i==1 || mod(i,10)==0 || i==option.iter )
        if option.dis
            disp(['Iterating >>>>>> ', num2str(i),'th']);
        end
        XfitThis=A*Y;
        fitRes=sqrt(sum(sum((XfitPrevious-XfitThis).^2))); % Computes the improvement in the fit
        XfitPrevious=XfitThis;
        curRes=norm(X-XfitThis,'fro');

        % Here some new inserted code for plotting intermediate modules and Moran's I values        
        
        % Plot evolution of residuals
        if( iterationsCounter > 1 )
            residualArray = [residualArray curRes];
            iterationsArray = [iterationsArray iterationsCounter];
        end;
        unit = cell(1,k);
        for kk=1:k
            unit{kk} = reshape(Y(kk,:),[Nx,Ny]);
        end;
        figure( hResidualsFig );
        subplot(3,1,1);
        plot(iterationsArray, residualArray);
        title('Evolution of total residual');

        % Plot current values of Moran's I
        subplot(3,1,2);
        for kk=1:k
            moransI(kk) = MoransI(unit{kk},Nx,Ny);
        end;
        plot( moransI, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'y' );
        title('Spatial autocorrelations of current modules');
        ylim([-0.5 1]);
            
        % Plot current set of modules
        figure( hModuleSearchFig );
        NperRow = ceil(k/2);  % number of subplots per row
        for kk=1:k
            subplot('Position', [mod((kk-1),NperRow)/NperRow+0.05/NperRow, 0.47*(1-floor((kk-1)/NperRow))+0.01, 0.9/NperRow, 0.45] );
            imagesc( unit{kk} );
            set(gca,'xtick',[]);
            set(gca,'ytick',[]);
            axis off;
            if(kk==ceil(k/4)) title( transformation_name ); end;
        end;
        drawnow;
        % End of code for plotting intermediate data        
        
        if option.tof>=fitRes || option.residual>=curRes || i==option.iter
            numIter=i;
            finalResidual=curRes;
            break;
        end
    end
end
tElapsed=toc(tStart);
end

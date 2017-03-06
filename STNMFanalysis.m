
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STNMFanalysis: this code performs semi-non-negative-matrix factorization (semi-NMF)
% on the spike-triggered stimulus ensemble
%
% This accompanies the manuscript be Liu et al.:
% "Inference of neuronal functional circuitry with spike-triggered non-negative matrix factorization".
%
% The spike-triggered stimulus ensemble is assumed to correspond to 2-D spatial stimuli and is read in from file.
% Time (if the stimulus was spatiotemporal) is considered to be integrated out, e.g., by computing a weighted sum over
% the spatiotemporal spike-triggered stimuli, where the weights may come from the cell's temporal STA.
%
% The semi-NMF analysis is based on functions adopted from the NMF MATLAB Toolbox
% by Yifen Li and Alioune Ngom (https://sites.google.com/site/nmftool/).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%
% Copyright (C) 2017 Tim Gollisch and Jian K. Liu
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


clear;
close all;

% Some relevant parameter for the analysis
K = 20;                   % number of modules considered in the analysis
Niter = 200;              % number of iterations between perturbations of the current best set of modules
Npert = 200;              % number of perturbations to be tried out
threshold_MoransI = 0.25; % threshold on Moran's I; above this, a module is considered localized
% These parameters are likely the most important ones for adjusting the
% analysis to the data at hand.
% K should be somewhat larger than the expected number of subunits.
% For Niter and Npert, obviously, the longer the better, but this is a
% matter of computational efficiency. 20 iterations may be sufficient if
% one goes for multiple runs with different initial conditions. More
% iterations and more perturbations give more robust results if one wants to
% rely on a single run.

% Specify whether final set of modules should be saved to disc
% Set to 1 for saving and to zero for not saving
% See end of code for destination of saved data
saveFinalResult = 0;

% Load data.
load('model_cell.mat');
% The data file should contain three variables:
% STE: Nspikes x Nstimdimensions matrix of the spike-triggered stimulus ensemble
% Nx and Ny: stimulus dimensions in space; this must be set so that Nx*Ny=Nstimdimensions
% Each spike-triggered stimulus is then interpreted as an Nx x Ny spatial stimulus

STE = STE/std(reshape(STE,[],1)); % Normalizing the spike-triggered ensemble to unit STD; not strictly necessary, but may be helpful.

% some global variables to facilitate plotting intermediate results
global residualArray;       % an array that holds the development of the total residual
global iterationsArray;     % an array that hold the iterations corresponding to the residualArray
global iterationsCounter;   % counts the total number of iterations
global hResidualsFig;       % figure handle
global hModuleSearchFig;    % figure handle
global transformation_name; % a string variable that holds the name of the latest perturbation for display in the figure

% preparing figures and their placement on the screen
hSTAfig = figure('Position', [50, 550, 400, 400]);  % for showing the STA
hResidualsFig = figure('Position', [50, 50, 400, 400]);   % for showing evoluation of residuals and Moran's I values
hModuleSearchFig = figure('Position', [500, 50, 1000, 400]); % for showing current state of module search
hCurrentSolutionFig = figure('Position', [500, 550, 1000, 400]); % for showing current best modules

% Plot STA
STA = reshape( mean( STE ), [Nx,Ny] );
figure( hSTAfig );
imagesc( STA );
title( 'STA' );
drawnow;

% Initial parameters for semi-NMF
option.iter = Niter; % setting the number of iterations between perturbations
option.beta = 0.1;   % sparseness parameter
option.dis = 0;      % 1 = periodic display of number of iterations; 0 = no display of number of iterations
option.tof = 1e-30;
option.residual = 1e-30;
% tof and residual are halting conditions on the latest change between the
% fit and on the total residual.
% They are set here so low that they essentially never apply, so we just
% let the algorithm run to the max number of iterations.

% Preparing the arrays for the residuals and iterations
residualArray = [];
iterationsArray = [];
iterationsCounter = 0;

% Initializing the modules matrix with random numbers between zero and unity
Ystart=rand( K, Nx*Ny );

% A cell array for holding the modules
unit = cell(1,K);

% Array of values of Morans'I for the modules
moransI = zeros(1,K);

transformation_name = ' ';

for q=1:Npert % The main loop, each time with a different perturbation of the modules (set after each cycle).
    
    % Call to the actual semi-NMF function, modified from the NMF MATLAB Toolbox.
    [weights,modules,fRes] = sparseseminmfnnls4STNMF(STE, Ystart, K, Nx, Ny, option);

    % If not the first run, evaluate whether the latest perturbation has led to an improvement
    if( q>1 && fRes>previous_fRes )  % No improvement
        modules = previous_modules;  % Revert to the previous best set of modules
        fRes = previous_fRes;        % Revert to the previous best value of total residuals
    else % Improvement (or first run)
        % Extract 2D-pattern of modules into cell "unit" from modules matrix
        % and compute Moran's I for each module
        for k=1:K
            unit{k} = reshape(modules(k,:),[Nx,Ny]);
            moransI(k) = MoransI(unit{k},Nx,Ny);
        end;
                
        % Sort modules inside modules matrix according to spatial
        % autocorrelation so that the most localized ones come first
        modulesbuffer = modules;
        [sortedMoransI, indices] = sort( moransI, 'descend' );
        for k=1:K
            modules(k,:) = modulesbuffer(indices(k),:);
            unit{k} = reshape(modules(k,:),[Nx,Ny]);
            moransI(k) = MoransI(unit{k},Nx,Ny);
        end;
        
        % Save into buffer for current best solution
        previous_fRes = fRes;
        previous_modules = modules;
    end;
    
    % plot current best set of modules
    figure( hCurrentSolutionFig );
    NperRow = ceil(K/2);  % number of subplots per row
    for k=1:K
        subplot('Position', [mod((k-1),NperRow)/NperRow+0.05/NperRow, 0.47*(1-floor((k-1)/NperRow))+0.01, 0.9/NperRow, 0.45] );
        imagesc( unit{k} );
        set(gca,'xtick',[]);
        set(gca,'ytick',[]);
        axis off;
        if(k==ceil(K/4)) title( 'Current best' ); end;
    end;
    drawnow;
    
    % plot Moran's I values for current best set of modules
    figure( hResidualsFig );
    subplot(3,1,3);
    plot( moransI, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r' );
    title('Spatial autocorrelations of currently best modules');
    ylim([-0.5 1]);

    % Perturbe the modules by performing a manipulation on a random module

    % First find the localized modules by comparing the Moran's I values to a threshold
    Igood = find(moransI>threshold_MoransI);
    Ngood = length( Igood );

    % Before performing perturbation, initialize the starting condition for
    % the next run with the current best set of modules
    Ystart = modules;
    
    perturbationmode = randi(5); % Draw at random, which perturbation is to be performed
    if( perturbationmode == 1 )
        % delete one localized module
        transformation_name = 'Latest perturbation: one localized module replaced with noise';
        if( Ngood>=1 )
            ix1 = randi(Ngood);  % select one localized module at random
            Ystart(ix1,:) = rand(1,Nx*Ny);
        end;
    elseif( perturbationmode == 2 )
        % duplicate a module and then add noise to both copies
        transformation_name = 'Latest perturbation: one localized module duplicated + noise';
        if( Ngood>=1 && Ngood<K )
            ix1 = randi(Ngood);  % select one localized module at random
            ix2 = randi(K-Ngood)+Ngood;  % select one non-localized module at random
            Ystart(ix2,:) = Ystart(ix1,:);
            Ystart(ix1,:) = Ystart(ix1,:) + rand(1,Nx*Ny)*sum(Ystart(ix1,:))/(Nx*Ny)*2;
            Ystart(ix2,:) = Ystart(ix2,:) + rand(1,Nx*Ny)*sum(Ystart(ix2,:))/(Nx*Ny)*2;
        end;
    elseif( perturbationmode == 3 )
        % split one localized module into halfs vertically
        transformation_name = 'Latest perturbation: one localized module split vertically';
        if( Ngood>=1 && Ngood<K  )
            ix1 = randi(Ngood);  % select one localized module at random
            ix2 = randi(K-Ngood)+Ngood;  % select one non-localized module at random
            temp = Ystart(ix1,:);
            [~,IX] = max(temp);
            Half_1 = zeros(1,Nx*Ny);
            Half_2 = zeros(1,Nx*Ny);
            Half_1(1:IX) = temp(1:IX);
            Half_2(1+IX:end) = temp(1+IX:end);
            Ystart(ix1,:) = Half_1;
            Ystart(ix2,:) = Half_2;
        end;
    elseif( perturbationmode == 4 )
        % split one localized module into halfs horizontally
        transformation_name = 'Latest perturbation: one localized module split horizontally';
        if( Ngood>=1 && Ngood<K  )
            ix1 = randi(Ngood);  % select one localized module at random
            ix2 = randi(K-Ngood)+Ngood;  % select one non-localized module at random
            temp = Ystart(ix1,:);
            temp = reshape(temp,Nx,Ny);
            temp = temp';
            [~,IX] = max(temp(:));
            Half_1 = zeros(1,Nx*Ny);
            Half_2 = zeros(1,Nx*Ny);
            Half_1(1:IX) = temp(1:IX);
            Half_2(1+IX:end) = temp(1+IX:end);
            Ystart(ix1,:) = reshape(reshape(Half_1',Ny,[])',1,[]);
            Ystart(ix2,:) = reshape(reshape(Half_2',Ny,[])',1,[]);
        end;
    else
        % re-seed non-localized modules
        transformation_name = 'Latest perturbation: non-localized modules re-seeded';
        if( Ngood<K )
            Ystart(Ngood+1:K,:) = rand(K-Ngood,Nx*Ny);
        end;
    end;
end;

% Save final modules and Moran's I values if so specified
if( saveFinalResult )
    save('final_modules_from_STNMF.mat','unit','moransI');
end;

disp( 'STNMF analysis finished!' );




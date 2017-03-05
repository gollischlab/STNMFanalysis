# STNMFanalysis
Code for neuronal subunit analysis with spike-triggered non-negative matrix factorization

Accompanying the manuscript
“Inference of neuronal functional circuitry with spike-triggered non-negative matrix factorization”
by Liu, Schreyer, Onken, Rozenblit, Khani, Krishnamoorthy, Panzeri, and Gollisch

Contact: Tim Gollisch, tim.gollisch@med.uni-goettingen.de
Website: http://www.retina.uni-goettingen.de/

### Description:
This code performs spike-triggered non-negative matrix factorization (STNMF) on a supplied set of spike-triggered stimuli. The stimuli are interpreted as being purely spatial, organized in a rectangular 2D layout, with Nx x Ny pixels.<br>
The algorithm performs a decomposition of this spike-triggered ensemble (STE) into a number K of non-negative modules, so that the spike-triggered stimuli are best approximated by linear combinations of the modules in a least-squares sense.<br>
The STNMF is based on a semi-NMF approach, that is, the stimulus elements can be positive or negative, and only the elements of the modules are constrained to be non-negative.<br>
The optimization of the modules occurs iteratively, and the number of applied iterations is specified in the program. A key element of the algorithm is that perturbations of the current best modules are inserted periodically after fixed numbers of iterations. If the perturbation does not yield improved modules (in terms of the reconstruction residual), it is discarded, the previous best modules are retained, and a new perturbation is tested. See the related manuscript for details.

### Usage:
Download the files in the repository into a single directory or into your Matlab path. To run, execute STNMFanalysis.m.

### Output:
The program generates 4 figures:<br>
1. The STA, displayed as a 2D color map.<br>
2. The evolution of the total residual, plotted versus the total number of performed iterations, as well as the set of spatial autocorrelation values (Moran’s I) for the currently investigated modules and for the current best set of modules.<br>
3. The set of currently investigated modules, displayed as 2D color maps. Here, the title of the figure indicates which perturbation was applied most recently to the modules (shown only after the first perturbation).<br>
4. The current best set of modules, displayed as 2D color maps, updated after every perturbation (also shown only after the first perturbation).

### Input/data:
The input to the analysis is the spike-triggered ensemble (STE). This is read from file in the beginning of STNMFanalysis.m. In the distributed version, the program loads data from the file model_cell.mat. A sample data file model_cell.mat is distributed with the code, so the program can just be run directly as is, using the supplied data as input.<br>
The data file contains the STE from a model simulation (5 overlapping subunits with threshold-quadratic subunit nonlinearities as in Fig. 2 of the corresponding manuscript; about 3000 simulated spikes).<br>
Loading the data file creates the following variables, which are later expected to exist by the STNMFanalysis program:<br>
Nx and Ny: integer numbers that correspond to the spatial dimensions of the spike-triggered stimuli.<br>
STE: an Nspikes x Nstimdimensions matrix that contains in each row the stimulus elements p for a spike-triggered stimulus. Each spike-triggered stimulus is given as a 1D vector, which concatenates all the sequences of elements in the y-dimensions: p(1,1), p(1,2), …, p(1,Ny), p(2,1), …, p(2,Ny), …, p(Nx,Ny).<br>
Nstimdimensions = Nx x Ny must be fulfilled!<br>
For running the program on your own data, provide a data file in the same format for generating the variables Nx, Ny, and STE. Then read this data file instead of model_cell.mat in STNMFanalysis.m.

### Hints for exploring the program:
A full run of the program may be quite long (many hours), but often, a good idea of subunits can be obtained already after 2-3 perturbations. This depends on the randomly chosen initializations, so it may be instructive to start and stop the program a few times to see different developments of the modules. That is, let the program run for a few perturbations only (typically a few hundred iterations), then stop (Ctrl+C under Windows or Linux; Command+. on a Mac) and start again, which will provide a new, random initialization for the search of subunits.<br>
For more detailed investigations of new data sets and for balancing computational speed versus robustness of results, it makes sense to think about or play around with the main parameters of the algorithm, which are set at the beginning of STNMFanalysis.m.<br>
These main parameters are K (the number of considered modules) and Niter (the number of iterations before a new perturbation is tested).<br>
You can also change Npert (the total number of perturbations to be tried out), which affects essentially how long the search continues before the program stops, and threshold_MoransI, which determines which modules are treated in the applied perturbations as localized subunits and which are treated as likely noise.

### License issues:
The code in this package is distributed under the GNU General Public License.<br>
The semi-NMF algorithm is based on modified code from the NMF MATLAB Toolbox by Yifeng Li, Alioune Ngom (https://sites.google.com/site/nmftool/, citation: Y. Li and A. Ngom, The non-negative matrix factorization toolbox for biological data mining, BMC Source Code for Biology and Medicine, vol 8, pp. 10), which is also distributed under the GNU General Public License. See accompanying license file and license notes in the source code files.

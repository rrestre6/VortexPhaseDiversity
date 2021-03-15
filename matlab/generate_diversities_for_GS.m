%%
% File: generate_diversities_for_GS.m
% Author: Santiago Echeverri and Carlos Cuartas
% Date:	  18/07/2016
%%
%	This program is used to generate phase masks to be used as phase
%	diversities in a GS PD scheme with multiple images.
% The masks are a combination of Zernike parametrized aberrations $k$
% and spiral phase masks of topological charge $l$.
% The aberrations are organized as follow (Wyant expansion)
% l = [defocus, astigmatism, astigmatism45°, comaX, comaY, spherical, trefoil, trefoil45°, ...
% secondAstigmatism, secondAstigmatism45°, secondComaX, secondComaY, ...],
% piston and tilts has been eliminated in the numeration
%
%%
%	IN
%		RVect : 1D matrix with radial coordinates within the circle
%                 to use Zernike.m
%		TVect:  1D matrix with azimuthal coordinates within the
%                 circle to use Zernike.m
%		mask :  120x120 matrix where mask(R < 1) = 1; mask(R > 1) = 0;
%		mask2 : 120x120 matrix where mask(R < rLimit) = 1; mask(R > rLimit) = NaN;
%		basicSpirals: Cell that contains 120x120 matrices that define phase spirals
%                         of OAM's -2,-1,0,1,2
%		L : Array with list of spiral diversity indices for each mask
%		K : Cell with arrays of aberration coefficients for each
%          Zernike diversity mask
%%
%  OUT
%		Diversities : Cell that contains 120x120 matrices that define phase masks
%                            with diversities
% 	AntiDiversities: Cell that contains 120x120 matrices that define phase masks
%															with the antidiversities
% octave-3.8.1.
%%
%
% Example:
% 	In this example six diversities are generated, the desired topological charges are 0, 1 and 2 and each one requires a diversty of 0.5lambda defocus
%
%		% Constants
%		warning('off')
%		N = 2; gridSize = 256;
%		dx = N / gridSize;
%		rLimit = 1;
%
%		% Mesh and mask
%		[X, Y] = meshgrid(-N / 2: dx: N/2 - dx);
%		[Theta, R] = cart2pol(X, Y);
%
%		mask = zeros(gridSize);
%		mask(R < rLimit) = 1;
%
%		mask2 = zeros(gridSize);
%		mask2(R > rLimit) = NaN;
%
%		rPupil = R / rLimit;
%		RVect = rPupil(~isnan(mask2));
%		TVect = Theta(~isnan(mask2));
%
%		basicSpirals = cell(5,1);
%		l = 1;
%		for order = [-2, -1, 0, 1, 2]
%		  basicSpirals(l) = spiralGenTo(ones(size(X)), order);
%		  l = l + 1;
%		end
%
%		% Diversities
%		lDiv = [0 0 1 1 2 2];
%		kDiv = {[0 0 0]', [0.5 0 0]', [0 0 0]', [0.5 0 0]', [0 0 0]', [0.5 0 0]'};
%
%		[Diversities, AntiDiversities] = generate_diversities_for_GS(RVect,TVect, mask,mask2,rLimit,basicSpirals,lDiv,kDiv);
%
%%

%% Function definition

function [Diversities, AntiDiversities] = generate_diversities_for_GS(RVect,TVect, ...
                                                  mask,mask2,rLimit,basicSpirals,L,K)

	%% Constants
	nZernikes  = length(K);
	gridSize = size(mask2,1);
	dA = (2 / gridSize / rLimit) ^ 2;

	%% Preload matrices
	Diversities = cell(nZernikes,1);
	AntiDiversities = cell(nZernikes,1);

	for lk = 1:nZernikes

		%% Get spiral mask of order l
		% topological charges are organized as basicSpirals = l{-2, -1, 0, +1, +2}
		spiral_mask = basicSpirals{ L(lk) + 3 };

		%% Create k Zernike diversities
		coeff=3:size(K{lk},1)+2;
		[waveAberationVector,~,~] = Zernike (coeff, RVect, TVect);

		% Normilize Zernike polynomials
		zProds = waveAberationVector.'* waveAberationVector * dA;
		Zernike_mask_1D = bsxfun(@rdivide, waveAberationVector, sqrt(diag(zProds).'));
		waveAberationVector = sum(bsxfun(@times, K{lk}.', Zernike_mask_1D), 2);

		% Create 2D Zernike polynomials
		Zernike_mask_2D = zeros(gridSize);
		Zernike_mask_2D(~isnan(mask2))=(waveAberationVector);

		%% Final diversity composed by the l and k diversities, and their respective antidiversity
		Diversities{lk} = spiral_mask .* exp(1i * Zernike_mask_2D);
		AntiDiversities{lk} = basicSpirals{ 3 - L(lk) } .* exp(-1i * Zernike_mask_2D);
    end
end

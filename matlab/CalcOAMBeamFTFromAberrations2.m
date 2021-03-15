%%
% File: CalcOAMBeamFTFromAberrations2.m
% Author: Santiago Echeverri, René Restrepo, Carlos Cuartas and Néstor
% Uribe
% Date:	  19/07/2016
% Modification:	  10/07/2020
% Vortex-enhanced coherent-illumination phase diversity for phase retrieval in coherent imaging systems 
% Optics Letters Vol. 41, Issue 8, pp. 1817-1820 (2016). https://doi.org/10.1364/OL.41.001817
%%
%  This program creates an optical vortex aberrated by using Zernike basis with a weight vector in a 4-f setup.
%
%   Obj  --- fft ---- aberrations+SLM ---- ifft ----- Ima
%
%   Obj: Entrace field.
%   fft: Fast Fourier transform of the entrace field.
%   aberration+SLM: In the Fourier plane aberrations and SLM is added.
%   ifft: Inverse fast Fourier transform of the complex field.
%   Ima: Image plane which contains the vortex and the aberrations.
%%
%  IN
%	FTSize: Size after creating a padarray with before fft's are computed.
%	gridSize: Pixels used to create grids.
%	gaussianC: Gaussian distribution waist.
%	rLimit: Maximum radius in Fourier field.
%	zernikeAberrations: Vector weights on Zernike basis.
%	dx: Vortex displacement in x direction.
%	dy: Vortex displacement in y direction.
%	order: Vortex topological charge.
%
%  OUT
%	amplitudeCenter: Final amplitude.
%	IntensityCenter: Final intensity (|A^2|).
%	ampIni: Initial amplitude at object plane.
%	waveAberation: Surface with final aberrations.
%	gridSize: Simulated images size.
%	colRangeSLM: Column range at SLM.
%	rowRangeSLM: Row range at SLM.
%
% octave-3.8.1.
%%

%% Function definition

function [amplitudeCenter, IntensityCenter, ampIni, SLM, gridSize, colRangeSLM, rowRangeSLM] = CalcOAMBeamFTFromAberrations2(FTSize, gridSize, gaussianC, rLimit, zernikeAberrations, dx, dy, order)

  warning('off');

	%% Check if aberrations is a row vector
	if size(zernikeAberrations, 2) > size(zernikeAberrations, 1)
		zernikeAberrations = zernikeAberrations';
	end

	%% Constants
	dx = 0; % gaussian beam displacement from center in X
	dy = 0; % gaussian beam displacement from center in Y
	[X, Y] =  meshgrid(-1:(2 / gridSize):(1 - 2 / gridSize), -1:(2 / gridSize):(1 - 2 / gridSize));
	X2 = X - dx;
	Y2 = Y - dy;
	[Theta, R] = cart2pol(X, Y); % pupil coordinates
	[Theta2, R2] = cart2pol(X2, Y2); % Vortex coordinates

	%% Create masks
	mask = zeros(gridSize); % Pupil mask
	mask(R < 1) = 1;
	maskp = zeros(gridSize); % Phase mask
	maskp(R < rLimit) = 1;
	mask2=zeros(gridSize); %Create a new pupil with in order to use Zernike polynomials
	mask2(R>rLimit)=NaN;
	rPupil = R / rLimit;
	RVect=rPupil(~isnan(mask2)); %Coordinates to use Zernike.m
	TVect=Theta(~isnan(mask2));

	%% Image Range Calculus
	colRangeSLM = (FTSize - gridSize) / 2 + 1:(FTSize + gridSize) / 2;
	rowRangeSLM = (FTSize - gridSize) / 2 + 1:(FTSize + gridSize) / 2;

	%% Entrace field
	%The aberrations are with respect to pupil (X, Y) and only the vortex with respect to X2, Y2.
	gaussAmp = exp(- (X2 .^ 2 + Y .^ 2) / gaussianC ^ 2);
	%Gaussian distribution is the initial simulated amplitude

	Uin = mask .* gaussAmp .* exp(-1i.*(zeros(gridSize))); %Entrace field
	UinPadded = padarray(Uin, ([FTSize FTSize] - gridSize) / 2, 0, 'both'); % Entrace field padded
	ampIni=UinPadded;
  ampIni = NaN;

	%% FFT of entrace field
	Uout= fftshift(fft2(fftshift(UinPadded)));

	%% Aberrations
	rPupil = R / rLimit; %Scaling of R in order to have a normalized aperture.
  dA = (2 / gridSize / rLimit) ^ 2;
  % This would be originally dA = dx * dy, with dx = dy = 2 / gridSize as
	% defined in the meshgrid above, but because of the change in rPupil needs
	% to be scaled too

	% 1D matrix with coordinates within the circle to use Zernike.m
	RVect=rPupil(~isnan(mask2));
	TVect=Theta(~isnan(mask2));
  coeff = 3:size(zernikeAberrations,1)+2; % Vector with weights
  [waveAberationVector,~,~]=Zernike(coeff,RVect,TVect); %Calculate Zernike aberrations
	% The interior product of all Zernikes with themselves is
	zProds = waveAberationVector.'* waveAberationVector * dA;
	% Looking at theis product we see that the polynomials retrieved by
	% "Zernike.m" are orthogonal, BUT NOT NORMALIZED. Checking
	% wikipedia,
	% https://en.wikipedia.org/wiki/Zernike_polynomials#Orthogonality, I
	% suspect we are missing a pi / (2n +2). Or just use this numerical
	% integration to normalize them:
	zNorm = bsxfun(@rdivide, waveAberationVector, sqrt(diag(zProds).'));
	% It is easy to see that there are normalization factors by looking at the
	% normalized polynomials here:
	% https://en.wikipedia.org/wiki/Zernike_polynomials#Zernike_polynomials
	% (pay attention to those radicals). However these in the wikipedia are
	% normalized to pi, so all that list needs to be divided by sqrt(pi).
	waveAberationVector = sum(bsxfun(@times, zernikeAberrations.', zNorm), 2);
  waveAberation = zeros(gridSize); %Reorganize Z. vects in a matrix
  waveAberation(~isnan(mask2))=(waveAberationVector); %Create aberrated pupil for the system waveAberation = wrapToPi(waveAberation);

  %% Spiral mask
  spiralMask=spiralGenTo(ones(size(X2)),order);
  SLM = spiralMask .* exp(1i * waveAberation);% aberrations + SLM
  SLMPadded=padarray(SLM, ([FTSize FTSize] - gridSize) / 2, 1, 'both');

  %% Fourier plane

  Uin2=Uout.*SLMPadded;

  %%Range Calculus

  colRangeSLM = (FTSize - gridSize) / 2 + 1:(FTSize + gridSize) / 2;
  rowRangeSLM = (FTSize - gridSize) / 2 + 1:(FTSize + gridSize) / 2;

  %% IFFT of complex field
  Uout2 = ifftshift(ifft2(ifftshift(Uin2)));

  %% Image field and final images

  amplitudeCenter=abs(Uout2); %Amplitude
  phaseCenter=angle(Uout2); %Phase
  IntensityCenter = abs(Uout2).^2;
  IntensityCenter=IntensityCenter-min(min(IntensityCenter));
  IntensityCenter=(IntensityCenter/max(max(IntensityCenter))); %Normalized intensity
end

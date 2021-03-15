%%
% File:	  RetrievePhaseFromPhaseDiversity1.m
% Author: Santiago Echeverri, Ren� Restrepo, Carlos Cuartas and N�stor
% Uribe
% Date:	  19/07/2016
% Modification:	  10/07/2020
% Vortex-enhanced coherent-illumination phase diversity for phase retrieval in coherent imaging systems 
% Optics Letters Vol. 41, Issue 8, pp. 1817-1820 (2016). https://doi.org/10.1364/OL.41.001817
%%
%  This program retrieves phase by using phase diversity algorithm procedure.
%%
%  IN
%	beamImage: Array with PD images.
%	FTSize: Size after creating a padarray with before fft's are computed.
%	beamDiameter: Pixels used to create grids.
%	gaussianC: Gaussian distribution waist.
%	rLimit: Maximum radius in Fourier field.
%	baseZAberrations: Vector array with base aberrations for PDs.
%	baseOAM: OAM guess (if unknow, otherwise = 0).
%	dx: Vortex displacement in x direction.
%	dy: Vortex displacement in y direction.
%	ZPD: Array with zernikes in each P.D.
%	OAMPD:  Vector with vortices topological charges.
%	minType: Minimization type.
%	minParms: Minimization parameter: {ftol,isz}; ftol: stop when result does not improve. isz: initial step size (expected distance from minimum).
%	verbosity: Verbosity level.
%%
%  OUT
%	zernikeAberrations: Final zernike aberrations.
%	ExpIntensity: Experimental images array.
%	colRangeIn: Column range at the entrace field.
%	rowRangeIn: Row range at the entrace field.
%	colRangeSLM: Column range at SLM.
%	rowRangeSLM: Row range at SLM.
%	dx: Vortex displacement in x direction.
%	dy: Vortex displacement in y direction.
%
% octave-3.8.1.
%%

  %% Function definition

function [zernikeAberrations, ExpIntensity, colRangeIn, rowRangeIn, ...
          colRangeSLM, rowRangeSLM, dx, dy, waveAberation, ...
          IntensityCenter] = ...
        RetrievePhaseFromPhaseDiversity1(beamImage, FTSize, beamDiameter, ...
                                         gaussianC, rLimit, ...
                                         baseZAberrations, baseOAM, ...
                                         dx, dy, ZPD, OAMPD, minType, ...
                                         minParms, verbosity) 
    
 % warning ("off"); %Get rid of multiple messages
 % pkg load optim;  %Load the required packages
 % pkg load image;

  %% Experimental images normalization
  
  nPDs = size(beamImage, 3); % Number of PD images
  
  for k = 1:nPDs
  beamImage(:,:,k) = beamImage(:,:,k) - min(min(beamImage(:,:,k)));
  beamImage(:,:,k) = (beamImage(:,:,k) / max(max(beamImage(:,:,k))));
  end

  %% Constants
  
  imageSize = size(beamImage, 1); % Input images size
  gridSize = beamDiameter; %Beam diameter used to make later grids
  colRangeIn = (FTSize - imageSize) / 2 + 1:(FTSize + imageSize) / 2; % Array size based in Fourier transform size and images size
  rowRangeIn = (FTSize - imageSize) / 2 + 1:(FTSize + imageSize) / 2;
  padSize=([FTSize FTSize 0] - [imageSize imageSize 0]) / 2;%Padarray size
  beamImagePadded = padarray(beamImage, padSize, 0, 'both');
  ExpIntensity = beamImagePadded; %Experimental images padded
  clear beamImagePadded; 

  %% First iteration
  
  [amplitudeCenter, IntensityCenter, ~, waveAberation, gridSize, colRangeSLM, rowRangeSLM] = CalcOAMBeamPhaseDiversityFTFromAberrations1(FTSize, beamDiameter, gaussianC, rLimit, baseZAberrations, baseOAM, dx, dy, ZPD, OAMPD, nPDs);
  
%% Create handles and initial figures

  maxExp = max(ExpIntensity(:));
  currentImage = 1;
  figure(currentImage);
  imagesc(abs(IntensityCenter(colRangeSLM, rowRangeSLM,1))),title('ignore'); colorbar; axis image;
  title('Simulated beam amplitude at SLM');
  currentImage=currentImage+1;
  for k = 1:nPDs
    figure(currentImage);
    hExpFT(k) = imagesc(ExpIntensity(colRangeIn, rowRangeIn, k), [0 maxExp]); colorbar; axis image;
    title(sprintf('Experimental FT %d', k));
    currentImage=currentImage+1;
    figure(currentImage);
    hAbsFT(k) = imagesc(abs(IntensityCenter(colRangeIn, rowRangeIn, k)), [0 maxExp]); colorbar; axis image;
    title(sprintf('Initial simulated FT - intensity %d', k));
    currentImage=currentImage+1;
  end
  figure(currentImage), hPhaseIFT = imagesc(angle(waveAberation(:,:,1))); colorbar; axis image; title('Initial simulated beam phase at SLM - PD 1');
drawnow;
currentImage=currentImage+1;

  %iterate = input("Enter to start fitting procedure..."); %Start fitting procedure

  %% Check OAM before minimization
  if (baseOAM == 0 && ~any(OAMPD))
    fitDxDy = false;
    x = baseZAberrations;
  else
    fitDxDy = true;
    % Get which dimension is the one that contains the data
    dim = 1 + (size(baseZAberrations, 2) > size(baseZAberrations, 1));
    % And concatenate in that dimension
    x = cat(dim, baseZAberrations, dx, dy);
  end

  %% Minimization

 if (strcmp(minType, 'nelder'))
    %[x_m, v_m, nevals] = minimize(@(x) ...
    %                             FTBeamPhaseDiversityErrorFunction1(FTSize, beamDiameter, gaussianC, rLimit, x, baseOAM, ZPD, OAMPD, nPDs, ExpIntensity, colRangeIn, rowRangeIn, hAbsFT, hPhaseIFT, fitDxDy, verbosity), x, "isz", minParms{1}, "verbose", "ftol", minParms{2});
    [x_m, v_m, nevals] = fminsearch(@(x) ... 
                                  FTBeamPhaseDiversityErrorFunction1(FTSize, beamDiameter, gaussianC, rLimit, x, baseOAM, ZPD, OAMPD, nPDs, ExpIntensity, colRangeIn, rowRangeIn, hAbsFT, hPhaseIFT, fitDxDy, verbosity), x, optimset ('Display', 'iter', 'TolFun', minParms{2}));  
  elseif (strcmp(minType, 'nonlin'))
	%[x_m, objf, cvg, outp] = nonlin_min( @(x) FTBeamPhaseDiversityErrorFunction1(FTSize, beamDiameter, gaussianC, rLimit, x, baseOAM, ZPD, OAMPD, nPDs, ExpIntensity, colRangeIn, rowRangeIn, hAbsFT, hPhaseIFT, fitDxDy, verbosity), x, optimset ('iter', true, 'TolFun', minParms{2}));
    [x_m, v_m, nevals] = fminunc(@(x) ... 
                                  FTBeamPhaseDiversityErrorFunction1(FTSize, beamDiameter, gaussianC, rLimit, x, baseOAM, ZPD, OAMPD, nPDs, ExpIntensity, colRangeIn, rowRangeIn, hAbsFT, hPhaseIFT, fitDxDy, verbosity), x, ...
                                  optimoptions ('fminunc', 'Display', 'iter', 'OptimalityTolerance', minParms{2}, 'MaxIterations', minParms{1}, 'MaxFunctionEvaluations', minParms{3}, ...
                                  'StepTolerance', minParms{2}, 'HessUpdate', 'bfgs'));  
 
 end

  
  %% Final results
  if fitDxDy
    baseZAberrations = x_m(1:(end - 2));
    dx = x_m(end - 1);
    dy = x_m(end);
  else
    baseZAberrations = x_m;
    dx = 0; dy = 0;
  end

  zernikeAberrations=baseZAberrations;
  fprintf('\ndx = %g, dy = %g, zernikeAberrations:\n', dx, dy);
  disp(zernikeAberrations) % Print final result
  
  %%Final calculus
  
 [amplitudeCenter, IntensityCenter, ~, waveAberation, gridSize, colRangeSLM, rowRangeSLM] = CalcOAMBeamPhaseDiversityFTFromAberrations1(FTSize, beamDiameter, gaussianC, rLimit, zernikeAberrations, baseOAM, dx, dy, ZPD, OAMPD, nPDs);

 
 % warning ("on", "Octave:broadcast"); % To get back to normal

  %% Draw figures of final minimum
  
  for k = 1:nPDs
    set((hAbsFT(k)), 'cdata', abs(IntensityCenter(colRangeIn,rowRangeIn, k)));
  end
  set(hPhaseIFT, 'cdata', angle(waveAberation(:, :,1)));
  drawnow;
  %figure(currentImage++),plot(nevals, v_m);
  %drawnow;
  
  %% End


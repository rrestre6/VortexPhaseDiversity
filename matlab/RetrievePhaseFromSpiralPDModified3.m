%%
% File: RetrievePhaseFromSpiralPDModified3.m
% Author: Carlos Cuartas, Santiago Echeverri and Rene Restrepo
% Date: 19/07/2016
%%
%  This program retrieves phase by using Gerchberg-Saxton like phase diversity algorithm
%  procedure proposed by Sharma et. al. See Ref[1].
%%
%  IN
%	BeamImage: Experimental image.
%	FTSize: Size after creating a padarray with before fft's are computed.
%	gridSize: Pixels used to create grids.
%	gaussianC: Gaussian waist.
%	rLimit: Maximum pupil radius in Fourier field.
%	order: Vortex topological charge.
%	minParms: Minimization parameters: max iterations
%	verbosity: Verbose level.
%
%%
%  OUT
%	Phase: Retrieved phase.
%	Intensity: Final reconstructed intensity.
%	maskp: Phase pupil empĺoyed
%
%%
% octave-3.8.1.
%
% [1] M. K. Sharma, C. Gaur, P. Senthilkumaran, and K. Khare. " Phase imaging
% using spiral-phase diversity". Appl. Opt. 54(13) 3979-3985, 2015.
%
%%

%% Function definition

function [Phase, intensity, ZBaseAberration] = RetrievePhaseFromSpiralPDModified3(beamImage, FTSize, gridSize, gaussianC, rLimit, lDiversities, kDiversities, minParms, verbosity)

	warning('off'); %Get rid of multiple messages

	%% Experimental images normalization
	nPD = size(beamImage,3);

	for k = 1:nPD
		beamImage(:,:,k) = beamImage(:,:,k) - min(min(beamImage(:,:,k)));
		beamImage(:,:,k) = beamImage(:,:,k) / max(max(beamImage(:,:,k)));
	end

	%% Constants
	dx = 0;
	dy = 0;
	imageSize = size(beamImage, 1); %Image size
	padSize=([FTSize FTSize 0] - [imageSize imageSize 0]) / 2; %Padarray size

	% Add zeros pad to experimental images, sqrt is taken as the field amplitude
	for k = 1:nPD
		beamImagePadded(:,:,k) = padarray(sqrt(beamImage(:,:,k)), padSize, 0, 'both');
	end

	ExpAmp = beamImagePadded; %Experimental images padded
	clear beamImagePadded;

	% Create grids
	[X, Y] =  meshgrid(-1:(2 / gridSize):(1 - 2 / gridSize), -1:(2 / gridSize):(1 - 2 / gridSize));
	X2 = X - dx;
	Y2 = Y - dy;
	[Theta, R] = cart2pol(X, Y);
	[Theta2, R2] = cart2pol(X2, Y2);

	%Create masks
	mask = zeros(gridSize); % Pupil mask
	mask(R < 1) = 1;
	maskp = zeros(gridSize); % Phase mask
	maskp(R < rLimit) = 1;
	mask2=zeros(gridSize); %Create a new pupil with in order to use Zernike polynomials
	mask2(R>rLimit)=NaN;

	%Scaling of R in order to have a normalized aperture.
    rPupil = R / rLimit;
	% 1D matrix with coordinates within the circle to use Zernike.m
    RVect=rPupil(~isnan(mask2));
    TVect=Theta(~isnan(mask2));

	%%Image Range Calculus
	colRangeSLM = (FTSize - gridSize) / 2 + 1:(FTSize + gridSize) / 2;
	rowRangeSLM = (FTSize - gridSize) / 2 + 1:(FTSize + gridSize) / 2;

	%% Entrace field

	%The aberrations are with respect to pupil (X, Y) and only the vortex with respect to X2, Y2.
	gaussAmp = exp(-((((X2).^2)/((gaussianC)^2))+(((Y2).^2)/((1*gaussianC)^2))));
	%Gaussian distribution is the initial simulated amplitude

	Uin = mask .* gaussAmp .* exp(-1i.*(zeros(gridSize))); %Entrace field
	UinPadded = padarray(Uin, ([FTSize FTSize] - gridSize) / 2, 0, 'both');
	ampIni=UinPadded;


    % The area element is related to the spacing, so
    dA = (2 / gridSize / rLimit) ^ 2;
    % This would be originally dA = dx * dy, with dx = dy = 2 / gridSize as
    % defined in the meshgrid above, but because of the change in rPupil needs
    % to be scaled too

    n_z = 37; % Number of elements for the basis. Also number of coefficients

    [Z,~,~]=Zernike(0:n_z,RVect,TVect); % by Jos� Alonso.

    % Now the interior product of all Zernikes with themselves is
    zProds = Z.'* Z * dA;

    % Or just use this numerical  integration to normalize them:
    zNorm = bsxfun(@rdivide, Z, sqrt(diag(zProds).'));

	%% FFT of entrace field
	Uout= fftshift(fft2(fftshift(UinPadded)));

	% Population of basicSpirals cell

	basicSpirals = cell(5,1);
	spiral = 1;
	for order = [-2, -1, 0, 1, 2]
        basicSpirals(spiral) = {spiralGenTo(ones(size(X2)),order)};
	    spiral = spiral+1;
	end

	% Initial phase guess
	phaseIni = exp(1i*zeros(FTSize));

	%% Create phase diversities and antidiversities
	[Diversities, AntiDiversities] = generate_diversities_for_GS(RVect, TVect, mask, mask2 , rLimit , basicSpirals , lDiversities, kDiversities);

	% Masks padded
	mask2 = padarray(mask2, padSize, NaN, 'both');
	mask = padarray(mask, padSize, 0, 'both');
	maskp = padarray(maskp, padSize, 0, 'both');

	% Diversities padded
	for k = 1:nPD
		DiversitiesPad{k} = padarray(Diversities{k}, padSize, 0, 'both').*mask;
		AntiDiversitiesPad{k} = padarray(AntiDiversities{k}, padSize, 0, 'both').*mask;
	end
	clear Diversities AntiDiversities

	%% Compute GS procedure

	MaxIterations = minParms;
	iteration = 1;

	ObjectPlane = ExpAmp(:,:,1) .* phaseIni;


	% Start loop
	while iteration <= MaxIterations % Stop if maximum iteration are reach

	ImagePlane = ifftshift( ifft2( fftshift( ObjectPlane ) ) ).* mask;

		for m = 1:nPD

			ImagePlane2 = ImagePlane .* DiversitiesPad{m} .* mask;
			ObjectPlane = fftshift( fft2( ifftshift ( ImagePlane2 ) ) );

			%% --> step 4
			% keep phase and apply A_2
			phaseIni = angle( ObjectPlane );

			phaseIni(FTSize/2+1,FTSize/2+1) = (sum(sum(phaseIni(-1+(FTSize/2+1):1+(FTSize/2+1), -1+(FTSize/2+1):1+(FTSize/2+1)),2),1) - phaseIni(FTSize/2+1,FTSize/2+1) ) / 8;

			phaseIni = exp(1i * phaseIni);

		%PSF retrieved
%		figure(1), imagesc(abs(ObjectPlane(colRangeSLM,rowRangeSLM)).^2), title('PSF recuperada');

%figure(3),imagesc(ExpAmp(colRangeSLM,rowRangeSLM,m));

		%Phase retrieved
%		figure(2),imagesc(angle(ImagePlane(colRangeSLM, rowRangeSLM)).* maskp(colRangeSLM, rowRangeSLM)), colorbar, title('Aberracion reconstruida'); colorbar; axis equal;
%drawnow

			ObjectPlane2 = ExpAmp(:,:,m) .* phaseIni;

			ImagePlane = ifftshift( ifft2( fftshift( ObjectPlane2 ) ) ).* mask .* AntiDiversitiesPad{m};

			phaseIni = angle( ImagePlane );

			phaseIni(FTSize/2+1,FTSize/2+1) = (sum(sum(phaseIni(-1+(FTSize/2+1):1+(FTSize/2+1), -1+(FTSize/2+1):1+(FTSize/2+1)),2),1) - phaseIni(FTSize/2+1,FTSize/2+1) ) / 8;

			phaseIni = exp(1i * phaseIni);

			ImagePlane = abs(ImagePlane) .* phaseIni;

%______________________________________________________________________________________________________________
			% Convert the aberrations from GS to a vector representation.
      %  AberrationGS_OAM_temp = AberrationGS_OAM(colRangeSLM,rowRangeSLM);
      %  clear('AberrationGS_OAM');
        %AberrationGS_OAM = exp(1i.*AberrationGS_OAM);
        phaseIni(isnan(mask2)) = 0;

        %eval(['Aberration',indian_variable_names{lDiv},'{',num2str(graphCounter),'}','=AberrationGS_OAM;']);
        phaseIni_vec = phaseIni(~isnan(mask2));

       % Projection of the aberrations into the Zernike basis gives the
       % Zernike coefficients that compose this wavefront
        ProjectionGS_OAM = ((angle(phaseIni_vec)) .' * zNorm) * dA;
        ZBaseAberration = ProjectionGS_OAM(4:18);
        order = 0;
        [~, ~, ~, phaseIni2, ~, ~, ~] = CalcOAMBeamFTFromAberrations2(FTSize, gridSize, gaussianC, rLimit, ZBaseAberration, dx, dy, order);
%______________________________________________________________________________________________________________

			%figure(3),  imagesc(angle(AntiDiversitiesPad{m}(colRangeSLM,rowRangeSLM)));
			%figure(4),  imagesc(angle(DiversitiesPad{m}(colRangeSLM,rowRangeSLM)));
		end
		phaseFin = angle(ImagePlane);
%         phaseIni=padarray(phaseIni2, padSize, 0, 'both');
%         ImagePlane = abs(ImagePlane) .* phaseIni;
		ObjectPlane = fftshift( fft2( ifftshift ( ImagePlane ) ) );
		phaseIni = angle(ObjectPlane);
		%phaseIni = flipud(fliplr(phaseIni));
		if verbosity ==1
			%figure(1), imagesc(flipud(fliplr(angle(ImagePlane(colRangeSLM,rowRangeSLM)).* maskp(colRangeSLM,rowRangeSLM))),[-pi pi]); title('Aberracion reconstruida');
            figure(1), imagesc(flipud(fliplr(angle(phaseIni2).* maskp(colRangeSLM,rowRangeSLM))),[-pi pi]); title('Aberracion reconstruida');
			figure(2), imagesc( abs(ObjectPlane(colRangeSLM,rowRangeSLM)).^2); title('PSF recuperada');
		end

		ObjectPlane = ExpAmp(:,:,1) .* exp(1i * phaseIni);

		iteration = iteration + 1;
	end
	%% Outputs
	Phase = flipud(fliplr(phaseFin(colRangeSLM, rowRangeSLM))).* maskp(colRangeSLM,rowRangeSLM);
    intensity = abs(ObjectPlane(colRangeSLM,rowRangeSLM)).^2;
end

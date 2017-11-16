%% Calculating the tensor from diffusion weighted images
%
% We calculate the tensor (Q) from diffusion weighted images.  We call the
% tensor Q because the tensor is a positive-definite Quadratic form.
%
% These are the key actors and relations.
%
% The 3x3 positive-definite quadratic form, Q, predicts the apparent
% diffusion data using the formula
%   
%   ADC = u' Q u
%
% where u is a unit vector that defines the direction.
%
% The diffusion signal is related to the ADC as
%
%   dSig = S0 exp(-b ADC)
%
% The diffusion ellipsoid (mean diffusion distance) is related to Q by
% finding the vectors that solve
%
%     1 = v' inv(Q) v
%
% (c) Stanford VISTA Team

%% Load diffusion weighted imaging data
clear all;
% The vistadata diffusion sample data are 40-directions.  The directory
% contains the dwi data as well as the bvals and bvecs. 
% dataDir = fullfile(mrvDataRootPath,'diffusion','sampleData');
% dwi     = dwiLoad(fullfile(dataDir,'raw','dwi.nii.gz'));
dataDir = {};
% HCP
% dataDir{1} = '/Users/gari/Documents/STANFORD_PROJECTS/dr/ANALYSIS/DWI/STRT20025/1000';
% dataDir{2} = '/Users/gari/Documents/STANFORD_PROJECTS/dr/ANALYSIS/DWI/STRT20025/2000';
% dataDir{3} = '/Users/gari/Documents/STANFORD_PROJECTS/dr/ANALYSIS/DWI/STRT20025/3000';
% dwi{ii}     = dwiLoad(fullfile(dataDir{ii},'dwi_aligned_trilin_noMEC.nii.gz'));

% BCBL
cd('/Users/gari/Documents/Documents - GariMBP/')
dataDir{1} = 'STANFORD_PROJECTS/dr/ANALYSIS/DWI/S005/1000';
dataDir{2} = 'STANFORD_PROJECTS/dr/ANALYSIS/DWI/S005/2500';
dataDir{3} = 'STANFORD_PROJECTS/dr/ANALYSIS/DWI/S005/MS';

dwi      = {};
fatrix   = {};
dttrix   = {};
for ii = [1,2,3]
    dwi{ii}     = dwiLoad(fullfile(dataDir{ii},'data_aligned_trilin_noMEC.nii.gz'));
    fatrix{ii}  = read_mrtrix(fullfile(dataDir{ii},'data_aligned_trilin_noMEC_fa.mif'));
    dttrix{ii}  = read_mrtrix(fullfile(dataDir{ii},'data_aligned_trilin_noMEC_dt.mif'));
end

% The dwiPlot routine enables visualization of simple objects
% dwiPlot(dwi,'bvecs');

%% Calculate the ADC from the diffusion-weighted data

% All the bvecs and bvals
% bvecs = dwi.bvecs;
% bvals = dwi.bvals;

% Find the non-difussion bvecs and bvals (b ~= 0).
for ii = [1,2,3]
    bvecs{ii} = dwiGet(dwi{ii},'diffusion bvecs');
    bvals{ii} = dwiGet(dwi{ii},'diffusion bvals');
end

% Pick a coordinate
% The points showing the ADC
% cCoords = [47 54 43];  % Circular
% dCoords = [44 54 43];  % Directional

% coords = [44 54 43];  % Directional

% MINI
% coords = [42 41 36];  % Selected one for the S005 MINI on the corpus callosum
% HCP
coords = [34, 41, 40];  % Selected one for the STRT20025 HCP on the corpus callosum


% It is possible to pick multiple, but not for this script
% coords = [47 54 43; 44 54 43]; % Both
% S0 = dwiGet(dwi,'b0 image',coords);
% GLU: I have several b0, 6 with one phase encoding, 5 with opposite phase
% encoding, I will create the average S0
% S0 = mean(S0);
for ii = [1,2,3]
    S0{ii} = dwiGet(dwi{ii},'b0 image',coords);
    S0{ii} = mean(S0{ii});

    % Diffusion data are nCoords x nImages, excludes the b=0 measures
    dSig{ii} = dwiGet(dwi{ii},'diffusion data image', coords);

    % dSig = S0 * exp(-b * ADC)
    % ADC = bVec*Q*bVec'
    ADC{ii} = - diag( (bvals{ii}).^-1 )  *  log(dSig{ii}(:)/S0{ii});  % um2/ms
end


figure;
subplot(3,3,1); 
xvector = [1:length(dSig{1}),1:length(dSig{2}),1:length(dSig{3})]';
yvector = [dSig{1}(:); dSig{2}(:); dSig{3}(:)];
colors  = [repmat([1,0,0],length(dSig{1}),1); ... % RED
           repmat([0,1,0],length(dSig{2}),1); ... % GREEN
           repmat([0,0,1],length(dSig{3}),1)];   % BLUE
scatter(xvector, yvector, 20, colors, 'filled');
title('DWI values: b1000=red, b2000=green, b3000=blue')


subplot(3,3,3); 
xvector = [1:length(dSig{1}),1:length(dSig{1}),1:length(dSig{1})]';
yvector = [ADC{1}(:); ADC{2}(:); ADC{3}(:)];
colors  = [repmat([1,0,0],length(dSig{1}),1); ... % RED
           repmat([0,1,0],length(dSig{1}),1); ... % GREEN
           repmat([0,0,1],length(dSig{1}),1)];   % BLUE
scatter(xvector, yvector, 20, colors, 'filled');
title('ADC values: b1000=red, b2000=green, b3000=blue')



% mrvNewGraphWin
% plot(ADC)
% xlabel('Direction list')
% ylabel('ADC (diffusion weighted)')

% In practice, don't compute as above.  But use the call
% ADC = dwiGet(dwi,'adc data image',coords);

% To see the points in 3D
% dwiPlot(dwi,'adc',ADC);


%% Solve for the quadratic, Q, that predicts the ADC values. 
%
% The ADC values are measured.  Each measurement has a direction and
% amplitude, m = bvecs*bval.  We want to find Q such that 
%
%    (bvec(:)'* Q * bvec(:) = ADC
%
% We express this as a quadratic equation.  Suppose the entries of Q are
% (sorry about this) qij.  Suppose that for each direction, bvec, we have a
% particular level of b.  Then we have for each direction,
%
%  ADC(:) = q11*b(:,1)^2 + ... 2*qij*b(:,i)*b(:,j) + ... q33*b(:,3)^2
%
% The coefficients qij are the same in every equation.  So, we can pull
% them out into a column vector.  We have a matrix, V, whose entries are
% these b values.
%
%   ADC = V*q
%

% Turn this into a function
%function Q = dwiTensor(dwi.bvecs,dwi.bvals,ADC);

% We compute using only the diffusion-weighted data.
%
% Because the bvecsP enter the equation on both sides of Q, we take the
% square root of the bvals and multiply it times the bvecs.
% b = bvecs .* repmat(bvals.^0.5,1,3);  % Scaled, non-zero bvecs
% b = bvecs;  % Do not scale.  We scale when we compute diffusion. 



% Here is the big matrix
for ii=[1,2,3]
    b{ii} = bvecs{ii};
    V{ii} = [b{ii}(:,1).^2, b{ii}(:,2).^2, b{ii}(:,3).^2, ...
             2* b{ii}(:,1).*b{ii}(:,2), 2* b{ii}(:,1).*b{ii}(:,3), ...
             2*b{ii}(:,2).*b{ii}(:,3)];

    % Now, we divide the matrix V by the measured ADC values to obtain the qij
    % values in the parameter, tensor
    tensor{ii} = V{ii}\ADC{ii};

    % We convert the format from a vector to a 3x3 Quadratic
    Q{ii} = dt6VECtoMAT(tensor{ii});  % eigs(Q)
    % svd(Q)
end
% end of function here 
% 
% To compare the observed and predicted, do this
for jj=[1,2,3]
    ADCest{jj} = zeros(size(ADC{jj}));
    for ii=1:size(bvecs{jj},1)
        u = bvecs{jj}(ii,:);
        ADCest{jj}(ii) = u(:)'*Q{jj}*u(:);
    end
end

% We would like to have a measure of how well the model does compared to a
% repeated measure.  Or, we would like a bootstrap of the estimate.  That
% is done in dtiRawFitTensor
% mrvNewGraphWin;
% for jj=[1,2,3]
%     subplot(2,3,jj);plot(ADCest{jj},ADC{jj},'o');
%     tck = (0:0.5:3); axis equal; grid on; set(gca,'xtick',tck,'ytick',tck)
%     xlabel('ADC (estimated)');ylabel('ADC (observed)');title(['b' num2str(jj) '000'])
% end

%% Here is the ellipsoid associated with the tensor.

% We have a prediction ADC = u'Qu.  The diffusion ellipsoid associated with
% Q are the scaled unit vectors such that v'inv(Q)v = 1.  The vectors v =
% sqrt(ADC)*u will be close to the ellipsoid defined by Q.
subplot(1,3,1); ellipsoidFromTensor(Q{1},[],[],0);title('(b1000) Diffusion dist. ellipsoid')
subplot(1,3,2); ellipsoidFromTensor(Q{2},[],[],0);title('(b2000) Diffusion dist. ellipsoid')
subplot(1,3,3); ellipsoidFromTensor(Q{3},[],[],0);title('(b3000) Diffusion dist. ellipsoid')


%% Here are the measured ADC values shown as vectors.
% These should be that peanut like shape
figure;
subplot(1,3,1); dwiPlot(dwi{1},'adc',ADC{1},Q{1});title('b1000')
subplot(1,3,2); dwiPlot(dwi{2},'adc',ADC{2},Q{2});title('b2000')
subplot(1,3,3); dwiPlot(dwi{3},'adc',ADC{3},Q{3});title('b3000')

% pts = diag(ADC(~b0).^0.5)*bvecs(~b0,:);
% plot3(pts(:,1),pts(:,2),pts(:,3),'b.');  % pts or pts/b?

%% A comparison of the estimated and observed ADC values
% mrvNewGraphWin;
% subplot(1,2,1)
% err = 100*(ADCest - ADC) ./ ADC;
% hist(err)
% xlabel('Percent Error')
% ylabel('N voxels')
% 
% mean(err)
% std(err)
% 
% err3 = diag(err)*bvecs;
% subplot(1,2,2)
% plot3(err3(:,1),err3(:,2),err3(:,3),'o')
% grid on
% axis equal
% 
% %% A comparison of the estimated and observed signal values
% %
% %    d = S0 exp(-b (u'Qu))
% %
% dSigEst = S0*exp(-bvals.* diag(bvecs*Q*bvecs'));
% 
% mrvNewGraphWin
% subplot(1,2,1)
% plot(dSigEst,dSig,'o')
% xlabel('Estimated diffusion sig')
% ylabel('Measured diffusion sig')
% 
% subplot(1,2,2), 
% p = 100*(dSigEst(:) - dSig(:)) ./ dSig(:);
% hist(p,20);
% xlabel('% Error')
% 
% %% The error in 3-space
% p3 = diag(p)*bvecs;
% plot3(p3(:,1),p3(:,2),p3(:,3),'o')
% grid on
% axis equal
% 
% %% End

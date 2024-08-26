function dataStruct = videoCurrentGen(stack,time,xy,vB,fkB,tWin,tStep,varargin)

%dataStruct = videoCurrentGen(stack, time, xy, vBounds, fkBounds, Twin, Tstep {,plotFlag})
%
%  Returns current from video stacks at times in t
%  by converting a block of the stack to f,ky space to find
%  the velocity of the advected surface structure(usually foam).
%  INPUTS
%    stack - grayscale image of timestack, size [MxN]
%    time - time line (starting from zero) of stack, size [Mx1]
%    xy - x,y position [Nx2] of each pixel in a dimension (should be equally spaced!)
%    Twin - the time length of the FFT window (in points)
%    Tstep - time length to step the window (in points)
%		 vBounds - [minV maxV], units m/s or vector of desired velocity steps, set this
%      to empty, [], to use defaults [-3 3]
%    fkBounds = [fmin fmax kmin kmax], vector of frequency and wavenumber bounds
%      energy out of side of these bounds will be set to 0.  Useful to eliminate some
%      of the wave contamination that leaks in.  Set this to empty, [], to use defaults.
%    {plotFlag} - optional, if true (~=0) will display a running plot of the
%      data processing
%  OUTPUT fields in dataStruct returned:
%    meanV - OCM guess at the mean current for the timestep
%    t - time index for meanV
%    ci - the 95% conf. interval around meanV
%    cispan -  the width of ci
%    prob - the probability of the model fit
%    QCspan - the 95th percentile minus the 50th percentile of the timestack
%             histogram, used to measure the amount of video "texture"
%    stdV - the width (std. dev.) of the energy in velocity spectrum
%    vAngle - orientation of the pixel array (radians)
%


%  This code requires:
%     Optimization Tollbox (lsqcurvefit.m)
%     Statistics and Machine Learning Toolbox (nlparci.m)


% take care of inputs and constants
dt = abs(mean(diff(time*24*3600)));
y = cumsum([0; 0; sqrt(sum(diff(xy).^2,2))]);
%vAngle = angle(diff(xy([1 end],:)).*[1 i]);
dy = abs(mean(diff(y)));
N = size(stack,2);
M = tWin;
L = dy*N;
T = dt*M;
taper = bartlett(M)*bartlett(N)';
if isempty(vB) % set the velocity bounds if empty
    vB = [-3 3];
end
dv = 0.05;
switch length(vB)
    case 2
        v = min(vB):dv:max(vB); %m/s, changed from +/-3 m/s
    otherwise
        v = vB;
end
sigma = 0.075;
% plotFlag = 0;
if ~isempty(varargin)
    plotFlag = varargin{1};
end
UB = [inf max(vB) 2 inf]; % upper bounds on search
LB = [0 min(vB) .01 0]; % lower bounds on search
opts = optimset('disp','off');

% make k vector
if rem(N,2) == 0
    k = (((-N/2)):((N/2)-1))/L;
else
    k = ((-(N-1)/2):((N-1)/2))/L;
end
gk = find((k>0)&(k<1/(2*dy))); %find + k's
k = k(gk);

% make f vector
if rem(M,2) == 0
    f = (((M/2)):-1:(-(M/2-1)))'/T;
else
    f = (((M-1)/2):-1:(-(M-1)/2))'/T;
end
gf = find((abs(f)~=0)&(abs(f)<1/(2*dt))); %find +/- f's
f = f(gf)/sign(mean(diff(y)));

% cut out low wavenumbers/frequencies
if isempty(fkB) % set the f and k bounds defaults if needed
    fkB = [0 inf min(k)*2 2];
end
K = repmat(k,length(f),1);
F = repmat(f,1,length(k)); 
FKind = abs(F)>fkB(1) &  abs(F)<fkB(2) & abs(K)>fkB(3) & abs(K)<fkB(4);
Smask = nan*ones(size(K));
Smask(FKind) = 1;
fkny = [max(abs(F(~isnan(Smask(:))))) max(abs(K(~isnan(Smask(:)))))];

% step window and compute S(f,ky) and transf. to S(f,v)
j = 1;
Nb = floor((size(stack,1)-(tWin-tStep))/tStep); %blocks
warning off

% set up output structure
dataStruct.meanI = nan(1,Nb);
dataStruct.QCspan = nan(1,Nb);
dataStruct.meanV = nan(1,Nb);
dataStruct.stdV = nan(1,Nb);
dataStruct.prob = nan(1,Nb);
dataStruct.ci = nan(Nb,2);
dataStruct.cispan = nan(1,Nb);
dataStruct.SNR = nan(1,Nb);
dataStruct.t = nan(1,Nb);
dataStruct.beta = nan(4,Nb); 

for window = 0:(Nb-1)
    % window data and construct
    cind = ((window*tStep)+1):((window*tStep)+tWin);
    block = stack(cind,:);
    meanStack2 = mean(block(:));
    dataStruct.meanI(j) = meanStack2;
    dataStruct.stdI(j) = std(block(:));
    block = block-repmat(mean(block),tWin,1); % block = stack minus mean
    
    % try another QC factor
    p95 = prctile(block(:),[95 50]);
    dataStruct.QCspan(j) = p95(1) - p95(2); 
    
    % calculate f-k spectrum
    stxfft = fft2(block.*taper);
    S_fk = 2*stxfft.*conj(stxfft)/(tWin*length(k));
    S_fk = fftshift(S_fk);
    S_fk = S_fk(gf,gk).*Smask;
    
    % calculate v-k spectrum 
    for ii = 1:length(k)
        S_vk(:,ii) = interp1(f/k(ii),S_fk(:,ii),v','linear');
    end
    
    % calculate S(v)
    V0 = sum(S_vk, 2, 'omitnan')';
    V1 = colfilt(V0,[1 5],'sliding',@mean);
    S_v = V1/max(V1);  %normalized spectrum
    [maxv, maxvind] = max(S_v); % maxv is the amplitude of the Gausian foam 
    mdV = v(maxvind); % mdV = v^bar
    if j == 1
        gind = find(~isnan(S_v));
    end
    
    % fit the S(v) to the model [Vamp Vmu Vsigma offset]
    % in the paper, beta0 parameters are referred to as [A_foam vbar sigma_foam A_noise]  
    beta0 = [maxv mdV 0.25 mean(S_v, 'omitnan')];
    jv = v(gind);
    vFun = @(beta0,jv) SofVwaveI(beta0,jv,fkny);
    try % try fit, if it fails you get no plots and get nans for output
        [beta,resnorm,resid,~,~,~,jacb] = lsqcurvefit(vFun,double(beta0),jv,double(S_v(gind)),LB,UB,opts);
        fitted = SofVwaveI(beta,v(gind),fkny);
        dataStruct.meanV(j) = beta(2);
        dataStruct.stdV(j) = beta(3);
        
        % find chi^2 of fit and statistical significance
        chi2(j) = resnorm/(sigma^2);
        dataStruct.prob(j) = 1 - chi2cdf(chi2(j),length(gind)-length(beta));
        
        % find 95% ci on the beta parameters
        cint = nlparci(beta,resid,jacb);
        dataStruct.ci(j,:) = cint(2,:);
        dataStruct.cispan(j) = abs(diff(cint(2,:)));
        dataStruct.SNR(j) = beta(1)/beta(4); % model guess at signal-to-noise
        dataStruct.beta(:,j) = beta; 
        dataStruct.t(j) = mean(time(cind));

        % Note to revisit
         if plotFlag
            figure
            subplot(221)
            imagesc(y,1:size(block,1)*dt,block)
            hold on
            ylabel('time (s)','fontsi',14)
            xlabel('y position (m)','fontsi',14)
            %plot([min(y) ],[])
            subplot(222)
            imagesc(f,k,log10(abs(S')))
            shading flat
            xlabel('frequency (Hz)','fontsi',14),ylabel('wavenumber (1/m)','fontsi',14)
            axis xy
            axis([-abs(fkny(1)) abs(fkny(1)) 0 fkny(2)])
            grid on
            colorbar
            subplot(223)
            imagesc(v,k,log10(abs(Sv')))
            shading flat
            colorbar
            xlabel('velocity (m/s)','fontsi',14),ylabel('wavenumber (1/m)','fontsi',14)
            axis xy
            axis([min(jv) max(jv) 0 fkny(2)])
            grid on
            subplot(224)
            plot(v,V,'linew',1) %semilogy(v,V,'linew',1)
            xlabel('velocity (m/s)','fontsi',14),ylabel('wavenumber (1/m)','fontsi',14)
            hold on
            plot(v(gind),fitted,'r--','linew',2)
            plot([0 0]+dataStruct.meanV(j),[0 1],'k','linew',2)
            plot([-1 1;-1 1]*dataStruct.stdV(j)+dataStruct.meanV(j),[0 0;1 1],'--k','linew',1)
            title(dataStruct.meanV(j))
            grid on
            legend('S(v)','S_{model}(v)')
            if plotFlag == 2
                pause
            else
                drawnow
            end
         end

        j = j+1;
        wherestep = [num2str(window + 1) ' of ' num2str(Nb)];
        fprintf(1,'	step %s		\r',wherestep);
    catch
        warning('Nonlinear fit failed - skipping this record')
    end
end
warning on


%%%%%%%%%%%%%%%%%%%%%%%%%
function S = SofVwaveI(bin,vin,fknyq)

% S = SofVwaveI(beta,V,fknyq)
%
% Produces the noise floor, S(v), given beta
%	([Vamp,Vmu,Vsigma,offset]) and the velocity vector, V.
%	fknyq = [knyq fnyq]

% set constants
fnyq = fknyq(1);
knyq = fknyq(2);
fL = 1/32;

% define betas
amp = bin(1); % >= 0
mu = bin(2); %  -3 <= mu <= 3
sigma = bin(3); % >= 0
offset = bin(4); % >= 0

% account for integration area in background fit
Vint = nan*ones(size(vin));
vLowInd = abs(vin) <= fnyq/knyq;
vHiInd = find(abs(vin) > fnyq/knyq);
Vint(vLowInd) = (knyq^2)/2;
Vint(vHiInd) = (fnyq^2)./(2*vin(vHiInd).^2);
Vint = Vint/max(Vint);

S = amp*exp(-((vin-mu)/sigma).^2) + offset.*Vint;
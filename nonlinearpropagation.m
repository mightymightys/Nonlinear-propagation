%*********************************************************
%   NONLINEAR PROPAGATION
%
%   Developement started: 10/2017
%
%   Author: Stefan Haessler haessler@ensta.fr;
%
%*********************************************************
%
%   Description: 1D NONLINEAR PROPAGATION of a short laser pulse
%	in a medium with dispersion and 3rd order nonlinearity	
%
%   Changelog: --
%
%*********************************************************
clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% COMPRESSION
      k2comp = -27; % GVD of compressor in fs^2
      k3comp = 25; % TOD of compressor in fs^3

    automaticcompression = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
close all

%switch on or off different effects
SPM=true;
SS=true;
DISP=true;
ATTN=true;


disp('****************************************************************** ')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT PULSE   
    E_0 = 6.5 ; %input pulse energy (mJ)
    tau_0 = 25; %input pulse duration FWHM in fs 
    lambda_0 = 800; %input carrier wavelength in nm
    omega_0 = 2*pi*299792458 / lambda_0 *1e-6;  % in PHz

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    k20 = 0; % input pulse GDD fs^2
    k30 = 0; % input pulse TOD fs^3

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    polariz=1; %0...linear; 1...circular (only reduce n2 by 2/3)
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calc parameters
    zMax = 2500;  % total propagation distance in mm   

    dz =  zMax/2500; % integration step in mm
    z = 0:dz:zMax; % propagation axis
    nPoints = 2^13;
    [nu,t] = ftAxis(nPoints,4); %M. Joffre's ftAxis(nPoints, nuMax) % nu up to 4 PHz
                                %(to sample well the field oscillations, just so you can make pretty plots)
                                %(and to resolve steep pulsefronts due to SS)
    omega = (2*pi)*nu;
    lambda = 299792458e-15 ./nu *1e9;  % WL in nm
    dt = t(2)-t(1);

    spcgrmwinwidth = tau_0*2/3;  % width of spectrogram temp. window (foot-to-foot) in units of t, i.e. here fs
                          % here you make the tradeoff btw. spec. and temp. resolution
    %spcgrmwinwidth = 60; % good for narrower spectra

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% PROPAGATION MEDIUM

    gas = mymenu('Which gas fills the fiber ?',{'Ar','Ne','He','as specified in script'},4);
    %gas = 3;
    gases = {'Ar','Ne','He','other'};
  
    fiberdiam = 0.536; %fiber inner diameter in mm

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %field attenuation constant in mm^-1 for lowest order fiber mode EH11
    if ATTN;  attn_alpha = (2.405/2/pi)^2 * (lambda_0*1e-6)^2 /(0.5*fiberdiam)^3 * .5*(1.5^2 +1) /sqrt(1.5^2-1);
    else      attn_alpha = 0;    
    end
    disp(['The fiber transmission for the pure EH11 mode is ',num2str(100*exp(-2*attn_alpha*zMax),'%10.2f'),' %.']);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Pressure profile / medium density
     p0 = 0.0; %gas pressure / density in bars at fiber entrance
     pL = 1.3; %gas input pressure / density in bars
     pZ = sqrt(p0^2 + z/zMax*(p0^2+pL^2));
     
%      p0 = 1.3 ;
%      pZ = p0*ones(size(z));
    
    if pZ(end) == pZ(1)
         infostring=['constant pressure of ',num2str(p0),' bars.'];
    else 
         infostring=['pressure gradient with ',num2str(pL),' bars at the end.'];
    end
    disp(['Filling gas is ',gases{gas},' with ',infostring]);
    
    if gas == 1 %Ar
        k2 = 0.018;   % Ar GVD fs^2/ mm / 1bar [APPLIED OPTICS 47, 4856 (2008)]
                %it's linear in pressure as long as density depends
                %linearly on pressure (ideal gas regime), i.e. at low pressures of up to 10 bars or so
        k3 = 0.012;   % Ar GVD fs^3/ mm / 1bar [APPLIED OPTICS 47, 4856 (2008)]
        %n20 = 19.4e-19; %cm^2 /W /bar, Ar [Opt. Express 18, 25847 (2010)]
        n20 = 1.27e-19; %cm^2 /W /bar, Ar [Phys. Rev. Lett., 106, 183902,(2011)]
        %n20 = 1e-19; %cm^2 /W /bar, Ar [Phys. Rev. Lett. 104, 103903 (2010)]
        Ibs = 2.5e14; %barrier suppression intensity in W/cm^2
    elseif gas == 2 %Ne
        k2= 0.0021;   % Ne GVD fs^2/ mm / 1bar [APPLIED OPTICS 47, 4856 (2008)]
        k3 = -0.006;   % Ne GVD fs^2/ mm / 1bar [APPLIED OPTICS 47, 4856 (2008)]
        %n20 = 1.8e-19; %cm^2 /W /bar, Ne [Opt. Express 18, 25847 (2010)]
        n20 = 0.13e-19; %cm^2 /W /bar, Ne [Phys. Rev. Lett., 106, 183902, (2011)]
        %n20 = 9e-21; %cm^2 /W /bar, Ne [APB 111, 447 (2013)]
        Ibs = 8.6e14; %barrier suppression intensity in W/cm^2
    elseif gas == 3 %He
        k2 = 0.0016 ;   % He GVD fs^2/ mm / 1bar pressure (or density unit) [refractiveindex.info],
        %k2 = 0.0008;   % He GVD fs^2/ mm / 1bar [APPLIED OPTICS 47, 4856 (2008)]
        k3 = -0.001;   % He GVD fs^3/ mm / 1bar [APPLIED OPTICS 47, 4856 (2008)]
        n20 = 3.6e-21; %cm^2 /W /bar, He 
        %n20 = 5e-21; %cm^2 /W /bar, He [Phys. Rev. Lett., 106, 183902, (2011)]
        Ibs = 14.6e14; %barrier suppression intensity in W/cm^2
    else  %other
        k2 = 36.3; % in FS, fs^2/ mmm
        k3 = 27.7 % in FS, fs^3/ mm
   
        %k2 = 45; % in BK7, fs^2/ mm
        %k3 = 32.3 % in BK7, fs^3/ mm
   
        %k2 = -0.007; %works for "self-compression" with really short spike on top of fat red pedestal
                % with 1.5bars He in 1m long, 250micron fiber, 1mJ, 24fs Fourier-limited input
        %k3 = 0.5*k2;
        n20 = 2.5e-16;   %cm^2 /W , FS @800nm 
        Ibs = 1e14;
    end
    
    if polariz==1
        n20=2/3*n20;
    end
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% set up initial electric field
carrier = exp(-1i*omega_0*t);
phase0 = k20/2*(omega).^2 + k30/6*(omega).^3; % Phase spectrale (dispersion)

    field_t = exp(-2*log(2)*t.^2/tau_0^2);
    field = ift(field_t) .*exp(1i*phase0);
    field_t = ft(field);
        tempintensity_in = field_t.*conj(field_t);
        [maxI,maxind] = max(tempintensity_in);
        highind = find(tempintensity_in>0.5*maxI,1,'last');
        lowind  = find(tempintensity_in>0.5*maxI,1,'first');

    tau_0 = interp1(tempintensity_in(highind:highind+1), t(highind:highind+1), 0.5*maxI)...
            - interp1(tempintensity_in(lowind-1:lowind), t(lowind-1:lowind), 0.5*maxI);
    
    % peak intensity at fiber input
    % find proportionality constant alpha such that  E_0 = alpha * \int abs(field_t).^2 dt; 
    % then power = dE/dt = alpha * abs(field_t).^2; and intensity = power/effective mode area (given by Vozzi as 0.48 times HCF core area
    I_0 = max(abs(field_t).^2) * E_0*1e-3/( sum(abs(field_t).^2)*(t(2)-t(1))*1e-15) /(0.48*pi*(.5*fiberdiam*0.1)^2);
    disp(['--']);
    disp(['The peak intensity at the fiber input (before coupling losses) is ',num2str(I_0,'%10.2e\n'),' W/cm^2,']);
    disp(['i.e. ',num2str(I_0/Ibs*100,'%10.0f\n'),'% if the barrier supression intensity of ',gases{gas},'.']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Autofocalization ?
    P_cr = (lambda_0*1e-7)^2 /2/pi/(n20* max(pZ))/1e9 ; %critical power for self-focusing in GW
    P_in = E_0*1e-3 /(tau_0*1e-15) /1e9; % power at input. Actually the peak power of a Gaussian is 0.94*that, but we're not that picky I suppose.
    disp(['--']);
    disp(['The peak power at the fiber input is ',num2str(P_in,'%10.0f\n'),' GW,']);
    disp(['i.e. ',num2str(P_in/P_cr,'%10.2f\n'),' times the cricital power for self-focusing for the gas pressure at the fiber exit of ',num2str(P_cr,'%10.0f\n'),' GW.']);
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Vozzi's ionization suppression ? [Vozzi et al., Appl. Phys. B 80, 285 (2005)]
    Avozzi = 2.62e-9 * [1.79, 1.14, 1, 2.08]; %Vozzi's prop. constant for estimating the min. fiber radius,
                                              %2.08 is for Kr
    amin = Avozzi(gas) * (tau_0*1e-15)^-.45 * (E_0*1e-3)^.51;
    disp(['--']);
    disp(['For this input pulse and gas, Catherina Vozzi advises a minimal fiber diameter of ',num2str(2*amin*1e6,'%10.0f\n'),' microns.']);
    disp(['Your fiber has a core diameter of ',num2str(fiberdiam*1e3,'%10.0f\n'),' microns.']);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Solve propagation diff. eq. with split-step method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Calculating propagation...');
fprintf('done: %03d percent',0)

map = zeros(length(z),nPoints); % initialize E-field(z,nu)
map_t = zeros(length(z),nPoints); % initialize E-field(z,t)
Bint = 0;
accNLphase_t = zeros(1,length(t));

field_t = field_t ;
field = field ;

EnergyCHECK = sum((abs(field_t)).^2); % this is representative of the pulse energy coupling into the fiber

map_t(1,:) = field_t.*carrier; % save initial temp. E-field, incl. carrier
map(1,:) = ift(field_t.*carrier); % save initial spec. E-field, incl. carrier

for iZ = 2:length(z)
  
  if rem(iZ,round(length(z)/100)) == 0;  fprintf('\b\b\b\b\b\b\b\b\b\b\b\b %03d percent', round(iZ/length(z)*100)); end
  
  % Time domain
  gamma = 2*pi/(lambda_0*1e-6)*pZ(iZ)*n20; %"nonlinear coefficient"
  
  if SPM
    phase_t_SPM = gamma*I_0*(abs(field_t)).^2*dz; % SPM temporal phase 
  else  
    phase_t_SPM = zeros(size(t));
  end
  
  if SS
    phase_t_SS  = 1i*I_0*gamma/omega_0 .*field_t.*ft(-1i*omega.*ift(conj(field_t)))*dz ...
                +2*1i*I_0*gamma/omega_0 .*conj(field_t) .*ft(-1i*omega.*ift(field_t))*dz; 
                % SS temporal phase, as in [Deiterding et al., Journal of Lightwave Technology 31, 2008 (2013)], Eq. (14)
  else          
      phase_t_SS = zeros(size(t));
  end
  
  Bint = Bint + max(phase_t_SPM) + max(phase_t_SS);  %B-integral as max accumulated nonlinear phase
  accNLphase_t = accNLphase_t +phase_t_SPM + phase_t_SS;
  
  field_t = field_t .*exp(1i*phase_t_SPM) .*exp(1i*phase_t_SS) *exp(-attn_alpha*dz); 
  
  % Spectral domain
  if DISP
    field = ift(field_t);
  
    phaseDISP = .5*k2*pZ(iZ)*dz*(omega).^2 + 1/6*k3*pZ(iZ)*dz*(omega).^3; % spectral phase (dispersion)
      
    field = field .*exp(1i*phaseDISP) ;
  
    field_t = ft(field); 
  end
  
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

  if rem(iZ,round(length(z)/100)) == 0 || iZ==length(z)
    % plot temporal intensity envelope
      figure(1)
      subplot(2,1,1)
          plot(t,(abs(field_t)).^2,'LineWidth',1.5);
          xlabel('Time (fs)');
          ylabel('E(t)');
          ylim([0 1]);
          xlim(2*[-tau_0 tau_0]);
        title(['z = ',num2str(z(iZ)),' mm']);
  % plot spectral intensity envelope
      subplot(2,1,2)
          plot(omega+omega_0,abs(field).^2,'LineWidth',1.5);
          xlabel('Frequency (PHz)');
          ylabel('|E(nu)|^2');
          xlim([1.5 3.5]);
          ylim([0 1.1*max(abs(field).^2)]);
  end
  drawnow

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

    if abs( sum((abs(field_t)).^2)*exp(2*attn_alpha*z(iZ)) - EnergyCHECK) > 1e-4*EnergyCHECK
      field_t = map_t(iZ-1,:); %go back one propagation step
      field = map(iZ-1,:);     %go back one propagation step
      map_t(end,:) = field_t .*carrier; % Save temp. field, incl. carrier wave
      map(end,:) = ift(field_t .*carrier); % Save spec. field, incl. carrier wave
      fprintf('\n');
      disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
      disp(['Energy is no longer conserved, your pulse is broken, probably due to wavebreaking.']);
      disp(['We STOP the propagation here at z = ',num2str(z(iZ)),' mm.']);
      disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
      break
    end 

  map_t(iZ,:) = field_t .*carrier; % Save temp. field, incl. carrier wave
  map(iZ,:) = ift(field_t .*carrier); % Save spec. field, incl. carrier wave
end;
fprintf('\n');
disp(['The accumulated B-integral is ',num2str(Bint,'%10.2f'),' rad.'])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% calculate final spectral width according to Eq. 9 in Pinault and Potasek, JOSA B 2, 1318 (1985), and multiply by 2 so you have an actual width
posspecpoints = (nPoints/2+10:nPoints);

spec = (abs(map(1,posspecpoints))).^2;
specint = sum(spec); 
inputrmsbandwidth = 2*sqrt( sum(omega(posspecpoints).^2 .*spec)/specint - (sum(omega(posspecpoints).*spec)/specint)^2); %in PHz

specint = sum(spec(2:end).*-diff(lambda(posspecpoints)));
inputrmsbandwidth_lambda = 2*sqrt( sum( lambda(posspecpoints(2:end)).^2.*spec(2:end).*-diff(lambda(posspecpoints)))/specint ...
                                - (sum( lambda(posspecpoints(2:end))   .*spec(2:end).*-diff(lambda(posspecpoints)))/specint)^2); %in nm


spec = (abs(map(end,posspecpoints))).^2;
specint = sum(spec);
specCM =  sum(spec.*omega(posspecpoints))/sum(spec); % center of mass of (asymmetrically) broadened spectrum
finalrmsbandwidth = 2*sqrt( sum(omega(posspecpoints).^2 .*spec)/specint - (sum(omega(posspecpoints).*spec)/specint)^2); %in PHz

specint = sum(spec(2:end).*-diff(lambda(posspecpoints)));
finalrmsbandwidth_lambda = 2*sqrt( sum( lambda(posspecpoints(2:end)).^2.*spec(2:end).*-diff(lambda(posspecpoints)))/specint ...
                                - (sum( lambda(posspecpoints(2:end))   .*spec(2:end).*-diff(lambda(posspecpoints)))/specint)^2); %in nm

clear spectrum
disp(['The final rms bandwidth is ',num2str(finalrmsbandwidth,'%10.2f'),' PHz or ',num2str(finalrmsbandwidth_lambda,'%10.2f'),' nm, after broadening by a factor ',num2str(finalrmsbandwidth/inputrmsbandwidth,'%10.2f'),'.'])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compression de l'impulsion en appliquant une GDD + TOD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if automaticcompression 
   [~,bar1]=min(abs(omega- (specCM-0.5*finalrmsbandwidth)));
   [~,bar2]=min(abs(omega- (specCM+0.5*finalrmsbandwidth)));
   P = polyfit(omega(bar1:bar2)-omega_0, unwrap(angle(map(end,bar1:bar2))),3);
   k2comp = -P(2)*2;
   k3comp = -P(1)*6;
end
% else
%     P(2) = -k2comp/2;
%     P(1) = -k3comp/6;
% end
compressedfield = map(end,:) .*exp(1i*.5*k2comp*(omega-omega_0).^2 + 1i/6*k3comp*(omega-omega_0).^3);
compressedfield_t(:) = ft(compressedfield);
compressedint(:) = (abs(compressedfield_t)).^2;
if max(compressedint(:))>0
    compressedtau = max(diff(t(abs(diff(sign(compressedint - 0.5*(max(compressedint)))))==2)));
else
    compressedtau = NaN;
end

disp(['After compression with ',num2str(k2comp),' fs^2 GDD and ',num2str(k3comp),' fs^3 TOD the pulse duration is = ',num2str(compressedtau),' fs FWHM.']);

FTlimitint = abs(ft(abs(map(end,:)))).^2;
if max(FTlimitint)>0
    FTlimittau = max(diff(t(abs(diff(sign(FTlimitint - 0.5*(max(FTlimitint)))))==2)));
    disp(['The Fourier limited pulse duration is = ',num2str(FTlimittau),' fs FWHM.']);
else
    FTlimittau = NaN;
end


    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plotting part...

    
%% make series of spectrograms
winchoice =2; % =1 for Gaussian, =2 for Hanning

nsteps =2;

steps = round(linspace(1,length(z),nsteps));
[spcgrm,spcgrmtime] = myspectrogram(map_t(1,:),t,spcgrmwinwidth,winchoice);
[u,v] = size(spcgrm);
spctrgrmmap = zeros(length(steps),u,v);
spctrgrmmap(1,:,:) = spcgrm;

    [spcgrm,spcgrmtime] = myspectrogram(map_t(iZ,:),t,spcgrmwinwidth,winchoice);
    spctrgrmmap(2,:,:) = spcgrm;
    spctrgrmmap=spctrgrmmap(1:2,:,:);


%%  
set(0,'defaultaxesfontsize',9); %Increase fontsize on figures.

freqlims=[1.5 4];
timelims=[-50 50];


%%
figure(5)

subplot(4,4,[1 9])
     area(omega,(abs(map(1,:))).^2,'FaceColor',[0.8 0.2 0.2]);
     xlim(freqlims)
     xlabel('Frequency (PHz)');
     ylim([0 1.1*max((abs(map(1,:))).^2)])
     set(gca,'ytick',[])
     ylabel('Intensity (arb.u.)');
     view(-90,90)
     axis ij
     
subplot(4,4,[2 12])
    imagesc(spcgrmtime,omega,squeeze(spctrgrmmap(1,:,:))');
    set(gca,'YDir','normal')
    ylim(freqlims)
    xlim(timelims);
    title('input pulse');

subplot(4,4,[14 16])
     area(t,(abs(map_t(1,:))).^2,'FaceColor',[0.8 0.2 0.2]);
     xlim(timelims)
     xlabel('Time (fs)');
     ylim([0 1.1*max((abs(map_t(1,:))).^2)])
     set(gca,'ytick',[])
     ylabel('Intensity (arb.u.)');
     

%%
figure(6)
%freqlims=[1.5 4];
%timelims=[-50 50];

subplot(4,4,[1 9])
     area(omega,(abs(map(end,:))).^2,'FaceColor',[0.8 0.2 0.2]);
     xlim(freqlims)
     xlabel('Frequency (PHz)');
     ylim([0 1.1*max((abs(map(end,:))).^2)])
     set(gca,'ytick',[])
     ylabel('Intensity (arb.u.)');
     set(gca,'YDir','normal')
     view(-90,90)
     axis ij
     
subplot(4,4,[2 12])
    clims=[-13 -8];
    imagesc(spcgrmtime,omega,squeeze(spctrgrmmap(end,:,:))');
    set(gca,'YDir','normal')
    ylim(freqlims)
    xlim(timelims);
    title('output pulse');
    
subplot(4,4,[14 16])
     area(t,(abs(map_t(end,:))).^2,'FaceColor',[0.8 0.2 0.2]);
     xlim(timelims)
     xlabel('Time (fs)');
     ylim([0 max([1, 1.1*max((abs(map(end,:))).^2)])])
     set(gca,'ytick',[])
     ylabel('Intensity (arb.u.)');


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[cmpspcgrm,spcgrmtime] = myspectrogram(compressedfield_t,t,spcgrmwinwidth,winchoice);

figure(4)
    imagesc(spcgrmtime,omega,cmpspcgrm');
    set(gca,'YDir','normal')
    ylim(freqlims)
    ylabel('Frequency (PHz)');
    xlim([-60 60]);
    xlabel('Time (fs)');
    title(['After compression with ',num2str(k2comp),' fs^2 GDD and ',num2str(k3comp),' fs^3 TOD.']);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
figure;
plot(t, compressedint/(sum(compressedint)*dt*1e-15)*E_0*1e-3*1e-12,'LineWidth',1.5);
 xlabel('Time (fs)');
  ylabel('Power (TW)');
  xlim([-5*compressedtau 5*compressedtau]);
  %ylim([0 2]);


  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
figure;

minwl = round((lambda_0 - 1.5*finalrmsbandwidth_lambda(end))/20)*20;
maxwl = round((lambda_0 + 1.5*finalrmsbandwidth_lambda(end))/20)*20;
[~,bar1]=min(abs(lambda-minwl));
[~,bar2]=min(abs(lambda-maxwl));
cmpspecphase = unwrap(angle(compressedfield));

P = polyfit(omega(bar2:bar1),cmpspecphase(bar2:bar1),1);
cmpspecphase = cmpspecphase - P(1)*omega - P(2);


[AX, H1, H2] = plotyy(lambda, 2*pi*299792458*lambda.^(-2).*(abs(compressedfield)).^2, lambda, cmpspecphase);  % scale dE/dlambda = dE/domega * domega/dlambda 
 set(H1,'LineWidth',1.5)
 set(H2,'LineWidth',1.5)
 xlabel('Wavelength (nm)');
 xlim(AX(2), [minwl maxwl]) 
 xlim(AX(1), [minwl maxwl])  
 ylim(AX(2), [-5 5]) 
 set(AX(2),'YTick',[-5 0 5]) 
 %set(AX(2),'YTick',[-15 0 15]) 
 
 ylabel(AX(1),'Spectrum') 
 ylabel(AX(2),'Phase spectrale (rad)') 


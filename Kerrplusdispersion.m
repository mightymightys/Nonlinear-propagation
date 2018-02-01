% Nonlinear propagation of a short laser pulse in a medium with dispersion and 3rd order nonlinearity	
close all
clear all

SPM=1;
SS=1;
DISP=1;
ATTN=1;

progressplots=1;
saveimgs=0; 
spcgrmmovie=0;

if saveimgs==0
    outfolder='C:\Users\haessler\Documents\SIMULATION results\Nonlinear propagation';
    if exist(outfolder, 'dir')==7
        warning('I will rename the old outfolder.')
        movefile(outfolder,[outfolder,'_old'])
    end
    mkdir(outfolder)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT PULSE
    %E_0 = 10 ; %input pulse energy (mJ)
    E_0 = 1.5;
    tau_0 = 25; %input pulse duration FWHM in fs 
    lambda_0 = 800; %input carrier wavelength in nm
    omega_0 = 2*pi*299792458 / lambda_0 *1e-6;  % in PHz
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    k20 = 0; % input pulse GDD fs^2
    %k20 = 400; % to demonstrate spec. narrowing in SPM with neg. input chirp
    %k30 = 0.5*k20; % input pulse TOD fs^3
    k30 = 0;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    polariz=1; %0...linear; 1...circular (only reduce n2 by 2/3)
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calc parameters
    %zMax = 2500;  % total propagation distance in mm
    zMax = 2000;    

    dz = .5; % integration step in mm
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
    gases = {'Ar','Ne','He','other'};
    %couplingeff = 0.7;
    couplingeff = 0.8;
    %fiberdiam = 0.536; %fiber inner diameter in mm
    %fiberdiam = 0.250; 
    fiberdiam = 0.320; 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %field attenuation constant in mm^-1 for lowest order fiber mode EH11
    attn_alpha = (2.405/2/pi)^2 * (lambda_0*1e-6)^2 /(0.5*fiberdiam)^3 * .5*(1.5^2 +1) /sqrt(1.5^2-1); 
    if ATTN==0;  attn_alpha = 0; end
    disp(['The fiber transmission for the pure EH11 mode is ',num2str(100*exp(-2*attn_alpha*zMax),'%10.2f'),' %.']);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Pressure profile / medium density
     p0 = 0.0; %gas pressure / density in bars at fiber entrance
     pL = 0.85; %gas input pressure / density in bars
     pZ = sqrt(p0^2 + z/zMax*(p0^2+pL^2));
    
     %p0 = 1.0;
     %pZ = p0*ones(size(z));
    
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
    elseif gas == 2 %Ne
        k2= 0.0021;   % Ne GVD fs^2/ mm / 1bar [APPLIED OPTICS 47, 4856 (2008)]
        k3 = -0.006;   % Ne GVD fs^2/ mm / 1bar [APPLIED OPTICS 47, 4856 (2008)]
        %n20 = 1.8e-19; %cm^2 /W /bar, Ne [Opt. Express 18, 25847 (2010)]
        n20 = 0.13e-19; %cm^2 /W /bar, Ne [Phys. Rev. Lett., 106, 183902, (2011)]
        %n20 = 9e-21; %cm^2 /W /bar, Ne [APB 111, 447 (2013)]
    elseif gas == 3 %He
        k2 = 0.0016 ;   % He GVD fs^2/ mm / 1bar pressure (or density unit) [refractiveindex.info],
        %k2 = 0.0008;   % He GVD fs^2/ mm / 1bar [APPLIED OPTICS 47, 4856 (2008)]
        k3 = -0.001;   % He GVD fs^3/ mm / 1bar [APPLIED OPTICS 47, 4856 (2008)]
        n20 = 3.6e-21; %cm^2 /W /bar, He 
        %n20 = 5e-21; %cm^2 /W /bar, He [Phys. Rev. Lett., 106, 183902, (2011)]
    else  %other
        k2 = 36.3; % in FS, fs^2/ mmm
        k3 = 27.7 % in FS, fs^3/ mm
   
        %k2 = 45; % in BK7, fs^2/ mm
        %k3 = 32.3 % in BK7, fs^3/ mm
   
        %k2 = -0.007; %works for "self-compression" with really short spike on top of fat red pedestal
                % with 1.5bars He in 1m long, 250micron fiber, 1mJ, 24fs Fourier-limited input
        %k3 = 0.5*k2;
        n20 = 2.5e-16;   %cm^2 /W , FS @800nm 
    end
    
    if polariz==1
        n20=2/3*n20;
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% COMPRESSION
    Avozzi = 2.62e-9 * [1.79, 1.14, 1, 2.08]; %Vozzi's prop. constant for estimating the min. fiber radius,
                                              %2.08 is for Kr
    amin = Avozzi(gas) * (tau_0*1e-15)^-.45 * (E_0*1e-3)^.51;
    disp(['For this input pulse and gas, Catherina Vozzi advises'])
    disp(['a minimal fiber diameter of ',num2str(2*amin*1e6),' microns.']);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% COMPRESSION
%      k2comp = -15; % GVD of compressor in fs^2
%      k3comp = 10; % TOD of compressor in fs^3
%     
%     k2comp = -20; % good for sim of our setup 
%     k3comp = 13; % good for sim of our setup 
    
    k2comp = 0; 
    k3comp = 0;

    
    
    
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% set up initial electric field
%T_0 = lambda_0 *1e-9 /299792458 *1e15;

carrier = exp(-1i*omega_0*t);

field_t = exp(-2*log(2)*t.^2/tau_0^2);
  phase0 = .5*k20*(omega).^2 + 1/6*k30*(omega).^3; % Phase spectrale (dispersion)
  field = ift(field_t) .*exp(1i*phase0);
field_t = ft(field); 

I_0 = max(abs(field_t).^2) * E_0*1e-3/sum(abs(field_t).^2) /(0.48*pi*(.5*fiberdiam*0.1)^2)/ (t(2)-t(1))*1e15;
     % peak intensity at fiber input
disp(['The peak intensity at the fiber input is ',num2str(I_0,'%10.2e\n'),' W/cm^2.']);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Solve propagation diff. eq. with split-step method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

map = zeros(length(z),nPoints); % initialize E-field(z,nu)
map_t = zeros(length(z),nPoints); % initialize E-field(z,t)
Bint = 0;
accNLphase_t = zeros(1,length(t));

field_t = field_t * sqrt(couplingeff);
field = field * sqrt(couplingeff);

map_t(1,:) = field_t.*carrier; % save initial temp. E-field, incl. carrier
map(1,:) = ift(field_t.*carrier); % save initial spec. E-field, incl. carrier

for iZ = 2:length(z)
  % Time domain
  gamma = 2*pi/(lambda_0*1e-6)*pZ(iZ)*n20; %"nonlinear coefficient"
  
  if SPM==1
    phase_t_SPM = gamma*I_0*(abs(field_t)).^2*dz; % Phase temporelle SPM 
  else  
    phase_t_SPM = zeros(size(t));
  end
  
  if SS==1
    phase_t_SS  = 1i*I_0*gamma/omega_0 .*field_t.*ft(-1i*omega.*ift(conj(field_t)))*dz ...
                +2*1i*I_0*gamma/omega_0 .*conj(field_t) .*ft(-1i*omega.*ift(field_t))*dz; 
                % as in Deiterding et al., Eq. (14)
  else          
      phase_t_SS = zeros(size(t));
  end
  
  Bint = Bint + max(phase_t_SPM) + max(phase_t_SS);  %B-integral as max accumulated nonlinear phase
  accNLphase_t = accNLphase_t +phase_t_SPM + phase_t_SS;
  
  field_t = field_t .*exp(1i*phase_t_SPM) .*exp(1i*phase_t_SS); 
  
  % Spectral domain
  field = ift(field_t);
  if DISP==1
    phaseDISP = .5*k2*pZ(iZ)*dz*(omega).^2 + 1/6*k3*pZ(iZ)*dz*(omega).^3; % spectral phase (dispersion)
  else
    phaseDISP = zeros(size(t));
  end
      
  field = field .*exp(1i*phaseDISP) *exp(-attn_alpha*dz);
  
  field_t = ft(field); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
if progressplots==1
  if rem(iZ,20)==0 || iZ==length(z)
    % plot temporal intensity envelope
      figure(1)
      subplot(2,1,1)
          plot(t,(abs(field_t)).^2,'LineWidth',1.5);
          %plot(t,real(field_t));
          xlabel('Temps (fs)');
          %ylabel('|u(t)|^2 ');
          ylabel('E(t)');
          ylim([0 1]);
          xlim(2*[-tau_0 tau_0]);
        title(['z = ',num2str(z(iZ)),' mm']);
  % plot spectral intensity envelope
      subplot(2,1,2)
          plot(omega+omega_0,abs(field).^2,'LineWidth',1.5);
          xlabel('Fréquence (PHz)');
          ylabel('|E(nu)|^2');
          xlim([1.5 4]);
          ylim([0 1.1*max(abs(field).^2)]);
  end
  drawnow
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

  map_t(iZ,:) = field_t .*carrier; % Save temp. field, incl. carrier wave
  map(iZ,:) = ift(field_t .*carrier); % Save spec. field, incl. carrier wave
end;

disp(['The accumulated B-integral is ',num2str(Bint,'%10.2f'),' rad.'])

 
  
  
%% make series of spectrograms
winchoice =2; % =1 for Gaussian, =2 for Hanning

if spcgrmmovie==0
    nsteps =2;
else
    nsteps=50;
end

steps = round(linspace(1,length(z),nsteps));
[spcgrm,spcgrmtime] = myspectrogram(map_t(1,:),t,spcgrmwinwidth,winchoice);
[u,v] = size(spcgrm);
spctrgrmmap = zeros(length(steps),u,v);
spctrgrmmap(1,:,:) = spcgrm;
if spcgrmmovie==1
    n=2;
    for iZ = steps(2:end)
      [spcgrm,spcgrmtime] = myspectrogram(map_t(iZ,:),t,spcgrmwinwidth,winchoice);
      spctrgrmmap(n,:,:) = spcgrm;
      n=n+1;
    end
else
    [spcgrm,spcgrmtime] = myspectrogram(map_t(iZ,:),t,spcgrmwinwidth,winchoice);
    spctrgrmmap(2,:,:) = spcgrm;
    spctrgrmmap=spctrgrmmap(1:2,:,:);
end

%%  
set(0,'defaultaxesfontsize',9); %Increase fontsize on figures.

freqlims=[1.5 4];
timelims=[-50 50];
% freqlims=[2 2.8];
% timelims=[-75 75];
%timelims=[-100 100];

if spcgrmmovie==1
for iZ = 1:length(steps)
    figure(3)
    subplot(4,4,[1 9],'replace')
     area(omega,(abs(map(steps(iZ),:))).^2,'FaceColor',[0.8 0.2 0.2]);
     xlim(freqlims)
     xlabel('Fréquence (PHz)');
     %ylim([0 max((abs(map(steps(iZ),:))).^2)])
     ylim([0 max((abs(map(steps(1),:))).^2)])
     set(gca,'ytick',[])
     ylabel('Intensité (arb.u.)');
     view(-90,90)
     axis ij
     
     subplot(4,4,[2 12],'replace')
        %clims=[-14 -9];
        %imagesc(spcgrmtime,omega,log(squeeze(spctrgrmmap(iZ,:,:))'),clims);
        imagesc(spcgrmtime,omega,squeeze(spctrgrmmap(iZ,:,:))');
        set(gca,'YDir','normal')
        load('MyGermanColormap','mycmap')
    %    load('MyBrightJetColormap','mycmap')
        set(gcf,'Colormap',mycmap)
        ylim(freqlims)
        %ylabel('Fréquence (PHz)');
        xlim(timelims);
        %xlabel('Temps (fs)');
         title(['z = ',num2str(z(steps(iZ))),' mm']);
        
     subplot(4,4,[14 16],'replace')
         area(t,(abs(map_t(steps(iZ),:))).^2,'FaceColor',[0.8 0.2 0.2]);        
         xlim(timelims)
         xlabel('Temps (fs)');
         %ylim([0 1.1*max((abs(map_t(steps(iZ),:))).^2)])
         ylim([0 1.1*max((abs(map_t(steps(1),:))).^2)])
         set(gca,'ytick',[])
         ylabel('Intensité (arb.u.)');
     %pause(0.01)
    drawnow
    
    if saveimgs==1 
        paperwidth=18;
        paperheight=14;
           set(gcf,'Renderer','painters',...
                  'PaperPositionMode', 'manual','PaperUnits','centimeters',...
                  'PaperSize',[paperwidth paperheight],'PaperPosition',[0.25 0.25 paperwidth-0.5 paperheight-0.5]);    

        outname = [outfolder,'/_spectrogramsout',num2str(iZ,'%0.2d')];
        print('-r300',[outname,'.png'],'-dpng');
    end
end
end


%%
figure(5)
%freqlims=[1.5 4];
%timelims=[-50 50];

subplot(4,4,[1 9])
     %area((abs(map(1,:))).^2,omega,'FaceColor',[0.8 0.2 0.2]);
     area(omega,(abs(map(1,:))).^2,'FaceColor',[0.8 0.2 0.2]);
     xlim(freqlims)
     xlabel('Fréquence (PHz)');
     ylim([0 1.1*max((abs(map(1,:))).^2)])
     set(gca,'ytick',[])
     ylabel('Intensité (arb.u.)');
     view(-90,90)
     axis ij
     
subplot(4,4,[2 12])
    %clims=[-13 -8];
    %imagesc(spcgrmtime,omega,log(squeeze(spctrgrmmap(1,:,:))'),clims);
    imagesc(spcgrmtime,omega,squeeze(spctrgrmmap(1,:,:))');
    set(gca,'YDir','normal')
    load('MyGermanColormap','mycmap')
%    load('MyBrightJetColormap','mycmap')
    set(gcf,'Colormap',mycmap)
    ylim(freqlims)
    xlim(timelims);
    %ylabel('Fréquence (PHz)');
    %xlabel('Temps (fs)');
    title('input pulse');

subplot(4,4,[14 16])
     area(t,(abs(map_t(1,:))).^2,'FaceColor',[0.8 0.2 0.2]);
     xlim(timelims)
     xlabel('Temps (fs)');
     ylim([0 1.1*max((abs(map_t(1,:))).^2)])
     set(gca,'ytick',[])
     ylabel('Intensité (arb.u.)');
     
     if spcgrmmovie==0 && saveimgs==1 
        paperwidth=18;
        paperheight=14;
           set(gcf,'Renderer','painters',...
                  'PaperPositionMode', 'manual','PaperUnits','centimeters',...
                  'PaperSize',[paperwidth paperheight],'PaperPosition',[0.25 0.25 paperwidth-0.5 paperheight-0.5]);    

        outname = [outfolder,'/_spectrogram_input'];
        print('-r300',[outname,'.png'],'-dpng');
    end
    %%
figure(6)
%freqlims=[1.5 4];
%timelims=[-50 50];

subplot(4,4,[1 9])
     %area((abs(map(end,:))).^2,omega,'FaceColor',[0.8 0.2 0.2]);
     area(omega,(abs(map(end,:))).^2,'FaceColor',[0.8 0.2 0.2]);
     xlim(freqlims)
     xlabel('Fréquence (PHz)');
     ylim([0 1.1*max((abs(map(end,:))).^2)])
     set(gca,'ytick',[])
     ylabel('Intensité (arb.u.)');
     set(gca,'YDir','normal')
     view(-90,90)
     axis ij
     
subplot(4,4,[2 12])
    clims=[-13 -8];
    %imagesc(spcgrmtime,omega,log(squeeze(spctrgrmmap(end,:,:))'),clims);
    imagesc(spcgrmtime,omega,squeeze(spctrgrmmap(end,:,:))');
    set(gca,'YDir','normal')
    load('MyGermanColormap','mycmap')
%    load('MyBrightJetColormap','mycmap')
    set(gcf,'Colormap',mycmap)
    ylim(freqlims)
    %ylim([1.65 3.05])
    %ylabel('Fréquence (PHz)');
    xlim(timelims);
    %xlim([-100 100]);
    %xlabel('Temps (fs)');
    title('output pulse');
    
subplot(4,4,[14 16])
     area(t,(abs(map_t(end,:))).^2,'FaceColor',[0.8 0.2 0.2]);
     xlim(timelims)
     xlabel('Temps (fs)');
     ylim([0 1.1*max((abs(map_t(end,:))).^2)])
     set(gca,'ytick',[])
     ylabel('Intensité (arb.u.)');

     if spcgrmmovie==0 && saveimgs==1 
        paperwidth=18;
        paperheight=14;
           set(gcf,'Renderer','painters',...
                  'PaperPositionMode', 'manual','PaperUnits','centimeters',...
                  'PaperSize',[paperwidth paperheight],'PaperPosition',[0.25 0.25 paperwidth-0.5 paperheight-0.5]);    

        outname = [outfolder,'/_spectrogram_output'];
        print('-r300',[outname,'.png'],'-dpng');
    end

%%
% figure;
%     plot(t,real(accNLphase_t),'LineWidth',1.5);
%     xlim([-50,50]) 
%     %xlim(2*[-tau_0 tau_0]);
%     foo=max(real(accNLphase_t));
%     ylim([0,1.1*foo]) 
%     xlabel('Temps (fs)');
%     ylabel('\phi_{nl} (rad)');
%     %title(['real value of accumulated nonlinear phase']);
%         paperwidth=9;
%         paperheight=4;
%            set(gcf,'Renderer','painters',...
%                   'PaperPositionMode', 'manual','PaperUnits','centimeters',...
%                   'PaperSize',[paperwidth paperheight],'PaperPosition',[0.25 0.25 paperwidth-0.5 paperheight-0.5]);    
% 
%         %outname = [outfolder,'/real_nonlinearphase'];
%         %print('-r300',[outname,'.png'],'-dpng');
% 
%     
% figure;
%     plot(t(1:end-1),-diff(real(accNLphase_t))/dt,'LineWidth',1.5);;
%     xlim([-50,50])
%     %xlim(2*[-tau_0 tau_0]);
%     foo=max(abs(diff(real(accNLphase_t))/dt));
%     ylim(1.1*[-foo foo]);
%     xlabel('Temps (fs)');
%     ylabel('\delta\omega (PHz)');
%     %title(['frequency shift due to SPM']);
%         paperwidth=9;
%         paperheight=4;
%            set(gcf,'Renderer','painters',...
%                   'PaperPositionMode', 'manual','PaperUnits','centimeters',...
%                   'PaperSize',[paperwidth paperheight],'PaperPosition',[0.25 0.25 paperwidth-0.5 paperheight-0.5]);    
% 
%         %outname = [outfolder,'/real_nonlinearphase'];
%         %print('-r300',[outname,'.png'],'-dpng');
%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compression de l'impulsion en appliquant une GDD + TOD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

compressedfield = map(end,:) .*exp(1i*.5*k2comp*(omega-omega_0).^2 + 1i/6*k3comp*(omega-omega_0).^3);
compressedfield_t = ft(compressedfield);
compressedint =(abs(compressedfield_t)).^2;
compressedtau = diff(t(abs(diff(sign(compressedint - 0.5*(max(compressedint)))))==2));

disp(['After compression with ',num2str(k2comp),' fs^2 GDD and ',num2str(k3comp),' fs^3 TOD the pulse duration is = ',num2str(compressedtau),' fs FWHM.']);

[cmpspcgrm,spcgrmtime] = myspectrogram(compressedfield_t,t,spcgrmwinwidth,winchoice);

FTlimitint = abs(ft(abs(map(end,:)))).^2;
FTlimittau = diff(t(abs(diff(sign(FTlimitint - 0.5*(max(FTlimitint)))))==2));
disp(['The Fourier limited pulse duration is = ',num2str(FTlimittau),' fs FWHM.']);


% %%
%   figure(1)
%   hold all
%   plot(t,(abs(map_t(1,:))).^2);
%   plot(t,(abs(map_t(end,:))).^2);
%   xlabel('Temps (fs)');
%   ylabel('I(t)');
%   xlim([-60 60]);
%   ylim([0 1]);
%   
% %%
%   figure(2)
%   hold all
%   plot(omega,(abs(map(1,:))).^2);
%   plot(omega,(abs(map(end,:))).^2);
%   %plot(nu,log10((abs(map(1,:))).^2));
%   %plot(nu,log10((abs(map(end,:))).^2));
%   xlabel('Fréquence (PHz)');
%   ylabel('|I(\omega)|^2');
%   xlim([1 4]);
%   %ylim([0 .03]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(4)
    imagesc(spcgrmtime,omega,cmpspcgrm');
    %clims=[-13 -8];
    %imagesc(spcgrmtime,nu,log(cmpspcgrm'),clims);
    set(gca,'YDir','normal')
    load('MyGermanColormap','mycmap')
%    load('MyBrightJetColormap','mycmap')
    set(gcf,'Colormap',mycmap)
    ylim(freqlims)
    %ylim([1.65 3.05])
    ylabel('Fréquence (PHz)');
    xlim([-60 60]);
    xlabel('Temps (fs)');
    title(['After compression with ',num2str(k2comp),' fs^2 GDD and ',num2str(k3comp),' fs^3 TOD.']);

if saveimgs==1 
        paperwidth=18;
        paperheight=14;
           set(gcf,'Renderer','painters',...
                  'PaperPositionMode', 'manual','PaperUnits','centimeters',...
                  'PaperSize',[paperwidth paperheight],'PaperPosition',[0.25 0.25 paperwidth-0.5 paperheight-0.5]);    

        outname = [outfolder,'/_cmp-spectrogram'];
        print('-r300',[outname,'.png'],'-dpng'); 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
plot(t, (abs(compressedfield_t)).^2,'LineWidth',1.5);
 xlabel('Temps (fs)');
  ylabel('I(t)');
  %xlim([-30 30]);
  xlim([-50 50]);
  %ylim([0 2]);

if saveimgs==1 
        paperwidth=14;
        paperheight=6;
        set(gcf,'Renderer','painters',...
                  'PaperPositionMode', 'manual','PaperUnits','centimeters',...
                  'PaperSize',[paperwidth paperheight],'PaperPosition',[0.25 0.25 paperwidth-0.5 paperheight-0.5]);    

        outname = [outfolder,'/_cmp-temp-intensity'];
        print('-r300',[outname,'.png'],'-dpng'); 
end

  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
minwl=450;
%minwl=550;
maxwl=1100;
%minwl=700;
%maxwl=900;

[~,bar1]=min(abs(lambda-minwl));
[~,bar2]=min(abs(lambda-maxwl));
%[~,bar3]=min(abs(lambda-lambda_0));
cmpspecphase = unwrap(angle(compressedfield));
%cmpspecphase = cmpspecphase - min(cmpspecphase(bar2:bar1));
%cmpspecphase = cmpspecphase - cmpspecphase(bar3);
P = polyfit(omega(bar2:bar1),cmpspecphase(bar2:bar1),1);
cmpspecphase = cmpspecphase - P(1)*omega - P(2);

%[AX, H1, H2] = plotyy(nu, (abs(compressedfield)).^2, nu, unwrap(angle(compressedfield)));
%  xlabel('Fréquence (PHz)'); 
%  xlim(AX(2), [0.2 0.6]) 
%  xlim(AX(1), [0.2 0.6])  
[AX, H1, H2] = plotyy(lambda, (abs(compressedfield)).^2, lambda, cmpspecphase);
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

if saveimgs==1 
        paperwidth=16;
        paperheight=10;
        set(gcf,'Renderer','painters',...
                  'PaperPositionMode', 'manual','PaperUnits','centimeters',...
                  'PaperSize',[paperwidth paperheight],'PaperPosition',[0.25 0.25 paperwidth-0.5 paperheight-0.5]);    

        outname = [outfolder,'/_cmp-spectrum'];
        print('-r300',[outname,'.png'],'-dpng'); 
end


%%
if saveimgs==1 
        header1 = 'parameters of this simulations are:';
        fid=fopen([outfolder,'/simparameters.txt'],'wt');
        fprintf(fid, '\r\n');
        fprintf(fid, [ header1 '\r\n']);
        fprintf(fid, ['SWITCHES: \r\n']);
        fprintf(fid, '     SPM on ? = %f \r\n', SPM);
        fprintf(fid, '     SS  on ? = %f \r\n', SS);
        fprintf(fid, '     Dispersion on ? = %f \r\n', DISP);
        fprintf(fid, '     attenuation on ? = %f \r\n', ATTN);
        fprintf(fid, 'zMax = %f mm \r\n', zMax);
        fprintf(fid, 'dz = %f mm \r\n', dz);
        fprintf(fid, '\r\n');
        fprintf(fid, '# INPUT PULSE \r\n');
        fprintf(fid, 'E_0 = %f mJ \r\n', E_0);
        fprintf(fid, 'tau_0 = %f fs \r\n', tau_0);
        fprintf(fid, 'lambda_0 = %f nm \r\n', lambda_0);
        fprintf(fid, 'k2_0 = %f fs^2 \r\n',  k20);
        fprintf(fid, 'k3_0 = %f fs^3 \r\n',  k30);
        if polariz == 0
            fprintf(fid, 'linear polarization \r\n');
        elseif polariz == 1
            fprintf(fid, 'circular polarization \r\n');
        end
        fprintf(fid, 'C.Vozzis mininmal fiberdiam = %f mm (limit ionization)\r\n', 2*amin*1e3);
        fprintf(fid, '\r\n');
        fprintf(fid, '# PROPAGATION MEDIUM \r\n');
        fprintf(fid, 'filling gas = %s \r\n',  gases{gas});
        fprintf(fid, 'fiberdiam = %f mm \r\n',  fiberdiam);
        fprintf(fid, 'fiber transmission for EH11 = %10.2f percent\r\n', 100*exp(-2*attn_alpha*zMax));
        fprintf(fid, 'peak intensity at fiber input = %10.2e W/cm^2 \r\n', I_0 );
        if pZ(end) == pZ(1)
         fprintf(fid, 'constant pressure = %f bar \r\n',  p0);
        else 
         fprintf(fid, 'pressure gradient from %f to %f bar \r\n', p0,pL);
        end
        fprintf(fid, 'GVD = %f fs^2 /mm /bar \r\n',  k2);
        fprintf(fid, 'TOD = %f fs^3/mm /bar \r\n',  k3);
        fprintf(fid, 'nonlinear index /bar n2 = %e /bar \r\n', n20);
        fprintf(fid, '=> accumulated B-integral = %10.2f +i(%10.2f) rad \r\n', real(Bint), imag(Bint));
        fprintf(fid, '\r\n');
        fprintf(fid, '# COMPRESSION \r\n');
        fprintf(fid, 'compressor GDD = %f fs^2 \r\n',  k2comp);
        fprintf(fid, 'compressor TOD = %f fs^3 \r\n',  k3comp);
        fprintf(fid, '=> compressed pulse duration = %f fs FWHM \r\n.', compressedtau);
        fclose(fid);
end

 
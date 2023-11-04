%script to process acoustic data, remove instrument apparaus response and
%derive source properties seismologically
%Louisa Brotherson - 29/7/21 - also edited by Ben Edwards

fclose('all');

clear
close all

% change file directory to wwhatever folder your files are in - bath
% processes a whole folder of tests at a time
folder = 'C:\Users\lbrot\Documents\Experiments\Stick-slip\acoustics\Roughness tests\Perspex_400grit_30_1.0_03_02_22\' ; %Roughness tests\Perspex_1200grit_50_1.0_31_03_22';%Perspex_30_1.0_02_12_20\';%'C:\Users\lbrot\Documents\Experiments\Horn River\acoustics\';%
% change to folder with the extra codes in
addpath('C:\Users\lbrot\Documents\Codes')
addpath(folder)

% rough = 1 adds amplifiers in (8 channels), = 0 without(4 channels) 
wave = 'P';
rough = 1;
gain = 3;
kappa_on = 0; %fixed or unfixed kappa in Brune model (Ch 3 or ask Ben)
vp = 5.65e3; 
vs = 3.16e3;

if wave == 'P' && rough == 1
    
[excelldata,exmeta,exfilepath,exname,exext]=readin_atf_acoustic(folder,4);
instrumentapparatus = load("C:\\Users\\lbrot\\Documents\\ball_drops_BE\\ball_drops_BE_seismic\\instrument_apparatus_median.mat");
instrumentapparatus_top = load("C:\\Users\\lbrot\\Documents\\ball_drops_BE\\balldroptop\\instrument_apparatus_median.mat");
instrumentapparatus_high = load("C:\\Users\\lbrot\\Documents\\ball_drops_BE\\ball_drops_BE_seismic\\instrument_apparatus_high.mat");
instrumentapparatus_high2 = load("C:\\Users\\lbrot\\Documents\\ball_drops_BE\\ball_drops_BE_seismic\\instrument_apparatus_high2.mat");
instrumentapparatus_low = load("C:\\Users\\lbrot\\Documents\\ball_drops_BE\\ball_drops_BE_seismic\\instrument_apparatus_low.mat");
instrumentapparatus_low2 = load("C:\\Users\\lbrot\\Documents\\ball_drops_BE\\ball_drops_BE_seismic\\instrument_apparatus_low2.mat");
amp_filter_file = readmatrix('C:\Users\lbrot\Documents\Experiments\Interpolate_tests.xlsx');
amp_filter = amp_filter_file(:,gain+1);
PAS_filter = amp_filter;
THCPAS_filter = amp_filter/2;
amp_freq = amp_filter_file(:,1);
cfm = vp;
divider = 4;

elseif wave == 'P' && rough == 0
    
[excelldata,exmeta,exfilepath,exname,exext]=readin_atf_acoustic(folder,4);
instrumentapparatus = load("C:\\Users\\lbrot\\Documents\\ball_drops_BE\\ball_drops_BE_seismic\\instrument_apparatus.mat");
instrumentapparatus_high = load("C:\\Users\\lbrot\\Documents\\ball_drops_BE\\ball_drops_BE_seismic\\instrument_apparatus_high.mat");
instrumentapparatus_high2 = load("C:\\Users\\lbrot\\Documents\\ball_drops_BE\\ball_drops_BE_seismic\\instrument_apparatus_high2.mat");
instrumentapparatus_low = load("C:\\Users\\lbrot\\Documents\\ball_drops_BE\\ball_drops_BE_seismic\\instrument_apparatus_low.mat");
instrumentapparatus_low2 = load("C:\\Users\\lbrot\\Documents\\ball_drops_BE\\ball_drops_BE_seismic\\instrument_apparatus_low2.mat");
amp_filter = load(['C:\Users\lbrot\Documents\Experiments\siggen_05_07_22\4V\' num2str(gain) '\THCPAS_responses_P.mat']);
THCPAS_filter = amp_filter.THCPAS_response_gm;
PAS_filter = amp_filter.PAS_response_gm;
amp_freq = amp_filter.exfas.f_ko_snrclip{1}{1};
cfm = vp;
divider = 4;

elseif wave == 'S' && rough == 1
    
[excelldata,exmeta,exfilepath,exname,exext]=readin_atf_acoustic_S(folder,4);
instrumentapparatus = load("C:\\Users\\lbrot\\Documents\\ball_drops_BE\\ball_drops_BE_seismic_S\\instrument_apparatus_median.mat");
instrumentapparatus_top = load("C:\\Users\\lbrot\\Documents\\ball_drops_BE\\balldroptop_S\\instrument_apparatus_median.mat");
instrumentapparatus_high = load("C:\\Users\\lbrot\\Documents\\ball_drops_BE\\ball_drops_BE_seismic_S\\instrument_apparatus_high_median.mat");
instrumentapparatus_high2 = load("C:\\Users\\lbrot\\Documents\\ball_drops_BE\\ball_drops_BE_seismic_S\\instrument_apparatus_high2_median.mat");
instrumentapparatus_low = load("C:\\Users\\lbrot\\Documents\\ball_drops_BE\\ball_drops_BE_seismic_S\\instrument_apparatus_low_median.mat");
instrumentapparatus_low2 = load("C:\\Users\\lbrot\\Documents\\ball_drops_BE\\ball_drops_BE_seismic_S\\instrument_apparatus_low2_median.mat");
amp_filter_file = readmatrix('C:\Users\lbrot\Documents\Experiments\Interpolate_tests.xlsx');
amp_filter = amp_filter_file(:,gain+1);
PAS_filter = amp_filter;
THCPAS_filter = amp_filter/2;
amp_freq = amp_filter_file(:,1);
cfm = vs;
divider = 4;

elseif wave == 'S' && rough == 0
    
[excelldata,exmeta,exfilepath,exname,exext]=readin_atf_acoustic_S(folder,4);
instrumentapparatus = load("C:\\Users\\lbrot\\Documents\\ball_drops_BE\\ball_drops_BE_seismic_S\\instrument_apparatus.mat");
instrumentapparatus_high = load("C:\\Users\\lbrot\\Documents\\ball_drops_BE\\ball_drops_BE_seismic_S\\instrument_apparatus_high.mat");
instrumentapparatus_high2 = load("C:\\Users\\lbrot\\Documents\\ball_drops_BE\\ball_drops_BE_seismic_S\\instrument_apparatus_high2.mat");
instrumentapparatus_low = load("C:\\Users\\lbrot\\Documents\\ball_drops_BE\\ball_drops_BE_seismic_S\\instrument_apparatus_low.mat");
instrumentapparatus_low2 = load("C:\\Users\\lbrot\\Documents\\ball_drops_BE\\ball_drops_BE_seismic_S\\instrument_apparatus_low2.mat");
THCPAS_filter = 1;
PAS_filter = 1;
amp_freq = amp_filter.exfas.f_ko_snrclip{1}{1};
cfm = vs;
divider = 4;

else 
    
    error('Wave type (P or S) or roughness (0 or 1) entered incorrectly. Try again.')

end

instrument_apparatus = [instrumentapparatus.inst_response_gmclip;instrumentapparatus_high.inst_response_gmclip_high;instrumentapparatus_high2.inst_response_gmclip_high2;instrumentapparatus_low.inst_response_gmclip_low;instrumentapparatus_low2.inst_response_gmclip_low2]/cfm;
instrument_apparatus_top = instrumentapparatus_top.inst_response_gmclip/cfm;

%get rid of nans 

resamp_instr_freq = repmat(instrumentapparatus.inst_response_gmclip_fr,5,1);
instrument_apparatus(isnan(resamp_instr_freq))=[];
instrument_apparatus = reshape(instrument_apparatus,5,[]);

instrument_apparatus_top(isnan(resamp_instr_freq(1,:)))=[];

instrumentapparatus_top.inst_response_gmclip_fr(isnan(instrumentapparatus.inst_response_gmclip_fr))=[];
instrumentapparatus.inst_response_gmclip_fr(isnan(instrumentapparatus.inst_response_gmclip_fr))=[];

THCPAS_filter(isnan(amp_freq))=[];
PAS_filter(isnan(amp_freq))=[];
amp_freq(isnan(amp_freq))=[];

% find unique frequancies ampl_fr_r
[ampl_fr_r,ia1,~] = unique(amp_freq);

% find corresponding values (just taking amplitude of first instance)
THCPAS_ampl_r=THCPAS_filter(ia1);
PAS_ampl_r=PAS_filter(ia1);

%interpolate ampl_r to same frequencies as instr
THCPAS_ampl_i=interp1(ampl_fr_r,THCPAS_ampl_r,instrumentapparatus.inst_response_gmclip_fr);
PAS_ampl_i=interp1(ampl_fr_r,PAS_ampl_r,instrumentapparatus.inst_response_gmclip_fr);
THCPAS_ampl_i_top=interp1(ampl_fr_r,THCPAS_ampl_r,instrumentapparatus_top.inst_response_gmclip_fr);
PAS_ampl_i_top=interp1(ampl_fr_r,PAS_ampl_r,instrumentapparatus_top.inst_response_gmclip_fr);

% if you want to plot instrument-apparatus responses
plotit = 0;
if plotit == 1

    for i = 1:size(instrument_apparatus,1)

        figure
        % original (no interpolation)
        semilogy(amp_freq,THCPAS_filter, 'LineWidth', 2), 
        hold on
        plot(instrumentapparatus.inst_response_gmclip_fr,THCPAS_ampl_i, 'x', 'LineWidth', 0.5,'MarkerSize', 5)
        plot(instrumentapparatus.inst_response_gmclip_fr,instrument_apparatus(i,:), 'LineWidth', 2)
        %interpolated (both plotted against instr_fr
        plot(instrumentapparatus.inst_response_gmclip_fr,instrument_apparatus(i,:), 'x', 'LineWidth', 0.5,'MarkerSize', 5)
        %combined
        plot(instrumentapparatus.inst_response_gmclip_fr,instrument_apparatus(i,:).*THCPAS_ampl_i, 'x', 'LineWidth', 1,'MarkerSize', 5)
        xlim([0 3.5e5])
        ylim([10^-7 10^4])    
        xlabel('Frequency (Hz)')
        ylabel('Gain (-)')
        legend('PAS Response','PAS Response, interpolated','\Psi^{int}','\Psi^{int}', 'Amplified \Psi^{int}', 'Location', 'best')
        saveas(gcf, ['C:\Users\lbrot\Documents\Experiments\Stick-slip\acoustics\resample_THCPAS_amp_' num2str(i)], 'png')    

        figure
        % original (no interpolation)
        loglog(amp_freq,PAS_filter, 'LineWidth', 2), 
        hold on
        loglog(instrumentapparatus.inst_response_gmclip_fr,PAS_ampl_i, 'x', 'LineWidth', 0.5,'MarkerSize', 5)
        loglog(instrumentapparatus.inst_response_gmclip_fr,instrument_apparatus(i,:), 'LineWidth', 2)
        %interpolated (both plotted against instr_fr
        loglog(instrumentapparatus.inst_response_gmclip_fr,instrument_apparatus(i,:), 'x', 'LineWidth', 0.5, 'MarkerSize', 5)
        %combined
        loglog(instrumentapparatus.inst_response_gmclip_fr,instrument_apparatus(i,:).*PAS_ampl_i, 'x', 'LineWidth', 1,'MarkerSize', 5)
        xlim([0 4e5])
        ylim([10^-7 10^4])
        xlabel('Frequency (Hz)')
        ylabel('FAS (V)')
        saveas(gcf, ['C:\Users\lbrot\Documents\Experiments\Stick-slip\acoustics\resample_PAS_amp_' num2str(i)], 'png')           
                
        if i == 1
        figure
        % original (no interpolation)
        loglog(amp_freq,PAS_filter, 'LineWidth', 2), 
        hold on
        loglog(instrumentapparatus.inst_response_gmclip_fr,PAS_ampl_i, 'x', 'LineWidth', 0.5,'MarkerSize', 5)
        loglog(instrumentapparatus.inst_response_gmclip_fr,instrument_apparatus_top, 'LineWidth', 2)
        %interpolated (both plotted against instr_fr
        loglog(instrumentapparatus.inst_response_gmclip_fr,instrument_apparatus_top, 'x', 'LineWidth', 0.5, 'MarkerSize', 5)
        %combined
        loglog(instrumentapparatus.inst_response_gmclip_fr,instrument_apparatus_top.*PAS_ampl_i, 'x', 'LineWidth', 1,'MarkerSize', 5)
        xlim([0 4e5])
        ylim([10^-7 10^4])
        xlabel('Frequency (Hz)')
        ylabel('FAS (V)')
        saveas(gcf, 'C:\Users\lbrot\Documents\Experiments\Stick-slip\acoustics\resample_PAS_amp_top', 'png')           
        end
    end
end


%% preparing seismogram for fft and fft-ing (with KO-smoothing)

taperwidth = 0.05;
dt = 1e-7;
freqs=logspace(2,6,1024);
k=1;

for i = 1:length(excelldata)
    %convert amp to volts as per header info
    excelldata{i}.data = excelldata{i}.data.*str2num(exmeta{1,6}{2});
    excelldata{i}.data = excelldata{i}.data - mean(excelldata{i}.data);
    
    %make taper
    extaper=[sin(single(1:int64(length(excelldata{i}.data)*taperwidth/2.)).*pi/(2*length(excelldata{i}.data)*taperwidth/2.)),ones(1, int64(length(excelldata{i}.data)*(1.-taperwidth))),cos(single(1:int64(length(excelldata{i}.data)*taperwidth/2.)).*pi/(2*length(excelldata{i}.data)*taperwidth/2.))]';
    %rounding correction:
    extaper=[extaper;zeros(length(excelldata{i}.data)-length(extaper),1)];

    %compute FAS
    exfas.c{k}{i}=fft(excelldata{i}.data.*extaper);   %complex specturm

    exfas.a{k}{i}=2.*abs(exfas.c{k}{i}).*dt;        %CFAS amplitude (units = input/Hz)
    exfas.a{k}{i}(int64(length(exfas.a{k}{i})./2)+1:end)=[];

    %compute frequencies
    exfas.f{k}{i}=double(1:int64(length(exfas.a{k}{i}))).*(1./(length(exfas.a{k}{i})))*0.5/dt;
    exfas.df{k}{i}=exfas.f{k}{i}(2)-exfas.f{k}{i}(1);

    %do the K-O smoothing
    [exfas.f_ko{k}{i}, exfas.ko{k}{i}] = ko(exfas.f{k}{i},exfas.a{k}{i},min(exfas.f{k}{i}),max(exfas.f{k}{i}), 1024,40,1e-3); %1024 points, b=40, 0.1% threshold for kernel weights
    
end

    exnoisedata=excelldata;
    exnoisemeta=exmeta;
    exnoiselen=0.2;
    
    %repeat the above for the noise data (assume first 20 % of data)
    
for i = 1:length(exnoisedata)
    exnoisedata{i}.data(int64(length(exnoisedata{i}.data)*exnoiselen):end)=0;
    % remove the zeros from the nsamp
    exnoisemeta{i,3}{2}=exnoisemeta{i,3}{2}-exnoisedata{i}.data(int64(length(exnoisedata{i}.data)*exnoiselen):end);
end

%% plot, pick start/end times and save time series

for i = 1:length(exnoisedata)
    
    round_volts(i) = round(max(abs(excelldata{i}.data)),1);
    if round_volts(i) ==0
        round_volts(i) = 0.1;
    end

    figure
    plot([2:length(excelldata{i}.data)]*dt*1000,excelldata{i}.data(2:end),'-k','LineWidth',1.5)
    hold on
    plot([2:length(exnoisedata{i}.data)*exnoiselen]*dt*1000,exnoisedata{i}.data(2:length(exnoisedata{i}.data)*exnoiselen),'-r','LineWidth',1.5)
    xlabel('Time (ms)')
    ylabel('Voltage (V)')
    legend('Signal','Noise','Location','southeast')
    xlim([0 length(excelldata{i}.data)*dt*1000])
    ylim([-1.2*round_volts(i) 1.2*round_volts(i)])
    
    if wave == 'P' 
           
    saveas(gcf, [exfilepath '\\stickslip_P_' num2str(i)], 'png')    
    
    elseif wave == 'S'
        
    saveas(gcf, [exfilepath '\\stickslip_S_' num2str(i)], 'png')    
    
    end
    
    close all
    
end

%% do the same for noise so you can work out SNR and then remove instrument-apparatus responses

k=1;

for i = 1:length(exnoisedata)

    %taper 
    exntaper=[sin(single(1:int64(length(exnoisedata{i}.data)*taperwidth/2.)).*pi/(2*length(exnoisedata{i}.data)*taperwidth/2.)),ones(1, int64(length(exnoisedata{i}.data)*(1.-taperwidth))),cos(single(1:int64(length(exnoisedata{i}.data)*taperwidth/2.)).*pi/(2*length(exnoisedata{i}.data)*taperwidth/2.))]';
    %rounding correction:
    if(length(exntaper)<length(exnoisedata{i}.data))
        exntaper=[exntaper;zeros(length(exnoisedata{i}.data)-length(exntaper),1)];
    elseif(length(exntaper)>length(exnoisedata{i}.data))
        exntaper(length(exnoisedata{i}.data)+1:end)=[];
    end

    %compute FAS, frequencies, etc
    exnfas.c{k}{i}=fft(exnoisedata{i}.data.*exntaper);
    exnfas.a{k}{i}=2.*abs(exnfas.c{k}{i}).*dt;
    exnfas.a{k}{i}(int64(length(exnfas.a{k}{i})./2)+1:end)=[];
    exnfas.f{k}{i}=double(1:int64(length(exnfas.a{k}{i}))).*(1./(length(exnfas.a{k}{i})))*0.5/dt;

    %K-O smoothing (setup as per signal)
    [exnfas.f_ko{k}{i}, exnfas.ko{k}{i}] = ko (exnfas.f{k}{i},exnfas.a{k}{i},min(exnfas.f{k}{i}),max(exnfas.f{k}{i}), 1024,40,1e-3);
    exnfas.ko{k}{i}=exnfas.ko{k}{i}.*sqrt(1./exnoiselen);
    
        % snr analysis
    
    snr=3; %threshold
    exfas.ko_snrclip{k}{i}=exfas.ko{k}{i};
    %remove data and freqs where SNR<snr
    exfas.ko_snrclip{k}{i}(exfas.ko{k}{i}./exnfas.ko{k}{i}<snr)=[];
    exfas.f_ko_snrclip{k}{i}=exfas.f_ko{k}{i};
    exfas.f_ko_snrclip{k}{i}(exfas.ko{k}{i}./exnfas.ko{k}{i}<snr)=[];
    
end
    
for p = 1:length(exnoisedata)
    k = 1;
    
    exfas.ko_snrclip_cor{k}{p}=zeros(1,length(freqs));
    exfas.ko_snrclip_cor_fr{k}{p}=zeros(1,length(freqs));
    exfas.ko_snrclip_cor_n{k}{p}=zeros(1,length(freqs));
    

    for i = 1:length(freqs)-1
        exfas.ko_snrclip_cor{k}{p}(i)=exfas.ko_snrclip_cor{k}{p}(i)+sum(log(exfas.ko_snrclip{k}{p}(exfas.f_ko_snrclip{k}{p}>=freqs(i)&exfas.f_ko_snrclip{k}{p}<freqs(i+1))));
        exfas.ko_snrclip_cor_fr{k}{p}(i)=exfas.ko_snrclip_cor_fr{k}{p}(i)+sum(log(exfas.f_ko_snrclip{k}{p}(exfas.f_ko_snrclip{k}{p}>=freqs(i)&exfas.f_ko_snrclip{k}{p}<freqs(i+1))));
        exfas.ko_snrclip_cor_n{k}{p}(i)=exfas.ko_snrclip_cor_n{k}{p}(i)+length( log(exfas.f_ko_snrclip{k}{p}(exfas.f_ko_snrclip{k}{p}>=freqs(i)&exfas.f_ko_snrclip{k}{p}<freqs(i+1))));
        
        
    end

    exfas.ko_snrclip_cor{k}{p}=exfas.ko_snrclip_cor{k}{p}./exfas.ko_snrclip_cor_n{k}{p};
    exfas.ko_snrclip_cor{k}{p}=exp(exfas.ko_snrclip_cor{k}{p});
    exfas.ko_snrclip_cor{k}{p}(exfas.ko_snrclip_cor{k}{p}==0)=NaN;
    exfas.ko_snrclip_cor{k}{p}(isnan(exfas.ko_snrclip_cor{k}{p}))=NaN;

    exfas.ko_snrclip_cor_fr{k}{p}=exfas.ko_snrclip_cor_fr{k}{p}./exfas.ko_snrclip_cor_n{k}{p};
    exfas.ko_snrclip_cor_fr{k}{p}=exp(exfas.ko_snrclip_cor_fr{k}{p});
    
     % find real spectra
     
     if k == 1
         
          exfas.ko_snrclip_cor_hold{k}{p} = exfas.ko_snrclip_cor{1}{p};
          
     end        
     
     % interpolate so frequencies are the same as instrument apparatus
     % response - LB 31/10/22   
     
     % remove NaNs 
     exfas.ko_snrclip_cor_hold{k}{p}(isnan(exfas.ko_snrclip_cor_fr{k}{p}))=[];
     exfas.ko_snrclip_cor_fr{k}{p}(isnan(exfas.ko_snrclip_cor_fr{k}{p}))=[];
     
     % do unique for noise vector
    [exnfas.ko_r{k}{p},ia2,~] = unique(exnfas.f_ko{k}{p});     
    exnfas.f_ko_r{k}{p} = exnfas.ko{k}{p}(ia2);
     % interpolate all - eliminate short vectors
     
     if length(exfas.ko_snrclip_cor_fr{k}{p})>50
         
         exfas.ko_snrclip_cor_i{k}{p} = interp1(exfas.ko_snrclip_cor_fr{k}{p},exfas.ko_snrclip_cor_hold{k}{p},instrumentapparatus.inst_response_gmclip_fr);
         exnfas.ko_i{k}{p} = interp1(exnfas.ko_r{k}{p},exnfas.f_ko_r{k}{p},instrumentapparatus.inst_response_gmclip_fr);     
         
     else
         
         exfas.ko_snrclip_cor_i{k}{p} = zeros(1,length(instrumentapparatus.inst_response_gmclip_fr));
         exnfas.ko_i{k}{p} = zeros(1,length(instrumentapparatus.inst_response_gmclip_fr));
         
     end 
     
     % if you change channels for THC etc. change these!

     for k = 1:size(instrument_apparatus,1)
         
         if p <= length(exnoisedata)/divider
             
            exfas.ko_snrclip_cor{k}{p}=exfas.ko_snrclip_cor_i{1}{p}./(THCPAS_ampl_i.*instrument_apparatus(k,:));
            exnfas.ko_snrclip_cor{k}{p}=exnfas.ko_i{1}{p}./(THCPAS_ampl_i.*instrument_apparatus(k,:));

                       
         elseif p > length(exnoisedata)/divider && p <= (length(exnoisedata))*(2/divider)
             
            exfas.ko_snrclip_cor{k}{p}=exfas.ko_snrclip_cor_i{1}{p}./(PAS_ampl_i.*instrument_apparatus(k,:)); 
            exnfas.ko_snrclip_cor{k}{p}=exnfas.ko_i{1}{p}./(PAS_ampl_i.*instrument_apparatus(k,:));
            
         elseif p > length(exnoisedata)*(2/divider) && p <= length(exnoisedata)*(3/divider)
             
            exfas.ko_snrclip_cor{k}{p}=exfas.ko_snrclip_cor_i{1}{p}./instrument_apparatus(k,:);   
            exnfas.ko_snrclip_cor{k}{p}=exnfas.ko_i{1}{p}./instrument_apparatus(k,:);   
            
         else
             
            exfas.ko_snrclip_cor{k}{p}=exfas.ko_snrclip_cor_i{1}{p}./instrument_apparatus_top;
            exnfas.ko_snrclip_cor{k}{p}=exnfas.ko_i{1}{p}./instrument_apparatus_top;
            
         end
      
     end 
       
end   
%% brune modelling and saving to file

clear stickslip_data stickslip_ndata m_file m_re

for i = 1:length(exnoisedata)
    
    for k = 1:size(instrument_apparatus,1)

    stickslip_data{k}{i}(:,1) = instrumentapparatus.inst_response_gmclip_fr;
    stickslip_data{k}{i}(:,2) = exfas.ko_snrclip_cor{k}{i};
    stickslip_data{k}{i}(any(isnan(stickslip_data{k}{i}), 2), :) = [];

    stickslip_ndata{k}{i}(:,1) = instrumentapparatus.inst_response_gmclip_fr;
    stickslip_ndata{k}{i}(:,2) = exnfas.ko_snrclip_cor{k}{i};
    stickslip_ndata{k}{i}(any(isnan(stickslip_ndata{k}{i}), 2), :) = [];    
    
    end
    
end

for i = 1:length(exnoisedata)
    
    for k = 1:size(instrument_apparatus,1)    
    
   figure
   
   if max(excelldata{i}.data)>0.2 && isempty(stickslip_data{k}{i})==0 && length(exfas.ko_snrclip_cor_fr{1}{i})>50
    
   m{k}{i} = brune_fit_LB(stickslip_data{k}{i}(:,1),stickslip_data{k}{i}(:,2),min(exfas.ko_snrclip_cor_fr{1}{i}(isfinite(exfas.ko_snrclip_cor_fr{1}{i}))), max(exfas.ko_snrclip_cor_fr{1}{i}(isfinite(exfas.ko_snrclip_cor_fr{1}{i}))),stickslip_ndata{1}{i}(:,1),stickslip_ndata{1}{i}(:,2),kappa_on); 
   
    if wave == 'P' && kappa_on == 1
        
        saveas(gcf, [exfilepath '\\brune_P_kappa_' num2str(k) '_' num2str(i)], 'png')
        
    elseif wave == 'P' && kappa_on == 0
        
        saveas(gcf, [exfilepath '\\brune_P_' num2str(k) '_' num2str(i)], 'png')
        
    elseif wave == 'S' && kappa_on == 1
        
        saveas(gcf, [exfilepath '\\brune_S_kappa_' num2str(k) '_' num2str(i)], 'png')        
      
    elseif wave == 'S' && kappa_on == 0
        
        saveas(gcf, [exfilepath '\\brune_S_' num2str(k) '_' num2str(i)], 'png')
    
    end
    
   else
       
     m{k}{i} = [0 0 0];
     
   end
    
    end
    
    close all   
    
    m_re{i} = cell2mat(m{1}(i));
    
end

%% reshape files

for k = 1:size(instrument_apparatus,1) 

    m_file_hold(k,:) = cell2mat(m{k});
    m_file{k} = reshape(m_file_hold(k,:),3,[]); 
    m_file{k} = m_file{k}';

end

m_file = m_file';

%% saving files

for k = 1:size(instrument_apparatus,1)    


    if wave == 'P' && kappa_on == 1

        fid = fopen([exfilepath '\sourceprop_P_kappa' num2str(k) '.txt'],'w');
        fprintf(fid,'%1.5d %1.5d %1.5d \n',m_file{k}');
        fclose(fid);
        save([exfilepath '\acoustic_process_P.mat'])    
        
    elseif wave == 'P' && kappa_on == 0
        
        fid = fopen([exfilepath '\sourceprop_P' num2str(k) '.txt'],'w');
        fprintf(fid,'%1.5d %1.5d %1.5d \n',m_file{k}');
        fclose(fid);
        save([exfilepath '\acoustic_process_P.mat'])           

    elseif wave == 'S' && kappa_on == 1

        fid = fopen([exfilepath '\sourceprop_S_kappa' num2str(k) '.txt'],'w');
        fprintf(fid,'%1.5d %1.5d %1.5d \n',m_file{k}');
        fclose(fid);
        save([exfilepath '\acoustic_process_S.mat'])
        
    elseif wave == 'S' && kappa_on == 0

        fid = fopen([exfilepath '\sourceprop_S' num2str(k) '.txt'],'w');
        fprintf(fid,'%1.5d %1.5d %1.5d \n',m_file{k}');
        fclose(fid);
        save([exfilepath '\acoustic_process_S.mat'])               

    end
    
end
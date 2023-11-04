function [f_ko, fas_kosmooth] = ko (f,fas,fmin,fmax, nsamp_ko,b,ko_threshold)

%
% Konno-Ohmachi smoothing
%

%b=30; % b of the Konno-Ohmachi smoothing
%nsamp_ko; % number of central frequencies to use (ie number of smoothed fas to generate)
%fmin=0.1; % minimum frequency to generate
%fmax=50; % max frequency to generate
%ko_threshold=1e-3; % threshold below which to assume kernel weights = 0
%                     1e-2 = 1%

if(size(fas,2)==1)
    fas=fas';
end
if(size(fas,1)~=1)
    disp 'Error in input signal - provide single column data'
    exit;
end

[~,fmin_loc] = min(abs(f-fmin)); % get location of min frequency (vector index of f)
[NULL,fmax_loc] = min(abs(f-fmax)); % get location of max frequency (vector index of f)

for h=1:nsamp_ko % loop over the nsamp_ko samples 
    f_target=exp(log(f(fmin_loc))+h*(log(f(fmax_loc))-log(f(fmin_loc)))/nsamp_ko); % calc the target frequency (from log-space nsamp_ko samples)
    [NULL,i] = min(abs(f-f_target)); % f_target's index in f
    f0=f(i); % the central frequency for the KO weights
    ko_weights=(sin(b.*log10(f./f0))./(b.*log10(f./f0))).^4; % KO weights
    ko_weights(i)=1; % KO-weight = 1 where frequecy == central frequency
    
    %speed up - threshold for weight > 1e-3
    
    
    fas_tmp=fas;
    fas_tmp(ko_weights<1e-6)=[];
    ko_weights(ko_weights<1e-6)=[];
    if(length(ko_weights)<1)
        fas_kosmooth(h)=0;
        f_ko(h)=f(i);
        continue;
    end
    
    fas_kosmooth(h)=exp(sum(log(fas_tmp(1:length(fas_tmp))).*ko_weights)./sum(ko_weights)); % apply the weights - note we apply to the log(fas), then return to the fas with exp (different to the spreadsheet - but better as ratios are normal distributions in log-space)
    %fas_kosmooth(h)=(sum((fas(1:length(ew)./2)).*ko_weights)./sum(ko_weights)); % this would be the same in the lin-space (as the spreadsheet example)
    f_ko(h)=f(i);
  
end
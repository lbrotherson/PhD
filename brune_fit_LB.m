
function [m] = brune_fit_LB(f, fas, flow, fhigh, nf, nfas)

% function to fit the Brune function to displacement FAS
% flow and fhigh are the fitting limit
% m(1) = M0*(long-period level)
% m(2) = f0
% m(3) = kappa
% Edited by Louisa Brotherson 24/06/21

r=log(fas);
nr = log(nfas);
loglog(f,exp(r),'r','LineWidth',1.5)                         % plot the input data
hold on
loglog(nf(nf>flow&nf<=f(end)),exp(nr(nf>flow&nf<=f(end))),'k','LineWidth',1.5)   

% this is the ratio of mainshock/aftershock
brune=@(m)(log((m(1).*(1./(1+(f(f>flow&f<fhigh)./m(2)).^2)).*exp(-pi.*f(f>flow&f<fhigh).*m(3)))));
% this is the misfit between input data (r) and egf_ratio
misfit=@(m)(r(f>flow&f<fhigh)-brune(m));

m0=[1e-2,1,0.1];                        %set the starting values for m
mlb=[1e-12 0.01 0.001];             %set the lower bound values for m
mub=[1 100 5];            %set the upper bound values for m
    
%plot starting model
hold on ; plot(f(f>flow&f<fhigh),exp(brune(m0)),'--b','LineWidth',1.5);

m=lsqnonlin(misfit,m0,mlb,mub);         % run the non-linear optimisation

plot(f(f>flow&f<fhigh),exp(brune(m)),'b','LineWidth',1.5);

xlabel('Frequency (Hz)')
ylabel('FAS (Nm)')
xlim([10^2 10^6])

legend('signal','noise','starting model (m_{0})', 'inverted model (m)','Location','eastoutside')

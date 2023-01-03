function [M1,M2,M3,stickslip_ind_start_new3,stickslip_ind_end,breaking,stickslip_ind_start_new4] = find_stickslip_v4(in,slength,smin,stol,llimit,ulimit,tmin)
%FIND_STICKSLIP finds indices of start of stickslip
%   Matlab function to find indice of start of stickslip and stress drop
%   for plotting and processing
%   Louisa Brotherson - 25/10/2021

tau = in(:,2); % stress in MPa
disp = in(:,1);

M1 = movmean(tau,slength);
M2 = movmean(tau,slength*4);
M3 = movmean(tau,slength*16);
M4 = movmean(tau,3);

% STA/LTA, excluding those outside of tolerance

for i = 2:length(tau)-1
    
    if tau(i+1)<=tau(i) && tau(i)>=tau(i-1) && tau(i)>=M1(i)*stol && tau(i)>=M2(i)*stol && M1(i)>=M2(i) && disp(i)>=llimit && disp(i)<=ulimit && tau(i)>=smin && tau(i)>=mean(tau)%&& axial(i+1)-axial(i)<=-1*smin%%(axial(i+1)-axial(i))/(disp(i+1)-disp(i))<=-1 && (axial(i)-axial(i-1))/(disp(i)-disp(i-1))>=1 && %%|| (axial(i+1)-axial(i))/(disp(i+1)-disp(i))>=20 || (axial(i)-axial(i-1))/(disp(i)-disp(i-1))<=-20%&& (M1(i+ceil(slength/2))-M1(i))/(disp(i+ceil(slength/2))-disp(i))<=0%(M2(i+1)-M2(i))/(disp(i+1)-disp(i))<=-1 && (M2(i+2)-M2(i+1))/(disp(i+2)-disp(i+1))>=200 && axial(i+1)>=axial(i) && axial(i+2)<=axial(i+1) && (M1(i)-M1(i-1))+smin <= 0 && (M1(i)-M1(i-2))+smin <= 0 
        
        stickslip_ind_start(i) = i;
        
    elseif tau(i+1)>=tau(i) && tau(i)>=tau(i-1) && tau(i)>=M1(i)*stol&& tau(i)>=M2(i)*stol && M1(i)>=M2(i)&& disp(i)>=llimit && disp(i)<=ulimit && tau(i)>=smin && tau(i)>=mean(tau)
        
        stickslip_ind_start(i) = i; 
 
    elseif tau(i+1)<=tau(i) && tau(i)<=tau(i-1) && tau(i)>=M1(i)*stol&& tau(i)>=M2(i)*stol && M1(i)>=M2(i)&& disp(i)>=llimit && disp(i)<=ulimit && tau(i)>=smin && tau(i)>=mean(tau)
        
        stickslip_ind_start(i) = i; 
        
    elseif tau(i+1)==tau(i) && tau(i)>=M1(i)*stol&& tau(i)>=M2(i)*stol && disp(i)>=llimit && disp(i)<=ulimit && tau(i)>=smin && tau(i)>=mean(tau)
        
        stickslip_ind_start(i) = i; 
        
    else
        
        stickslip_ind_start(i) = 0;
        
    end
    
    
end

stickslip_ind_start(stickslip_ind_start==0)=[];

% discarding inflection points (maxima only)

for i = 2:length(stickslip_ind_start)-1
    
    if M2(stickslip_ind_start(i+1))*1.000001 >= M2(stickslip_ind_start(i)) || M2(stickslip_ind_start(i-1))*1.000001 >= M2(stickslip_ind_start(i))
        
        stickslip_ind_start_new(i) = 0;
        
    else
        
         stickslip_ind_start_new(i) = stickslip_ind_start(i);
         
    end
     
    
end

stickslip_ind_start_new(stickslip_ind_start_new==0)=[];

%remove outliers (false maxima)

for i = 1:length(stickslip_ind_start_new)-1
    
    if M3(stickslip_ind_start_new(i)+slength*16)>= M3(stickslip_ind_start_new(i))
        
        stickslip_ind_start_new2(i) = 0;
        
    else
        
         stickslip_ind_start_new2(i) = stickslip_ind_start_new(i);
         
    end 
    
end

stickslip_ind_start_new2(stickslip_ind_start_new2==0)=[];

stickslip_ind_start_new4 = stickslip_ind_start_new2;

%get rid of duplicates

for i = 2:length(stickslip_ind_start_new2)
    
    if i == 2
        
        stickslip_ind_start_new3(i-1) = stickslip_ind_start_new2(1);        
    
    elseif stickslip_ind_start_new2(i)-stickslip_ind_start_new2(i-1)>=110
        
        stickslip_ind_start_new3(i-1) = stickslip_ind_start_new2(i-1);
              
    else
        
        stickslip_ind_start_new3(i-1) = 0;        
        
    end    
    
end

stickslip_ind_start_new3(stickslip_ind_start_new3==0)=[];

%find the ends of stickslips

for i = 1:length(stickslip_ind_start_new3)
    
    for j = 1:1000
        
        stickslip_ind_end(i) = stickslip_ind_start_new3(i)+j;
        
        if  M1(stickslip_ind_end(i)-1) >= M1(stickslip_ind_end(i))*1.001 && M1(stickslip_ind_end(i)+1) >= M1(stickslip_ind_end(i))*1.001 && M1(stickslip_ind_end(i)-2) >= M1(stickslip_ind_end(i))*1.001 && M1(stickslip_ind_end(i)+2) >= M1(stickslip_ind_end(i))*1.001 && M1(stickslip_ind_end(i)-3) >= M1(stickslip_ind_end(i))*1.001 && M1(stickslip_ind_end(i)+3) >= M1(stickslip_ind_end(i))*1.001 && M1(stickslip_ind_end(i)-4) >= M1(stickslip_ind_end(i))*1.001 && M1(stickslip_ind_end(i)+4) >= M1(stickslip_ind_end(i))*1.001 && tau(stickslip_ind_end(i)+50) >= tau(stickslip_ind_end(i)) %tau(stickslip_ind_end(i)+1) >= tau(stickslip_ind_end(i))*1.0000001 %&& && M1(stickslip_ind_end(i)-2) >= M1(stickslip_ind_end(i))*1.0000001 && M1(stickslip_ind_end(i)+2) >= M1(stickslip_ind_end(i))*1.0000001% && tau(stickslip_ind_end(i)-1) >= tau(stickslip_ind_end(i))*1.0001 
            
            breaking{i} = j;
            
            break 
            
        else
            
            stickslip_ind_end(i) = stickslip_ind_start_new3(i)+200;           
                               
        end
        
        
    end
    
    
end

for i = 1:length(stickslip_ind_end)
    
    dtau(i) = tau(stickslip_ind_start_new3(i))-tau(stickslip_ind_end(i));
    
    if dtau(i) <= tmin
        
        stickslip_ind_start_new3(i) = 0;
        stickslip_ind_end(i) = 0;
        
    end
    
    
end

stickslip_ind_start_new3(stickslip_ind_start_new3==0)=[];
stickslip_ind_end(stickslip_ind_end==0)=[];


end
function HB=calcHB(SpO2,RespEvents,SleepStage,PLOTFLAG)
%
%Calculates hypoxic burden
%21 July 2022
%Author: Philip de Chazal
%Version. 1.0
%This code was used in 
%Kate Sutherland, Nadi Sadr, Yu Sun Bin, Kristina Cook, Hasthi U. Dissanayake, Peter A. Cistulli and Philip de Chazal
%"Comparative associations of oximetry patterns in Obstructive Sleep Apnea with incident cardiovascular disease", Sleep, 2022
%
%Inputs
%------
%
%SpO2 = 
%struct with fields:
%   Sig: [L×1 double] %Sample values of the SpO2 signal
%     SR: 1 %Sample rate in Hertz
%
%RespEvents = 
%struct with fields:
%            Type: {1×M cell} %string code describing event type e.g. 'H','OA','C' 
%            Start: [1×M double] %Time in seconds from the start of the recording for the start of the event
%         Duration: [1×M double] %Duration of the event in seconds
%
%SleepStage = 
%struct with fields:
%    Annotation: [1×N double] %The code of the sleep stage annotations e.g. [0,0,1,2,1,2,....] 
%          Codes: [0 1 2 3 4 5 9] %List of all possible codes
%    Description: {'Wake','Stage 1','Stage 2','Stage 3','Stage 4','REM','Indeterminant'} %List of explantaion of codes
%             SR: 1 %Sample rate of sleep stage annotations
%
%PLOTFLAG: 0-no plots, 1-plots
%
%
%Outputs
%-------
%
%HB: Hypoxic burden value. units are percent desaturation minutes per hour of sleep. NaN if no value can be calculated
%

%To do %%%%%%
%1. Test that code is independent of SpO2 and SleepStage sample rate (currently only accepts 1 Hz)
%
%Known issues
%The average event gap should be the average offset to next event onset gap. Currently it is average time between event onsets 

if nargin==3
    PLOTFLAG=false;
end

if SpO2.SR~=1 || SleepStage.SR~=1
    error('SpO2 signal and Sleep stage annotations need to sampled at 1Hz.Code is setup for other sample rates but not tested.Particularly need to check filter code!')
end

if isempty(RespEvents.Type)
    HB=nan; % (percent desaturation minutes of SDB) per (hour sleep)
else
    %Set physiologically unreasonable SpO2 values to NaN
    SpO2phys=SpO2.Sig;
    SpO2phys(SpO2phys<50)=nan;
    SpO2phys(SpO2phys>100)=nan;
    
    NumEvents=length(RespEvents.Type); %Process all respiratory events
    LenSpO2=length(SpO2phys);
    
    
%%%%%%%%%%%
%HB calculation
%%%%%%%%%%%%%
%Steps
%   
% 1.	Determine the timing of the average event from the event files
%       a.	Calculate the average event duration (DurAvg)
%       b.	Calculate the average event gap (AvgEventGap)
% 2.	Overlay the oxygen saturation signals (SpO2) associated with all respiratory events. 
%       The signals are synchronized at the termination of the respiratory events (TimeZero) 
%       and the sampling window includes the SpO2 signal from TimeZero -120 seconds to TimeZero +120 seconds. 
%       a.	They are averaged to form the ensemble averaged SpO2 curve.
%       b.	The ensemble averaged SpO2 curve is filtered with a 0.03Hz low pass filter to form the filtered ensemble averaged SpO2 curve.
%       c.	The filtered averaged SpO2 signal is truncated to span the average onset point to the minimum of the average onset of the next event and 90 seconds. 
%           It is truncated to TimeZero-DurAvg to TimeZero + the minimum of 90 seconds and AvgEventGap. The resulting signal is referred to as the SpO2 response.
% 3.	Determine the start and end point of the search window from the SpO2 response.
%       a.	Find minimum point of SpO2 response (Nadir)
%       b.	Find maximum difference between start of truncated averaged SpO2 signal and Nadir (MaxDesatOnset). 
%       c.	Find last peak at least 75% of amplitude of MaxDesatOnset before the time occurrence of Nadir. This is the start point of the search window (WinStart). 
%       d.	Find maximum difference between Nadir and the end of the SpO2 response (MaxDesatOffset). 
%       e.	Find first peak at least 75% of amplitude of MaxDesatOnset after the time occurrence of Nadir. This is the end point of the search window (WinFinish). 
% 4.	For each event do the following
%       a.	Find the pre-event baseline saturation which is defined as the maximum SpO2 during the 100 seconds prior to the end of the event.
%       b.	Find the area between pre-event baseline, the SpO2 curve, and WinStart and WinEnd of the search window.
%           i.	If any of the SpO2 curve is above the pre-event baseline, then do not add this negative area
%           ii.	If event search window overlaps the next event, then do not add the area twice.
% 5.	The Hypoxic Burden is defined as the sum event areas divided by the sleep time and has units of %minutes per hour.
    
 
    %1a. Calculate the average event duration
    DurAvg=ceil(nanmean(RespEvents.Duration(:)));
    
    %1b. Calculate the average event gap
    AvgEventGap=ceil(nanmean(diff(RespEvents.Start(:))));

    %2a. SpO2 signals for each event are averaged to form the ensemble averaged SpO2 curve.
    SpO2event=nan(240*SpO2.SR+1,NumEvents);
    for event=1:NumEvents
        Finish=round((RespEvents.Start(event) +...
            RespEvents.Duration(event))*SpO2.SR+0.5);
        if Finish-120*SpO2.SR>0 && Finish+120*SpO2.SR<LenSpO2
            SpO2event(:,event)=SpO2phys(Finish-120*SpO2.SR:Finish+120*SpO2.SR); %percent minutes
        end
    end
    SpO2avg=nanmean(SpO2event,2);
    
    %2b. The ensemble averaged SpO2 curve is filtered with a 0.03Hz low pass filter to form the filtered ensemble averaged SpO2 curve.
    if SpO2.SR==1
        %     %Code to design the FIR lowpass filter with cutoff at 1/30Hz
        %     N   = 30;        % FIR filter order
        %     Fp  = 1/30;       % 1/30z passband-edge frequency
        %     Fs  = 1;       % 1Hz sampling frequency
        %     Rp  = 0.00057565; % Corresponds to 0.01 dB peak-to-peak ripple
        %     Rst = 1e-5;       % Corresponds to 100 dB stopband attenuation
        %
        %     B = firceqrip(N,Fp/(Fs/2),[Rp Rst],'passedge'); % eqnum = vec of coeffs
        %     %fvtool(B,'Fs',Fs,'Color','White') % Visualize filter
        
        
        B=[0.000109398212241,   0.000514594526374,   0.001350397179936,   0.002341700062534,...
            0.002485940327008,   0.000207543145171,  -0.005659450344228,  -0.014258087808069,...
            -0.021415481383353,  -0.019969417749860,  -0.002425120103463,   0.034794452821365,...
            0.087695691366900,   0.144171828095816,   0.187717212244959,   0.204101948813338,...
            0.187717212244959,   0.144171828095816,   0.087695691366900,   0.034794452821365,...
            -0.002425120103463,  -0.019969417749860,  -0.021415481383353,  -0.014258087808069,...
            -0.005659450344228,   0.000207543145171,   0.002485940327008,   0.002341700062534,...
            0.001350397179936,   0.000514594526374,   0.000109398212241];
    else
        N   = 30*SpO2.SR;        % FIR filter order
        Fp  = 1/30;       % 1/30z passband-edge frequency
        Fs  = SpO2.SR;       % sampling frequency
        Rp  = 0.00057565; % Corresponds to 0.01 dB peak-to-peak ripple
        Rst = 1e-5;       % Corresponds to 100 dB stopband attenuation
        
        B = firceqrip(N,Fp/(Fs/2),[Rp Rst],'passedge'); % eqnum = vec of coeffs
        %fvtool(B,'Fs',Fs,'Color','White') % Visualize filter
    end
    SpO2avgfilt=filtfilt(B,1,SpO2avg);
    
    %The filtered averaged SpO2 signal is truncated to span the average onset point to the minimum of the average event gap and 90 seconds
 
    s=120*SpO2.SR+1-DurAvg*SpO2.SR;
    f=120*SpO2.SR+1+min(90,AvgEventGap)*SpO2.SR; %Limit to 90 seconds
    m=DurAvg;
    
    SpO2avgfilt=SpO2avgfilt(s:f);
    
    %%%%
    %Plots
    %%%%%%%
    if PLOTFLAG
        clf
        plot(SpO2avgfilt)
        hold on
        plot(SpO2avg(s:f))
        plot(m,SpO2avgfilt(m),'b*')
    end
    
    NadirIdx=[]; WinStart=[]; WinFinish=[];
    [Nadir,NadirIdx]=findpeaks(-SpO2avgfilt);
    if ~isempty(Nadir)
        [~,idx]=max(Nadir);
        NadirIdx=NadirIdx(idx);
        %3a. Find minimum point of SpO2 response (Nadir)
        Nadir=-Nadir(idx);
        if NadirIdx>=3 && NadirIdx<=length(SpO2avgfilt)-2
            
            if PLOTFLAG
                plot(NadirIdx,SpO2avgfilt(NadirIdx),'bv')
            end
            
             %3b.	Find maximum difference between start of truncated averaged SpO2 signal and Nadir (MaxDesatOnset).
            [peaks,ind]=findpeaks(SpO2avgfilt(1:NadirIdx));
            if ~isempty(ind)
                MaxDesatOnset=max(peaks);
                %3c.	Find last peak at least 75% of amplitude of MaxDesatOnset before the time occurrence of Nadir. This is the start point of the search window (WinStart).
                idx=find(peaks-Nadir > 0.75*(MaxDesatOnset-Nadir));
                WinStart=ind(idx(end));
                if PLOTFLAG
                    plot(WinStart,SpO2avgfilt(WinStart),'b<')
                end
                %3d.	Find maximum difference between Nadir and the end of the SpO2 response (MaxDesatOffset).
                [peaks,ind]=findpeaks(SpO2avgfilt(NadirIdx:end));
                
                if ~isempty(ind)
                    MaxDesatOffset=max(peaks);
                    %3e.	Find first peak at least 75% of amplitude of MaxDesatOnset after the time occurrence of Nadir. This is the end point of the search window (WinFinish).
                    idx=find(peaks-Nadir > 0.75*(MaxDesatOffset-Nadir));
                    WinFinish=ind(idx(1))+NadirIdx-1;
                    if PLOTFLAG
                        plot(WinFinish,SpO2avgfilt(WinFinish),'b>')
                    end
                end
            end
        end
    end
    
     if isempty(NadirIdx) || isempty(WinStart) || isempty(WinFinish)
        disp('used population defaults')
        %If unable to find window start and/or end points use defaults
        %Defaults determined from average event response from all SHHS recordings
        WinStart=m-5*SpO2.SR; %Window start is 5 seconds before event end
        WinFinish=m+45*SpO2.SR; %Window finish is 45 seconds after event end
    end
    WinStart=WinStart-m;
    WinFinish=WinFinish-m;
    
    PercentMinsDesat=0; 
    Limit=0; %Prevents double counting of areas when events are within window width of each other
    
    for event=1:NumEvents
        Finish=round((RespEvents.Start(event) +...
            RespEvents.Duration(event))*SpO2.SR +0.5);
        
        if Finish-100*SpO2.SR>0 && Finish+WinFinish<LenSpO2
            
            %Double count and negative area correction
            %       4a.	Find the pre-event baseline saturation which is defined as the maximum SpO2 during the 100 seconds prior to the end of the event.
            %       4b.	Find the area between pre-event baseline, the SpO2 curve, and WinStart and WinEnd of the search window.
            %           4bi.	If any of the SpO2 curve is above the pre-event baseline, then do not add this negative area
            %           4bii.	If event search window overlaps the next event, then do not add the area twice.
            PercentMinsDesat=PercentMinsDesat +  nansum( max(nanmax(SpO2phys(Finish-100*SpO2.SR:Finish))-SpO2phys(max(Finish+WinStart,Limit):Finish+WinFinish),0) )/(60*SpO2.SR); %percent minutes
            Limit=Finish+WinFinish;
        end
    end
    
    HourSleep=sum(~isnan(SpO2phys(SleepStage.Annotation>0 & SleepStage.Annotation<9)))/(3600*SleepStage.SR); %hours
    
    % 5.	The Hypoxic Burden is defined as the sum event areas divided by the sleep time and has units of %minutes per hour.
    HB=PercentMinsDesat/HourSleep; % (percent desaturation minutes of SDB) per (hour sleep)
    
    if PLOTFLAG
        pause(0.1);
        title(sprintf('Blue trace: ensembled SpO2 signal\n',...
                       'Red trace: filtered ensembled SpO2 signal\n', ...
                       '<|: window start\n', ...
                       '|>:window end\n', ... 
                       '*:Time Zero\n', ...
                       'v: Nadir point'))
        title(sprintf('Red trace: ensembled SpO2 signal\nBlue trace: filtered ensembled SpO2 signal\n<|: search window start, |>: search window end\n*: Time Zero (average event offset point), v: Ensembled nadir point'))
        ylabel('Percent saturation')
        xlabel('Time from ensembled event onset point (seconds)')
    end
end

end
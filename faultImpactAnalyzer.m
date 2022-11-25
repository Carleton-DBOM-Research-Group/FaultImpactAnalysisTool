function [energy, znMdlPrmtr] = faultImpactAnalyzer(seasonalAvailability,sOaMin,dataDirectory,start,stop,wkndOp,htgSp,clgSp,htgSpUnocc,clgSpUnocc)
    
    %   author: Burak Gunay [burak.gunay@carleton.ca]
    %   date: October 31, 2022
    
    %   faultImpactAnalyzer function estimates the existing and corrected
    %   (fault-free) heating and cooling load of a multiple zone VAV AHU system
    
    %   Function calculates these based on air-side measurements of BAS
    %   trend data
      
    
    %   inputs:

    %   seasonalAvailability: 0 means heating and cooling available
    %   year-round. 1 means heating is available from October to April and
    %   cooling is available from May to September
    
    %   sOaMin: minimum outdoor air damper position setpoint (typically 0.1
    %   to 0.5)
    
    %   dataDirectory: the folder in which data from a multiple zone 
    %   VAV AHU system is located
    %   example: C:\Users\burak\Dropbox\New folder\data
    
    %   start: morning start time for the AHU fans (typically 4 to 8)
    
    %   stop: evening stop time for the AHU fans (typically 16 to 20)
    
    %   wkndOp: 1 if AHU system is available in the weekends or 0 if the
    %   system is off in the weekends.
    
    %   htgSp: heating temperature setpoint (typically 20 to 23 degC)
    
    %   clgSp: cooling temperature setpoint (typically 22 to 25 degC)
    
    %   htgSpUnocc: unoccupied mode heating setpoint (typically 15 to 20
    %   degC)
    
    %   clgSpUnocc: unoccupied mode cooling setpoint (typically 25 to 30
    %   degC)
        
    

    
    %   file name format for VAV zones is "zone_XXXXXX.xlsx" where XXXXXX
    %   numeric values indicating the controller ID.
    %   example: zone_431282.xlsx
    
    %   file name format for the AHU serving the VAVs listed in the folder
    %   is "ahu_XXXXXX.xlsx" where XXXXXX numeric values indicating the 
    %   controller ID.
    %   example: ahu_431200.xlsx   
    
    %   file name for the weather data is "weather.csv" 

    %   time series data in each file must have identical start/stop
    %   dates
    %   hourly intervals and a full calendar year (8760 h) are recommended
    
    %   Zone data files contain time series data in the following format
    %       column 1 - time strings (yyyy-mm-dd hh:mm)
    %       column 2 - indoor temperature (degC)
    %       column 3 - vav airflow rate (L/s)
    %       column 4 - vav airflow setpoint (L/s)
    %       column 5 - damper position (%) 
    %       column 6 - perimeter heater valve position (%) (if available)
    
    %   ahu data file contains time series data in the following format
    %       column 1 - time strings (yyyy-mm-dd hh:mm)
    %       column 2 - supply air temperature (degC)
    %       column 3 - return air temperature (degC)
    %       column 4 - heating coil valve (%)
    %       column 5 - cooling coil valve (%)
    %       column 6 - outdoor air damper (%)
    %       column 7 - fan state (%) 
    %       column 8 - supply air pressure (Pa)
    
    %   weather data file contains time series data in the following format
    %       column 1 - time strings (yyyy-mm-dd hh:mm)
    %       column 2 - outdoor air temperature (degC)
    %       column 3 - global horizontal solar irradiance (W/m2)

    %       data required can be retrieved from NASA Power App
    
    %   outputs:
    
    %   energy 
	
        %   energy.existing.ahuHtg: ahu heating energy use in existing state [kWh]
        %   energy.existing.ahuClg: ahu cooling energy use in existing state [kWh]
        %   energy.existing.zoneHtg: zone perimeter heating energy use in existing state [kWh]
			
        %   energy.corrected.ahuHtg: ahu heating energy use in corrected state [kWh]
        %   energy.corrected.ahuClg: ahu cooling energy use in corrected state [kWh]
        %   energy.corrected.zoneHtg: zone perimeter heating energy use in corrected state [kWh]		
		
	    %	inpect figures Existing Energy Use.png and Existing Operation.png
		%   Corrected Energy Use.png and Corrected Operation.png to visualize
		%	and interpret the results.			
		
	% 	znMdlPrmtr
		
		% 	a table of zone model parameters
			% column 1 - Zone name identifier
			% column 2 - C thermal capacitance (kJ/degC)
			% column 3 - U thermal transmittance (kW/degC)
			% column 4 - Qprm perimeter heater capacity (kW)

	%   inpect figures generated for each zone for (C,U, Qprm, disturbances Qmisc)
	%   for supplementary information

    % load data
    currentDirectory = pwd;
    cd(dataDirectory);

    files=dir('*ahu*');
    [num,txt,raw] =xlsread(files(1).name,'data');
    t = datenum(txt(2:end,1)); % matlab date time
    tSa = num(:,1); % supply air temperature (degC)
    tRa = num(:,2); % return air temperature (degC)
    sHc = num(:,3); % heating coil valve (%)
    sCc = num(:,4); % cooling coil valve (%)
    sOa = num(:,5); % mixing box damper (%)
    sFan = num(:,6); % fan status (%)
    sFanOn = sFan > 10; % binary fan operating status
    pSa = num(:,7); % supply air pressure (Pa)

    if seasonalAvailability == 1
        htgAvailSt = day(datetime(datestr(t)),'dayofyear') < 130| day(datetime(datestr(t)),'dayofyear') > 283;
    else
        htgAvailSt = ones(size(t));
    end

    % scheduled operating hours as time-series
    if wkndOp == 1
        occupiedHours = (hour(t) > start & hour(t) < stop); 
    else
        occupiedHours = ((hour(t) > start & hour(t) < stop) & (weekday(t) > 1 & weekday(t) < 7));  
    end

    files=dir('*zone*');
    for k=1:length(files)
       zone(k).data =xlsread(files(k).name,'data');
       tIn(:,k) = zone(k).data(:,1); % indoor temperature (degC)
       qFlo(:,k) = zone(k).data(:,2); % airflow rate (L/s)
       qFloSp(:,k) = zone(k).data(:,3); % airflow rate setpoint (L/s)
       sDmp(:,k) = zone(k).data(:,4); % vav terminal damper (%)
       temp = size(zone(k).data);
       if temp(2) == 5
           sRad(:,k) = zone(k).data(:,5); % perimeter heater state (%)
       else
           sRad(:,k) = NaN(size(sDmp(:,1))); % perimeter heater is not available
       end
    end

    [data,txt,raw] = xlsread('weather.csv');
    tWt = datenum(txt(2:end,1)); % time for the weather file
    tOa = data(:,1); % outdoor temperature (degC)
    qSol = data(:,2); % direct normal solar irradiance (W/m2)

    tOa = interp1(tWt,tOa,t);
    qSol = interp1(tWt,qSol,t);

    cd(currentDirectory)

    % Train ahu models
    ind = sHc == 0 & sCc == 0 & sFan > 0; 
    tMa = tSa - (pSa./840);
    fOa = (tMa(ind) - tRa(ind))./(tOa(ind)-tRa(ind)); % outdoor air fraction
    handle = sOa(ind)/100; % normalize mixing box damper [0 1]
    ft = fittype('c/(1+exp(a*x+b))');
    options = fitoptions(ft);
    options.Lower = [-20  -20 -20];
    options.Upper = [20    20  20];
    [mdlDmp,gofDmp] = fit(handle(fOa < 1 & fOa > 0),fOa(fOa < 1 & fOa > 0),ft,options);

        
    % Train zone models
    u = [];
    u(:,2)  = tOa; % outdoor air temperature
    u(:,3)  = tSa; % supply air temperature
    u(:,4)  = qSol; % solar radiation
    u(:,5)  = sFanOn; % fan on/off status
    u(:,6) = htgAvailSt; % heating plant availability status

    x_znMdl = []; % matrix containing the parameters of the zone models
    disturbances = [];
    for m = 1:length(files)    
        u(:,1) = tIn(:,m);
        u(:,7) = qFlo(:,m);
        u(:,8) = sRad(:,m);
        %u(:,8) = movmean(u(:,8),[1 0]);

        [C,U,qRad,qMisc] = zoneModel(t,u,files(m).name);
        x = [C,U,qRad];
        x_znMdl = [x_znMdl;x];
        disturbances = [disturbances,qMisc];
    end

    % estimate existing heating and cooling energy use

    % 1) totalize total airflow received by zones
        qFloTotal = sum(qFlo,2)/1000;
    % 2) calculate temperature rise across heating coil & cooling coil
        tMa = (mdlDmp(sOa./100).*tOa + (1-mdlDmp(sOa./100)).*tRa);
        dT = tSa - tMa - (pSa./840);
    % 3) compute energy use by AHU heating and cooling coils
    if seasonalAvailability == 1
        qHtgAhu = (dT > 0).* dT.*1.2.*qFloTotal.*((htgAvailSt == 1).*(sHc > 0 & sFanOn));
        qClgAhu = (dT < 0).*-dT.*1.2.*qFloTotal.*((htgAvailSt == 0).*(sCc > 0 & sFanOn));
        for i = 1:length(files)
            qHtgZone(:,i) = x_znMdl(i,3).*sRad(:,i)./100.*htgAvailSt;
        end
        qHtgZone = sum(qHtgZone,2,'omitnan');
    else
        qHtgAhu = (dT > 0).*dT.*1.2.*qFloTotal.*(sHc > 0 & sFanOn);
        qClgAhu = (dT < 0).*-dT.*1.2.*qFloTotal.*(sCc > 0 & sFanOn);
        for i = 1:length(files)
            qHtgZone(:,i) = x_znMdl(i,3).*sRad(:,i)./100;
        end
        qHtgZone = sum(qHtgZone,2,'omitnan');
    end


    fig = figure('units','inch','position',[0,0,4,6]);
        labels = ["AHU heating (kW)","Zone heating (kW)","AHU cooling (kW)"]; 
        h = stackedplot(timetable(datetime(datestr(t)),qHtgAhu,qHtgZone,qClgAhu),"DisplayLabels",labels);
        ax = findobj(h.NodeChildren, 'Type','Axes');
        set([ax.YLabel],'Rotation',90,'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom')
        ax(1).TickDir = 'out';
        ax(2).TickDir = 'out';
        ax(3).TickDir = 'out';
        text(ax(3),5,max(qHtgAhu), strcat(num2str(mean(qHtgAhu)*8760,6),' kWh'))
        text(ax(2),5,max(qHtgZone), strcat(num2str(mean(qHtgZone)*8760,6),' kWh'))
        text(ax(1),5,max(qClgAhu), strcat(num2str(mean(qClgAhu)*8760,6),' kWh'))
        h.LineProperties(1).Color = [0.2 0.0 0.0];
        h.LineProperties(2).Color = [0.2 0.2 0.2];
        h.LineProperties(3).Color = [0.0 0.0 0.2];
    print(fig,'Existing Energy Use.png','-r600','-dpng');

    energy.existing.ahuHtg = qHtgAhu;
    energy.existing.ahuClg = qClgAhu;
    energy.existing.zoneHtg = qHtgZone;

    fig = figure('units','inch','position',[0,0,4,6]);
        labels = ["AHU airflow(L/s)","Outdoor airflow (L/s)","25th percentile","75th percentile","Supply air temp (" + char(176) + "C)"];
        TT = array2timetable([qFloTotal.*1000,(qFloTotal.*1000).*mdlDmp(sOa./100),prctile(tIn,25,2),prctile(tIn,75,2),tSa],...
        'RowTimes', datetime(datestr(t)), ...
        'VariableNames', labels);
        h = stackedplot(TT,{1,2,[3,4],5}); 
        h.AxesProperties(3).LegendLocation = 'southwest';
        ax = findobj(h.NodeChildren, 'Type','Axes');
        ax(2).YLim = [15 30];
        set([ax.YLabel],'Rotation',90,'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom')
        set([ax(2).YLabel],'String',"Zone temp (" + char(176) + "C)")
        h.LineProperties(1).Color = [0.2 0.2 0.2];
        h.LineProperties(2).Color = [0.2 0.2 0.2];
        h.LineProperties(4).Color = [0.2 0.2 0.2];
        ax(1).TickDir = 'out';
        ax(2).TickDir = 'out';
        ax(3).TickDir = 'out';
		ax(4).TickDir = 'out';
    print(fig,'Existing operation.png','-r600','-dpng');

    % estimate corrected heating and cooling energy use

    % initialize variables
        tInSim(1,:) = tIn(1,:); % initialize the indoor temperature
        qFloSim(1,:) = qFlo(1,:); % initialize the VAV discharge air volume
        tSaSim = 18; % initialize the supply air temperature
        zonesWthNoHtr = isnan(sRad(1,:)); % zones with no heater
        qPrmHtgSim(1,:) = zeros(size(qFloSim(1,:))); % initialize perimeter heater loads
        qFloSpMin = min(qFloSp(occupiedHours,:)); % min airflow setpoint
        qFloSpMax = max(qFloSp(occupiedHours,:)); % max airflow setpoint
        qFloSpMin((qFloSpMin<qFloSpMax.*0.1)) = qFloSpMax((qFloSpMin<qFloSpMax.*0.1)).*0.1;
        wakeupRequest = 0; % individual wakeup requests from zones
        wakeup = 0; % wakeup call for the unoccupied state to engage setup/setback modes

    for j = 1:length(t) % j is the hourly timesteps
        for i = 1:length(files) % loop through zones, i zone index
            for k = 1:60 % divide the timestep into 1-min increments
                if occupiedHours(j) == 1 
                    if tInSim(j,i) < htgSp
                        htgClgDmd(j,i) = htgSp - tInSim(j,i); % positive for heating
                    elseif tInSim(j,i) > clgSp
                        htgClgDmd(j,i) = clgSp - tInSim(j,i); % negative for cooling
                    else
                        htgClgDmd(j,i) = 0; % in deadband
                    end
                else
                    if tInSim(j,i) < htgSpUnocc % unoccupied too cold threshold
                        htgClgDmd(j,i) = htgSpUnocc - tInSim(j,i); % positive for heating
                        wakeupRequest = wakeupRequest + 1;
                    elseif tInSim(j,i) > clgSpUnocc % unoccupied too hot threshold
                        htgClgDmd(j,i) = clgSpUnocc - tInSim(j,i); % negative for cooling
                        wakeupRequest = wakeupRequest + 1;
                    else
                        htgClgDmd(j,i) = 0; % in deadband
                    end
                end

            qPrmHtgSim(j,i) = 0;
            % decide perimeter heating use
            if htgClgDmd(j,i) > 0
                if zonesWthNoHtr(i) == 0 & seasonalAvailability == 1 % heating available in zone i
                    if htgAvailSt(j) == 1
                        qPrmHtgSim(j,i) = 100*min(htgClgDmd(j,i),1); % ramp-up heating
                    end
                elseif zonesWthNoHtr(i) == 0
                    qPrmHtgSim(j,i) = 100*min(htgClgDmd(j,i),1);
                end
            end

            % decide discharge VAV airflow rate
            if occupiedHours(j) == 1
                qFloSim(j,i) = qFloSpMin(i); % by default keep it at minimum spt
                if htgClgDmd(j,i) < 0 % if cooling is needed
                    if tSaSim(j,1) < tInSim (j,i) % if AHU SAT is cool enough
                        qFloSim(j,i) = qFloSpMin(i) + (qFloSpMax(i) - qFloSpMin(i))* min(abs(htgClgDmd(j,i)),1); % ramp up airflow rate
                    end
                end
            elseif wakeup == 1
                qFloSim(j,i) = 0;
                if htgClgDmd(j,i) < 0 % if cooling is needed
                    if tSaSim(j,1) < tInSim (j,i) % if AHU SAT is cool enough
                        qFloSim(j,i) = qFloSpMax(i)* min(abs(htgClgDmd(j,i)),1); % ramp up airflow rate
                    end
                end
            else 
                qFloSim(j,i) = 0;
            end
                qVavSim(j,i) = (tSaSim(j,1) - tInSim(j,i))*1.2*qFloSim(j,i)/1000;
                sFanOnSim(j,1) = (occupiedHours(j) == 1 | wakeup == 1);
                qEnvSim(j,i) = (tOa(j,1) - tInSim(j,i))*x_znMdl(i,2);
                
                qNet = disturbances(j,i) + qEnvSim(j,i) + qPrmHtgSim(j,i)*x_znMdl(i,3)/100 + qVavSim(j,i);
                dT = qNet/x_znMdl(i,1);
                tInSim(j,i) = dT*60 + tInSim(j,i);
                qFloSimTemp(k) = qFloSim(j,i);
                qPrmHtgSimTemp(k) = qPrmHtgSim(j,i);
            end

            qFloSim(j,i) = mean(qFloSimTemp);
            qPrmHtgSim(j,i) = mean(qPrmHtgSimTemp);
            tInSim(j+1,i) = tInSim(j,i);

        end

        % estimate return air temperature based on weighted average of air
        % supplied to each zone
        if occupiedHours(j) == 1
            tRaSim(j,1) = sum((qFloSim(j,:).*tInSim(j,:)))./sum(qFloSim(j,:));  
        else
            tRaSim(j,1) = mean(tInSim(j,:));
        end

        totAirRequest(j) = sum(qFloSim(j,:)); % compute total airflow requested at time j

        % define SAT bounds based on outdoor temperature
        lowerBound = max(min(14.5 - 0.21*tOa(j,1),16),12);
        upperBound = max(min(18.3 - 0.28*tOa(j,1),20),14);

        if occupiedHours(j) == 1 | wakeup == 1 
            numZnNeedClg(j) = sum(tInSim(j,:) > clgSp + 1); % compute number of zones that need cooling at time j
            numZnNeedHtg(j) = sum(tInSim(j,:) < htgSp - 1); % compute number of zones that need heating at time j
            if numZnNeedHtg(j) < numZnNeedClg(j) % if more zones need cooling than heating
                tSaSim(j+1,1) = min(max(tSaSim(j,1) - 0.5,lowerBound),upperBound);
            elseif numZnNeedHtg(j) > numZnNeedClg(j) % if more zones need heating than cooling
                tSaSim(j+1,1) = max(min(tSaSim(j,1) + 0.5,upperBound),lowerBound);
            else
                tSaSim(j+1,1) = tSaSim(j,1);
            end
        else
            tSaSim(j+1,1) = tSaSim(j,1);
        end

        wakeup = 0; % wakeup call for the unoccupied state to engage setup/setback modes
        if wakeupRequest > 60
            wakeup = 1;
        end
        wakeupRequest = 0; % reset individual wakeup requests from zones

        fOaMinOa = max(5.562/(1+exp(-4.708*sOaMin/100+6.226)),0.1);
        tMaSimMinOa(j,1) = fOaMinOa*tOa(j,1) + (1-fOaMinOa)*tRaSim(j,1); % mixed air temperature at sOaMin
        if tOa(j,1) < tSaSim(j,1) % economizer plausible range
            if tMaSimMinOa(j,1) < tSaSim(j,1) % mixed air temperature at sOaMin is lower than tSaSim
                htgClgDmdAhu(j,1) = tSaSim(j,1) - tMaSimMinOa(j,1); % temperature rise required across heating coil
                mode(j,1) = 1;
                fOaSim(j,1) = fOaMinOa; % minimum outdoor air fraction
                sOaSim(j,1) = sOaMin;
            else
                fOaSim(j,1) = min((tSaSim(j,1) - tRaSim(j,1))/(tOa(j,1) - tRaSim(j,1)),1); % outdoor air fraction
                z = 1;
                handle = [];
                for s = sOaMin:5:100
                    handle(z,1) = 5.562/(1+exp(-4.708*s/100+6.226));
                    z = z + 1;
                end
                [ind,ind] = min(abs(handle - fOaSim(j,1)));
                s = sOaMin:5:100;
                sOaSim(j,1) = s(ind);
                htgClgDmdAhu(j,1) = tSaSim(j,1) - (fOaSim(j,1)*tOa(j,1) + (1-fOaSim(j,1))*tRaSim(j,1)); % economizer state (no temperature drop required across) - if damper broken, may need cooling
                mode(j,1) = 2;
            end
        elseif tRaSim(j,1) > tOa(j,1) & tSaSim(j,1) < tOa(j,1) & tOa(j,1) < 20
            htgClgDmdAhu(j,1) = tSaSim(j,1) - tOa(j,1); % temperature drop required across cooling coil (economizer with cooling mode (100% outdoor air))
            fOaSim(j,1) = 1;
            sOaSim(j,1) = 100;
            mode(j,1) = 3;
        elseif tMaSimMinOa(j,1) > tSaSim(j,1)
            htgClgDmdAhu(j,1) = tSaSim(j,1) - tMaSimMinOa(j,1); % temperature drop required across cooling coil (cooling mode)
            fOaSim(j,1) = fOaMinOa; % minimum outdoor air fraction
            sOaSim(j,1) = sOaMin;
            mode(j,1) = 4;
        else
            htgClgDmdAhu(j,1) = 0;
            fOaSim(j,1) = fOaMinOa; % minimum outdoor air fraction
            sOaSim(j,1) = sOaMin;
            mode(j,1) = 5;
        end

    end

    % estimate corrected heating and cooling energy use

    % 1) totalize total airflow received by zones
        qFloSimTotal = sum(qFloSim,2)/1000;
    % 2) compute energy use by AHU heating and cooling coils
        if seasonalAvailability == 1
            qHtgAhuSim = (htgClgDmdAhu > 0).* htgClgDmdAhu.*1.2.*qFloSimTotal.*(htgAvailSt);
            qClgAhuSim = (htgClgDmdAhu < 0).*-htgClgDmdAhu.*1.2.*qFloSimTotal.*(htgAvailSt == 0);
            for i = 1:length(files)
                qHtgZoneSim(:,i) = x_znMdl(i,3).*qPrmHtgSim(:,i)./100.*htgAvailSt;
            end
            qHtgZoneSim = sum(qHtgZoneSim,2,'omitnan');
        else
            qHtgAhuSim = (htgClgDmdAhu > 0).* htgClgDmdAhu.*1.2.*qFloSimTotal;
            qClgAhuSim = (htgClgDmdAhu < 0).*-htgClgDmdAhu.*1.2.*qFloSimTotal;
            for i = 1:length(files)
                qHtgZoneSim(:,i) = x_znMdl(i,3).*qPrmHtgSim(:,i)./100;
            end
            qHtgZoneSim = sum(qHtgZoneSim,2,'omitnan');
        end


    fig = figure('units','inch','position',[0,0,4,6]);
        labels = ["AHU heating (kW)","Zone heating (kW)","AHU cooling (kW)"];
        h = stackedplot(timetable(datetime(datestr(t)),qHtgAhuSim,qHtgZoneSim,qClgAhuSim),"DisplayLabels",labels);
        ax = findobj(h.NodeChildren, 'Type','Axes');
        set([ax.YLabel],'Rotation',90,'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom')
        ax(1).TickDir = 'out';
        ax(2).TickDir = 'out';
        ax(3).TickDir = 'out';
        text(ax(3),5,max(qHtgAhuSim), strcat(num2str(mean(qHtgAhuSim)*8760,6),' kWh'))
        text(ax(2),5,max(qHtgZoneSim), strcat(num2str(mean(qHtgZoneSim)*8760,6),' kWh'))
        text(ax(1),5,max(qClgAhuSim), strcat(num2str(mean(qClgAhuSim)*8760,6),' kWh'))
        h.LineProperties(1).Color = [0.2 0.0 0.0];
        h.LineProperties(2).Color = [0.2 0.2 0.2];
        h.LineProperties(3).Color = [0.0 0.0 0.2];
    print(fig,'Corrected Energy Use.png','-r600','-dpng');

    energy.corrected.ahuHtg = qHtgAhuSim;
    energy.corrected.ahuClg = qClgAhuSim;
    energy.corrected.zoneHtg = qHtgZoneSim;

    fig = figure('units','inch','position',[0,0,4,6]);
        labels = ["AHU airflow(L/s)","Outdoor airflow (L/s)","25th percentile","75th percentile","Supply air temp (" + char(176) + "C)"];
        TT = array2timetable([qFloSimTotal(1:end,:).*1000,(qFloSimTotal(1:end,:).*1000).*fOaSim(1:end,:),prctile(tInSim(1:end-1,:),25,2),prctile(tInSim(1:end-1,:),75,2),tSaSim(1:end-1,:)],...
        'RowTimes', datetime(datestr(t)), ...
        'VariableNames', labels);
        h = stackedplot(TT,{1,2,[3,4],5}); 
        h.AxesProperties(3).LegendLocation = 'southwest';
        ax = findobj(h.NodeChildren, 'Type','Axes');
        ax(2).YLim = [15 30];
        set([ax.YLabel],'Rotation',90,'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom')
        set([ax(2).YLabel],'String',"Zone temp (" + char(176) + "C)")
        h.LineProperties(1).Color = [0.2 0.2 0.2];
        h.LineProperties(2).Color = [0.2 0.2 0.2];
        h.LineProperties(4).Color = [0.2 0.2 0.2];
        ax(1).TickDir = 'out';
        ax(2).TickDir = 'out';
        ax(3).TickDir = 'out';
		ax(4).TickDir = 'out';
    print(fig,'Corrected operation.png','-r600','-dpng');
    
    % generate a plot defining the behaviour
    
    fig = figure('units','inch','position',[0,0,2,2]);
        ind = sHc == 0 & sCc == 0 & sFan > 0; 
        handle = sOa(ind)./100;
        scatter(handle(fOa < 1 & fOa > 0),fOa(fOa < 1 & fOa > 0),12,[0 0 1],...
            'filled','o','MarkerFaceAlpha',.1,'MarkerEdgeAlpha',.1)
        hold on
        plot((0.1:0.01:1)',min(mdlDmp(0.1:0.01:1),1),'r','LineWidth',2)
        % plot((0.3:0.01:1)',max(5.562./(1+exp(-4.708.*(((30:1:100)')./100)+6.226)),0.1),'k','LineWidth',2)
        
        upperLimit = [0;23;46;68;85;93;97;99;100;100]./100;
        lowerLimit = [0;1;3;5;9;15;25;39;64;100]./100;
        x = linspace(0,1,length(upperLimit));
        
        limits = patch([x' fliplr(x')], [lowerLimit fliplr(upperLimit)], 'g');
        limits.FaceAlpha = 0.1;
        limits.EdgeColor = 'w';
        limits.EdgeAlpha = 0;
        
        ylabel('Outdoor air fraction')
        xlabel('Mixing box damper signal')
        xlim([0 1])
        xticks(0:0.2:1)
        ylim([0 1])
        yticks(0:0.2:1)
        set(gca,'TickDir','out');
        box off
    print(fig,'outdoor air fraction to damper.png','-r600','-dpng');
    
    fig = figure('units','inch','position',[0,0,6,3]);
        subplot(1,2,1)
        ind = sFan > 0; 
        scatter(tOa(ind),mdlDmp(sOa(ind)./100),12,[0 0 1],...
        'filled','o','MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.5)
        xlabel('Outdoor temperature (^{o}C)')
        ylabel('Outdoor air fraction')
        ylim([0 1])
        xlim([-30 30])
        xticks(-30:10:30)
        subplot(1,2,2)
        scatter(tOa,fOaSim,12,[0 0 1],...
        'filled','o','MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.5)
        ylim([0 1])
        xlim([-30 30])
        xticks(-30:10:30)
        xlabel('Outdoor temperature (^{o}C)')
        ylabel('Outdoor air fraction')
    print(fig,'state of operation.png','-r600','-dpng');
    
   
        function [C,U,qRad,qMisc] = zoneModel(t,u,filename)

        % C Thermal mass (kJ/degC) (Scalar)
        % U Thermal transmittance (kW/degC) (Scalar)
        % qRad Perimeter heater capacity (kW) (Scalar)
        % qMisc Random heat sinks/sources (kW) (Array) 

        % u(:,1) indoor temperature (degC)
        % u(:,2) outdoor temperature (degC)
        % u(:,3) supply temperature (degC)
        % u(:,4) direct normal solar radiation (W/m2)
        % u(:,5) fan status (1 or 0)
        % u(:,6) heating availability status (1 or 0)
        % u(:,7) discharge airflow rate (L/s)
        % u(:,8) perimeter heater valve (%)

        % 1) Identify the thermal capacitance (kJ/K) at cooling switch on
        % instances
        dt = median(diff(t))*3600*24; % timestep size by seconds
        C = [];
        qVav = ((u(:,3) - u(:,1))*1.2).*(u(:,7)/1000).*u(:,5);
        clgSwitchOnInd = find(diff(qVav,1) < -0.2 & diff(u(:,5)) > 0 & u(1:end-1,2) > 10); % find instances with a change
        for i = 1:length(clgSwitchOnInd)
            slope_1 = (u(clgSwitchOnInd(i),1)   - u(clgSwitchOnInd(i)-1,1))/dt;
            slope_2 = (u(clgSwitchOnInd(i)+1,1) - u(clgSwitchOnInd(i),1))/dt;
            C(i,1) = (qVav(clgSwitchOnInd(i)+1) - qVav(clgSwitchOnInd(i)))/(slope_2 - slope_1);
        end
        C = C(C > 10^2 & C < 10^4);

        fig = figure('units','inch','position',[0,0,2,2]);
            h = boxplot(C,'Notch','on','Labels',{''},'Colors',[0 0 0]);
            set(gca,'xtick')
            set(h,'LineWidth', 1);
            ylabel({'Thermal mass \itC\rm (kJ/kg)'})
            set(gca,'TickDir','out');
            set(gca,'xtick',[])
            ylim([0 10000])
            box off
        print(fig,strcat('C_',filename,'.png'),'-r600','-dpng');

        C = prctile(C,50);

        if isnan(C) | C == 0
            C = 1000;
        end

        % 2) Identify the thermal transmittance (kW/degC) during free floating
        % in winter
        if sum(isnan(u(:,8))) == 0
            wntrFreeFloatInd = u(1:end-1,2) < 0 & u(1:end-1,4) < 20 & u(1:end-1,5) == 0 & u(1:end-1,6).*u(1:end-1,8) == 0 & diff(u(:,1)) < 0;
            start1 = strfind(wntrFreeFloatInd',[0 1 1 1 1 1]);

            if length(start1) > 5
                for i = 1:length(start1)
                    dQ(i,1) = (u(start1(i)+5,1) - u(start1(i)+1,1))*C/(dt*4); % heat released by the thermal mass over 4 timesteps
                    dTemp(i,1) = mean((u(start1(i)+1:start1(i)+5,2) - (u(start1(i)+1:start1(i)+5,1)))); % temperature difference between indoors and outdoors
                end

                mdl = fitlm(dTemp,dQ,'RobustOpts','on');
                U = mdl.Coefficients.Estimate(2);

                    fig = figure('units','inch','position',[0,0,2,2]);
                            h = scatter(dTemp,dQ,2,[0.3 0.3 0.3]);
                            hold on
                            plot((-50:10:-20)',predict(mdl,((-50:10:-20)')),'k','LineWidth',2)
                            text(-45,max(dQ),strcat('R^{2} = ',num2str(mdl.Rsquared.Adjusted,2)));
                            xlabel({'T_{oa} - T_{in} (^{o}C)'})
                            ylabel({'Heat released by mass (kW)'})
                            set(gca,'TickDir','out');
                            xlim([-50 -20])
                            box off
                    print(fig,strcat('U_',filename,'.png'),'-r600','-dpng');

                if U < 0
                    U = 0;
                elseif U > 0.5
                    U = 0.5;
                end
            else
                U = 0;
            end
        else
            U = 0;
        end

        % 3) Compute radiator capacity at heating switch on instances

        if sum(isnan(u(:,8))) == 0
        sRadTemp = u(:,6).*u(:,8) > 10;

        start1 = intersect(strfind(sRadTemp',[0 0 0 1 1]),find(u(:,2) < 0)); % find instances when radiator was off then turned on for at least two hours in Winter

            if length(start1) > 5
                for i = 1:length(start1)
                    if sum(isnan(u(:,8))) == 0
                        slope_1 = (u(start1(i)+2,1) - u(start1(i)+1,1))*C/dt - qVav(start1(i)+2);
                        slope_2 = (u(start1(i)+4,1) - u(start1(i)+2,1))*C/(dt*2) - mean(qVav(start1(i)+3:start1(i)+4));
                        qRad(i,1) = slope_2 - slope_1;
                        sRadFiltered(i,1) = mean(u((start1(i)+3:start1(i)+4),8));
                        a(i,1) = qRad(i,1)/sRadFiltered(i,1).*100;
                    end       
                end
                    a = a(a > 0 & a < 5);
                    fig = figure('units','inch','position',[0,0,2,2]);
                        h = boxplot(a,'Notch','on','Labels',{''},'Colors',[0 0 0]);
                        set(gca,'xtick')
                        set(h,'LineWidth', 1);
                        ylabel({'Q_{prm} (kW)'})
                        set(gca,'TickDir','out');
                        set(gca,'xtick',[])
                        ylim([0 5])
                        box off
                    print(fig,strcat('Qrad_',filename,'.png'),'-r600','-dpng');              
                    qRad = median(a);

            else
                qRad = 0.5;
            end
        else
            qRad = 0;
        end

        % 4) Estimate misc heat sources/sinks
        if sum(isnan(u(:,8))) == 0
        qMisc = diff(u(:,1)).*(C/dt) + ... % heat storage/release rate (kW)
              - qVav(2:end,1) + ... % heat added/removed by the VAV (kW) 
              - u(2:end,6).*u(2:end,8).*qRad/100 + ... % heat added by the radiator (kW)
              -(u(2:end,2) - u(2:end,1)).*U; % heat exchange from the envelope (kW)
        else
        qMisc = diff(u(:,1)).*(C/dt) + ... % heat storage/release rate (kW)
              - qVav(2:end,1) + ... % heat added/removed by the VAV (kW) 
              -(u(2:end,2) - u(2:end,1)).*U; % heat exchange from the envelope (kW)    
        end

        qMisc = [qMisc;qMisc(end,1)];  
        qMisc = filloutliers(qMisc,'linear');
        % qMisc(qMisc < 0) = 0;

        timeOfDay = hour(t);
        dayOfWeek = weekday(t);

        for i = 1:7
            for j = 0:23
                ind = timeOfDay == j & dayOfWeek == i;
                qMisc25th(j+1,i) = prctile(qMisc(ind,1),25);
                qMisc75th(j+1,i) = prctile(qMisc(ind,1),75);
            end
        end

        fig = figure('units','inch','position',[0,0,5,2]);
            subplot(1,2,1)
            plot((0:23)',mean(qMisc25th(:,2:6),2),'r')
            hold on
            plot((0:23)',mean(qMisc75th(:,2:6),2),'k')
            ylim([0 max(max(qMisc75th))])
            xlabel({'Time of day'})
            ylabel({'Q_{misc} (kW)'})
            set(gca,'TickDir','out');
            xlim([0 24])
            xticks([0:3:24])
            legend('25^{th} percentile','75^{th} percentile','Location','northoutside')
            legend('boxoff')
            text(1,max(max(qMisc75th)),'(a)')
            box off

            subplot(1,2,2)
            plot((0:23)',mean(qMisc25th(:,[1,7]),2),'r')
            hold on
            plot((0:23)',mean(qMisc75th(:,[1,7]),2),'k')
            ylim([0 max(max(qMisc75th))])
            xlabel({'Time of day'})
            ylabel({'Q_{misc} (kW)'})
            set(gca,'TickDir','out');
            xlim([0 24])
            xticks([0:3:24])
            legend('25^{th} percentile','75^{th} percentile','Location','northoutside')
            legend('boxoff')
            box off
            text(1,max(max(qMisc75th)),'(b)')
        print(fig,strcat('qMisc',filename,'.png'),'-r600','-dpng');
        close all   
    end
	
	znMdlPrmtr = table(string({files.name})', x_znMdl);
    
end
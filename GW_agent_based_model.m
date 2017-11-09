%% Randomization seeds
load('seed1.mat')
load('seed2.mat')
load('seed3.mat')
rng(seed1)
%rng('shuffle')

%% Simulation paramters
% Number of simulation years
simYears = 10;
% Total number of dogs
nbDogs = 30000;
% Maximum number of worms in water
maxGW = nbDogs/3;
% Parameter for the infection function
% Recommended range: 0.0001 - 0.001
para = 0.0005;
% Number of lizards in simulation
nbLiz = 30000;
% Initial average infection rate
infProb0 = 0.001;
infProbMin = 0.0001;
% Number of iterations (of simulation)
nbIt = 1; %10
% Unit infection rate in water per worm
maxRate = 0.2;
rainySeason = 0.5;
% Exit day of each worm, distributed uniformaly between 10 and 14 months
daysExLb = 12 * 30;
daysExUb = 14 * 30;
% Plus a 10-15 days until development into L3 larvae
daysInfLb = 10;
daysInfUb = 15;
GWWaterLive = 30;
% Number of simulation days
simDays = 360 * simYears + daysExUb + daysInfUb;
% Number of days in rainy season June - Oct each year
rainyStart = 5 * 30 + 1;
rainyEnd = 10 * 30;
initialRainyProb = 0/3;
% Number of worms per infection is distributed exponentially with rate 1
avgNbGW = 1;
% Probability of worm contacting water after exude from dog
contactProb = 1;
% Initializing variables
avgDogInf = 0;
avgNewInf = 0;
avgNewDogInf = 0;
avgGWEx = 0;
avgInfRate = 0;
% Type of infection
waterInfect = 1;
lizInfect = 0;
daysLizLb = 5 * 30;
daysLizUb = 36 * 30;

%% Intervention
% Tethering dogs that exude worms
tt = 1;
% Probability of discovering a dog exuding worm
ttProb = [0, 0.4, 2/3, 2/3, 0.76]; %first few months have lower probability
targetTtProb = 0.95;
ttProb = [ttProb, ones(1, simYears - length(ttProb) + 1) * targetTtProb];
% ttProb = ones(1, simYears + 1) * 0.95; %constant probability for tethering
if length(ttProb) ~= simYears + 1
    error('Lenght of ttProb should equal to simYears + 1');
end
% Length of days for each thethering
ttLength = 30;

% Water treatment
ABATE = 1;
ABATEProb = [0.25, 0.25, 0.25, 0.25, 0.25];
targetABATEProb = 0.95;
ABATEProb = [ABATEProb, ones(1, simYears - length(ABATEProb) + 1) * targetABATEProb];
ABATELength = 30;

%% Simulation runs
for iIt = 1:nbIt
    %% Initialization of variables
    nbGWWater = 0; %number of worms in the water
    
    nbGWEx = zeros(simDays, 1); %number of worms exuding
    nbDogInf = zeros(simDays, 1); %total dog infections
    nbNewDogInf = zeros(simDays, 1); %new dog infections
    nbNewInf = zeros(simDays, 1); %new HUMAN infections
    rateWaterByDay = zeros(simDays, 1); %??
    dogInf0 = binornd(nbDogs, infProb0); %number of initial dog infections
    dogsInf = round(nbDogs * rand(dogInf0, 1), 0); %??
    
    nbLizInf = 0; %number of lizards infected
    rateLizByDay = zeros(simDays, 1); %rate of lizard (infections?) per day
    lizInf = zeros(nbLiz, 1); %number of lizards infected
    lizRecDay = zeros(nbLiz, 1); %??
    
    % Array indicating the index of dogs each worm lives in
    GWInDog = [];
    % Array indicating which day worms die in water
    GWDieDay = [];
    % Binary array indicating whether a dog is thethered
    ttDogs = zeros(nbDogs, 1); %1=tethered
    % Indicate when thethering ends
    ttDogsEnd = ttDogs;
    ABATEEnd = 0;
    for iDogs = 1:length(dogsInf)
        %randomized number of worms per dog
        GWInDog = [GWInDog; repmat(dogsInf(iDogs), ...
            floor(exprnd(avgNbGW)) + 1, 1)];
    end
    % Exude day for each worm
    GWExudeDay = round((daysExUb - daysExLb) * ...
        rand(length(GWInDog), 1) + daysExLb + daysInfLb + ...
        (daysInfUb - daysInfLb) * rand(length(GWInDog), 1), 0);
    for iWorm = 1:length(GWExudeDay)
        if rand() <= 1 - initialRainyProb
            if rand() <= (rainyStart - 1) / (360 - (rainyEnd - rainyStart))
                GWExudeDay(iWorm) = GWExudeDay(iWorm) + ...
                    round(rand() * (rainyStart - 1), 0);
            else
                GWExudeDay(iWorm) = GWExudeDay(iWorm) + ...
                    round(rand() * (360 - rainyEnd - 1) + rainyEnd + 1, 0);
            end
        else
            GWExudeDay(iWorm) = GWExudeDay(iWorm) + ...
                round(rand() * (rainyEnd - rainyStart) + rainyStart, 0);
        end
    end
    
    %% Simulate each day
    for iDays = 1:simDays
        currentYear = ceil(iDays/360);
        if currentYear > length(ttProb)
            currentYear = length(ttProb);
        end
        if mod(iDays, 100) == 0
            fprintf("Day %4d... ",iDays);
        end
        if mod(iDays, 1000) == 0
            fprintf("\n");
        end
        % Undo thethering in dogs
        ttDogs(ttDogsEnd < iDays) = 0;
        % Worm die in water on iDays
        if sum(GWDieDay == iDays) > 0
            nbGWWater = nbGWWater - sum(GWDieDay == iDays);
            GWDieDay(GWDieDay == iDays) = [];
        end
        if length(GWDieDay) ~= nbGWWater
            error('wat1');
        end
        
        %         if sum(lizRecDay == iDays) > 0
        %             ind = find(lizRecDay == iDays);
        %             nbLizInf = nbLizInf - length(ind);
        %             lizRecDay(ind) = 0;
        %         end
        
        if sum(GWExudeDay == iDays) > 0
            
            ind = find(GWExudeDay == iDays);
            nbGWEx(iDays) = length(ind);
            
            if tt == 1 %tethering intervention
                % index of non thethered dogs
                for iDogs = 1:length(ind)
                    if rand() <= ttProb(currentYear)
                        ttDogs(GWInDog(ind(iDogs))) = 1; %tethered
                        ttDogsEnd(GWInDog(ind(iDogs))) = iDays + ttLength;
                    end
                end
                nonTtInd = find(ttDogs == 0); %indices for nontethered dogs
                ind1 = ind(ismember(GWInDog(ind), nonTtInd)); %indices of at-risk dogs
            else
                ind1 = ind; %all dogs with exuding worms are at-risk
            end
            % determine how many new worms are released into water
            if isempty(ind1)
                nbNewGWWater = 0;
            elseif ABATE == 1 && ABATEEnd >= iDays
                nbNewGWWater = binornd(length(ind1), contactProb * (1 - ABATEProb(currentYear)));
            else
                nbNewGWWater = binornd(length(ind1), contactProb);           
            end     
            nbGWWater = nbGWWater + nbNewGWWater;
            % delete worms from dogs
            GWExudeDay(ind) = [];
            GWInDog(ind) = [];
            % Set when worms die in water
            GWDieDay = [GWDieDay; ...
                (iDays + GWWaterLive) * ones(nbNewGWWater, 1)];   
            if nbGWWater ~= length(GWDieDay)
                error('wat2');
            end
            if ABATE == 1 && nbNewGWWater > 0 && ABATEEnd < iDays
                ABATEEnd = iDays + ABATELength;
                nbGWKilled = binornd(nbGWWater, ABATEProb(currentYear));
                if nbGWKilled >= nbGWWater
                    nbGWWater = 0;
                    GWDieDay = [];
                    if nbGWWater ~= length(GWDieDay)
                        error('wat3');
                    end
                elseif nbGWKilled > 0
%                     tempGWDieDay = GWDieDay;
%                     for iKill = 1:nbGWKilled
%                         if rand() <= ABATEProb(currentYear)
%                             tempGWDieDay(iKill) = 0;
%                             nbGWWater = nbGWWater - 1;
%                         end
%                     end
%                     GWDieDay = GWDieDay(tempGWDieDay > 0);
                    nbGWWater = nbGWWater - nbGWKilled;
                    for iKilled = 1:nbGWKilled
                        ind = round(rand() * (length(GWDieDay) - 1), 0) + 1;
                        GWDieDay(ind) = [];
                    end
                    if nbGWWater ~= length(GWDieDay)
                        error('wat4');
                    end      
                end
            end         
        end
        if nbGWWater ~= length(GWDieDay)
            error('Inconsistency on number of worms in water');
        end
        % Decay of worms in dry season
        if mod(iDays, 360) <= rainyEnd && mod(iDays, 360) >= rainyStart
            rainyFactor = rainySeason;
        else
            rainyFactor = 1;
        end
        
        if nbGWWater == 0
            infRateWater = 0;
        else %infection function with 4 parameters
            infRateWater = ((sigmf(min(nbGWWater, maxGW), [para maxGW/2]) - ...
                sigmf(0, [para maxGW/2]) + infProbMin) / ...
                (sigmf(maxGW, [para maxGW/2]) - sigmf(0, [para maxGW/2])) ...
                * (maxRate - 2 * infProbMin) + infProbMin) * rainyFactor;
        end
        % nbGWWater = 1:maxGW;
        % infRateWater = (sigmf(min(nbGWWater, maxGW), [para maxGW/2]) - ...
        %    sigmf(0, [para maxGW/2]))/ ...
        %    (sigmf(maxGW, [para maxGW/2]) - sigmf(0, [para maxGW/2])) ...
        %    * maxRate
        %
        lizInfInd = round(nbLiz * rand(binornd(nbLiz, infRateWater), 1), 0);
        if lizInfInd > 0
            lizInf(lizInfInd) = 1;
            nbLizInf = sum(lizInf);
            lizRecDay(lizInfInd) = round((daysLizUb - daysLizLb) * ...
                rand(length(lizInfInd), 1) + daysLizLb);
        end
        
        rateWaterByDay(iDays) = infRateWater;
        
        %if min(GWExudeDay) <= iDays
        %    break
        %end
        % Calculate number of dogs infected each day
        % distributed exponentially with rate equals to number of total dogs times
        % infection rate
        newGWInDog = [];
        if waterInfect == 1
            dogsInf = round((nbDogs - 1) * ...
                rand(binornd(nbDogs, infRateWater), 1), 0) + 1;
            dogsInf = dogsInf(ttDogs(dogsInf) == 0);
            nbNewInf(iDays) = length(dogsInf);
            newDogs = dogsInf(~ismember(dogsInf, unique(GWInDog)));
            % Find out new worms infected
            for iDog = 1:length(dogsInf)
                newGWInDog = [newGWInDog; repmat(dogsInf(iDog), ...
                    floor(exprnd(avgNbGW)) + 1, 1)];
            end
        end
        if lizInfect == 1
            lizEat = round(nbLiz * rand(round(binornd(nbLiz, maxRate)), 1), 0);
            lizEat(lizEat == 0) = 1;
            dogInfIndLiz = round(nbDogs * rand(sum(lizInf(lizEat)), 1), 0);
            lizInf(lizEat) = 0;
            nbNewInf(iDays) = length(dogInfIndLiz) + nbNewInf(iDays);
            newDogs = dogInfIndLiz(~ismember(dogInfIndLiz, unique(GWInDog)));
            for iDog = 1:length(dogInfIndLiz)
                newGWInDog = [newGWInDog; repmat(dogInfIndLiz(iDog), ...
                    floor(exprnd(avgNbGW)) + 1, 1)];
            end
        end
        
        nbNewDogInf(iDays) = length(newDogs);
        % Update worms and GWExudeDay based on new infections
        GWInDog = [GWInDog; newGWInDog];
        GWExudeDay = [GWExudeDay; ...
            round((daysExUb - daysExLb) * rand(length(newGWInDog), 1) + ...
            daysExLb + daysInfLb + iDays + ...
            (daysInfUb - daysInfLb) * rand(length(newGWInDog), 1), 0)];
        % Update the number of infected dogs with worms
        nbDogInf(iDays) = length(unique(GWInDog));
        %if length(wormInDog) >= 50000
        %    break
        %end
        %nbGWWater
    end
    
    nbDogInfYear = reshape(nbDogInf(1:360 * simYears), 360, simYears);
    nbNewInfYear = reshape(nbNewInf(1:360 * simYears), 360, simYears);
    nbNewDogInfYear = reshape(nbNewDogInf(1:360 * simYears), 360, simYears);
    nbWormExYear = reshape(nbGWEx(1:360 * simYears), 360, simYears);
    infRateYear = reshape(rateWaterByDay(1 : 360 * simYears), 360, simYears);
    nbDogInfMonth = zeros(12, simYears); %***
    nbNewInfMonth = zeros(12, simYears);
    nbNewDogInfMonth = zeros(12, simYears);
    nbWormExMonth = zeros(12, simYears);
    infRateMonth = zeros(12, simYears);
    
    for iMonth = 1:12
        nbDogInfMonth(iMonth, :) = max(nbDogInfYear(30 * (iMonth - 1) + 1: iMonth * 30, :), [], 1);
        nbNewInfMonth(iMonth, :) = sum(nbNewInfYear(30 * (iMonth - 1) + 1: iMonth * 30, :), 1);
        nbNewDogInfMonth(iMonth, :) = sum(nbNewDogInfYear(30 * (iMonth - 1) + 1: iMonth * 30, :), 1);
        nbWormExMonth(iMonth, :) = sum(nbWormExYear(30 * (iMonth - 1) + 1: iMonth * 30, :), 1);
        infRateMonth(iMonth, :) = mean(infRateYear(30 * (iMonth - 1) + 1: iMonth * 30, :), 1);
    end
    avgDogInf = (avgDogInf * (iIt - 1) + nbDogInfMonth)/iIt; %***
    avgNewInf = (avgNewInf * (iIt - 1) + nbNewInfMonth)/iIt;
    avgNewDogInf = (avgNewDogInf * (iIt - 1) + nbNewDogInfMonth)/iIt;
    avgGWEx = (avgGWEx * (iIt - 1) + nbWormExMonth)/iIt; %***
    avgInfRate = (avgInfRate * (iIt - 1) + infRateMonth)/iIt;
end

%% Plotting
figure
plot(1:12,avgGWEx(:,1))
hold on
for year=2:4 %simYears
    plot(1:12,avgGWEx(:,year))
end

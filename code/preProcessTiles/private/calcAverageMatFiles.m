function calcAverageMatFiles(imStack,tileIndex,thisDirName,illumChans,lowValue)

    % function calcAverageMatFiles(imStack,tileIndex,thisDirName,illumChans,lowValue)
    %
    % Purpose
    % Calculates average tile for a section directory and saves as a structure in a .mat file.
    % One file per optical section. This function is called by preProcessTiles
    %
    % Inputs
    % imStack - image stack
    % tileIndex - matrix defining the the locations of the tiles in the array
    % thisDirName - 
    % illumChans - 
    % lowValue - an array the same size as imStack. 
    %
    % In both cases, rows are channels and columns are optical sections
    %
    %

    fprintf('Building and saving average images in %s\n',thisDirName)

    %We calculate all the possible combinations of channels and layers. 
    %This makes it possible to have a parfor loop if needed. 
    chans = repmat(illumChans(:),1,size(imStack,2));
    layers = repmat(1:size(imStack,2), length(illumChans), 1);
    combinations = [chans(:),layers(:)];
    correctionType = 'bruteAverageTrimmean'; %The name of this correction type

    for ii=1:size(combinations,1) 
        thisChan=combinations(ii,1);
        thisLayer=combinations(ii,2);
        fprintf('  Doing channel %d, layer %d/%d\n', thisChan, thisLayer,max(combinations(:,2)) )
        if thisChan==0, continue, end %Necessary?

        aveFname = fullfile(thisDirName,'averages',sprintf('%d%s%02d_%s.mat',thisChan,filesep,thisLayer,correctionType));

        %Create new directories as needed
        aveDirName = fullfile(thisDirName,'averages',sprintf('%d', thisChan));

        if ~exist(aveDirName,'dir')
            mkdir(aveDirName)
        end

        if exist(aveFname), delete(aveFname), end
        imSize = size(imStack{thisChan,thisLayer},1);

        thisStack = imStack{thisChan,thisLayer};

        if isempty(thisStack)
            %Catch the unusual situation of missing data for a particular layer (optical plane) or channel
            fprintf('%s finds that chan %d layer %d is empty. **SKIPPING**\n',mfilename,thisChan,thisLayer)
            continue
        end

        % We will get rid of *really* dim tiles by removing tiles with a mean lower then lowValue
        % The following ensures that we choose reasonable numbers based on the amp offsets
        mu = squeeze(mean(mean(thisStack)));
        lowVals = find(mu<lowValue(thisChan,thisLayer));

        %Fail gracefully if tile index is not complete
        if isempty(tileIndex{thisChan,thisLayer})
            fprintf(' **** WARNING **** Encountered missing data in %s chan: %d layer: %d. SKIPPING\n', ...
                thisDirName, thisChan, thisLayer)
            continue
        end
        %If we have got rid of over 90% of the tiles then don't calculate an average. Likely something is wrong
        row=tileIndex{thisChan,thisLayer}(:,5);

        propRemoved=(length(lowVals)/length(row));
        if propRemoved > 0.85
            fprintf('%s: removed %d%% of tiles from illumination correction. SKIPPING.\n', ...
                mfilename, round(propRemoved*100))
            continue
        end

        row(lowVals)=[];
        thisStack(:,:,lowVals)=[];
        if size(thisStack,3)<2
            fprintf('** WARNING stack size for generating average images is %d. SKIPPING THIS SECTION\n',size(thisStack,3))
            continue
        end

        %Calculate trimmed mean for the even rows (TODO: sub-function? this is duplicate code)
        defaultTrim=1;
        f=find(~mod(row,2));
        if isempty(f)
            error('Operation f=find(~mod(row,2)) has returned empty. Something is very wrong!')
        end

        trimQuantity=round((2/length(f))*100); %Defines the degree of trimming 
        if trimQuantity>=100 | trimQuantity<=0 | isnan(trimQuantity)
            fprintf('ERROR: even data trimQuantity is %d but should be between 0 and 100. Setting to %d\n', ...
                trimQuantity,defaultTrim);
            trimQuantity=defaultTrim;
        end

        evenData = trimmean(thisStack(:,:,f),trimQuantity ,3);
        evenN = length(f);

        %Calculate trimmed mean for the odd rows
        f=find(mod(row,2));
        if isempty(f)
            error('Operation f=find(~mod(row,2)) has returned empty. Something is very wrong!')
        end
        trimQuantity=(2/length(f))*100; %Defines the degree of trimming 
        if trimQuantity>100 | trimQuantity<0
            fprintf('ERROR: odd data trimQuantity is %d but should be between 0 and 100. Setting to %d\n', ...
                trimQuantity,defaultTrim);
            trimQuantity=defaultTrim;
        end
        oddData = trimmean(thisStack(:,:,f),trimQuantity,3);
        oddN = length(f);


        % Build average tile structure and save it
        avData.evenRows = single(evenData);
        avData.oddRows = single(oddData);
        avData.pooledRows = [];
        avData.evenN = evenN;
        avData.oddN = oddN;
        avData.poolN = [];
        avData.correctionType = correctionType;
        avData.channel = thisChan;
        avData.layer = thisLayer;
        avData.details.trimQuantity = trimQuantity;
        avData.details.lowVals = lowVals;

        save(aveFname,'avData')

    end

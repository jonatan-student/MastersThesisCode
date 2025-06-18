function results = load_data(basePath, Folders, concTest)
    filterConc = false;
    if nargin >= 3 && ~isempty(concTest)
        filterConc = true;
        % ensure .mat extension
        if ~endsWith(concTest, '.mat')
            concTest = [concTest, '.mat'];
        end
    end

    results = struct();  % initialize output

    for k = 1:numel(Folders)
        fld     = Folders{k};
        fldPath = fullfile(basePath, fld);

        % get all .mat filenames
        matInfo  = dir(fullfile(fldPath, '*.mat'));
        matNames = {matInfo.name};

        % split into conc_* vs. velocity files
        isConc   = startsWith(matNames, 'conc_');
        velFiles = {'flowCache.mat','tmpvelresults.mat'};
        isVel    = ismember(matNames, velFiles);

        concNames = matNames(isConc);
        velNames  = matNames(isVel);

        % if user asked for a specific concTest, filter it
        if filterConc
            matchIdx = strcmp(concNames, concTest);
            if ~any(matchIdx)
                warning('No concentration file "%s" found in folder %s', concTest, fld);
                concNames = {};  % nothing to load
            else
                concNames = concNames(matchIdx);
            end
        end

        % load concentration data
        concStruct = struct();
        for i = 1:numel(concNames)
            fprintf('\r[%s] Loading concentration %d/%d: %s', ...
                    fld, i, numel(concNames), concNames{i});
            data = load(fullfile(fldPath, concNames{i}));
            f    = genvarname(concNames{i}(1:end-4));  % strip “.mat”
            concStruct.(f) = data;
        end
        fprintf('\n');

        % load velocity data
        velStruct = struct();
        for i = 1:numel(velNames)
            fname = velNames{i};
            data  = load(fullfile(fldPath, fname));
            f     = genvarname(fname(1:end-4));
            velStruct.(f) = data;
        end

        % safe field‐name for this folder
        safeFld = genvarname(fld);

        % store into results
        results.(safeFld) = struct( ...
            'Concentration', concStruct, ...
            'Velocity',    velStruct ...
        );
    end
end

function varargout=generate_MaSIV_list(stitchedDir,overide)
% Make MaSIV files listing the stitched file locations
%
% function status=generate_MaSIV_list(stitchedDir,overide)
%
%
% Purpose
% Make stitched image file lists from data stitched with tvMat to
% enable import into MaSIV.
%
%
% Inputs
% stitchedDir - path to stitched data directory. e.g. 'stitchedImages_100'
% overide - [optional, 0 by default] if 1, build channel lists even if sections are missing
% 
% 
% Outputs
% status - [optional]
%           -1 - failed to build anything
%            0 - made partial list but there are missing sections
%            1 - made complete list
%
%
% Examples
% generate_MaSIV_list('stitchedImages_100')
% 
% Rob Campbell
%
% Notes: 
% Produces file with unix file seps on all platforms. Windows
% MATLAB seems OK about using these paths. Windows fileseps 
% mess up the fprintf.

if nargin<2
    overide=0;
end

if ispc
    fprintf('Fails on Windows machines. Not fixed yet\n')
end

if ~exist(stitchedDir,'dir')
    fprintf('Directory %s not found\n',stitchedDir)
    return
end



%Generate the file names we will use 
p=readMetaData2Stitchit;
stitchedFileListName=[p.sample.ID,'_ImageList_'];
MaSIV_DIR = [p.sample.ID,'_MaSIV'];
MaSIV_YML = [p.sample.ID,'_Meta.yml'];


if ~exist(MaSIV_DIR,'dir')
    mkdir(MaSIV_DIR)
end

%Write the YML file
YML.VoxelSize.x = 1; %p.voxelSize.x; %TODO: change this when MaSIV is fixed
YML.VoxelSize.y = 1; %p.voxelSize.y;
YML.VoxelSize.z = p.voxelSize.Z;
YML.stackName = p.sample.ID;
YML.imageBaseDirectory = '../';

stitchit.tools.writeSimpleYAML(YML,fullfile(MaSIV_DIR,MaSIV_YML))



%find the channels
chans = dir(stitchedDir);

for ii=1:length(chans)
    chanName = chans(ii).name;
    if startsWith(chanName, '.')
        % skip private and relative folder
        continue
    end
    tifDir = fullfile(stitchedDir, chanName);
    if ~isdir(tifDir)
        % skip if there are local files. Continue only for directories
        continue
    end

    tifs=dir(fullfile(tifDir,'*.tif'));

    if isempty(tifs)
        fprintf('No tiffs in %s. Skipping\n',tifDir)
        continue
    end

    missing=findMissingSections(tifs);
    if missing && ~overide
        fprintf(['\nMissing sections. Not building the image lists.\n',...
            'Please fix your data or overide this warning (help %s), if you know what you''re doing.\n\n'], ...
            mfilename)
        if nargout>0
            varargout{1}=-1;
        end
        return
    end
    if missing && overide
        fprintf('\n BUILDING THE LISTS WITH MISSING SECTIONS\n\n')

    end


    fprintf('Making channel %s file\n',chanName)

    thisFname = fullfile(MaSIV_DIR, sprintf('%sCh_%s.txt',stitchedFileListName,chanName) );
    fid=fopen(thisFname,'w+');
    for thisTif = 1:length(tifs)         
            fprintf(fid,[tifDir,'/',tifs(thisTif).name,'\n']);

    end
    fclose(fid);



end


if nargout>0
    varargout{1}=~missing;
end




% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
function missing=findMissingSections(tifs)
    %Look for missing sections 
    sections = zeros(length(tifs),1);
    optical = zeros(length(tifs),1);    


    for ii=1:length(tifs)
        tok=regexp(tifs(ii).name,'section_(\d+)_(\d+)','tokens');
        if isempty(tok)
            error('regexp failed')
        end

        tok=tok{1};
        if length(tok)~=2
            error('Failed to find two tokens')
        end

        sections(ii) = str2num(tok{1});
        optical(ii) = str2num(tok{2});              
    end

    sections=unique(sections);
    optical=unique(optical);

    missing=0;

    %now check if any sections arem missing (a bit brute-force, but it'll work)
    for sct=1:length(sections)
        for opt=1:length(optical)

            f=find(sections==sct);
            if isempty(f)
                fprintf('\t ** Missing physical section %d, optical section %d **\n',sct,opt)
                missing=1;
            end

        end
    end



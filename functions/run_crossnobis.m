RUNSUBS            = 1;
rundevice = 'bear'; % 'Desktop'

subjects           = {};
s                  = {'SS_080795','SP_190590','JP_310391','SM_210492','AB_110492',...
                      'SG_031196','AM_190496','AD_140696','TC_021196','NM_210598',...
                      'ST_260381','YH_080891','KV_290695','ZL_260395','KO_090596',...
                      'IG_240693','AM_240196','RA_181089','YC_180787','AA_291197'};


switch rundevice
    case 'Desktop'
        gpfs  = fullfile('c:\Science\');
        rdsfs = fullfile('Y:\');
        gpfssubjectdir     = fullfile(gpfs,'irsaeeg','Data');
        fileatts = '-h +w';
        subjectdir         = fullfile(rdsfs,'irsaeeg','Data');
        addpath(genpath(fullfile(gpfs,'software/fieldtrip-20160105')));
        copyf          = false;
    case 'bear'
        gpfs  = '/gpfs/bb/charesti/';
        rdsfs = '/rds/projects/2016/charesti-01/';
        gpfssubjectdir     = fullfile(gpfs,'nbu','irsaeeg','Data');
        fileatts = '+w';
        subjectdir         = fullfile(rdsfs,'irsaeeg','Data');
        addpath(genpath(fullfile(gpfs,'software/fieldtrip-20160105')));
        copyf          = true;    
end


for i = 1:numel(s)
    c         = pwd;
    cd(fullfile(subjectdir, s{i}))
    eval(sprintf('subjects{i} = %s_info;',s{i}))
    cd(c)
end

cfg = {};    
for i = RUNSUBS
    
    sub          = subjects{i};
    if ~exist(fullfile(gpfssubjectdir, sub.sID, 'classified'),'dir')
        mkdir(fullfile(gpfssubjectdir, sub.sID, 'classified'))
        fileattrib(fullfile(gpfssubjectdir, sub.sID, 'classified'),fileatts,'','s');
    end
    cfg2   = [];
    cfg2.ch=1:128;
    cfg2.conds=1:72;
    cfg2.saveto         = fullfile(gpfssubjectdir, sub.sID, 'classified','rdms_xnobis');
    cfg2.dataset        = fullfile(gpfssubjectdir, sub.sID, 'preprocessed','ALL_preprocessed_NoAR_concat');
    cfg2.makefig        = 1;
    cfg    = [cfg cfg2];
end
cellfun(@crossnobis, cfg);
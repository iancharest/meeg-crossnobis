

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
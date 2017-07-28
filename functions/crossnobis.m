function crossnobis(cfg)

try
    
    load(cfg.dataset)
    ch = cfg.ch;
    nch = numel(ch);
    conds  = cfg.conds;
    nconds = numel(conds);
    C = indicatorMatrix('allpairs',conds);
    
    times = data.time{1};
    ntimes = numel(times);
    ntrials = numel(squeeze(IC_selectTrials(data,1,1,1,1)));
    
    
    % now simulate a design matrix
    Xtrain = [kron(eye(nconds),ones(ntrials-1,1)),ones(nconds*(ntrials-1),1)];
    Xtest  = [eye(nconds), ones(nconds,1)];
    
    trains_ind = flipud(nchoosek(1:ntrials,ntrials-1))';
    
    RDMs = nan(ntrials,ntimes,size(C,1));
    
    thisDat=[];
    
    for condI=conds
        thisDat = cat(4,thisDat,squeeze(IC_selectTrials(data,condI,1,ch,1:ntimes)));
    end
    %
    % thisDat = permute(thisDat,[1 2 4 3]);
    %
    % % z transform across trials
    % gavg = repmat(nanmean(thisDat,4), 1,1,1,size(thisDat,4));
    % gstd = repmat(std(thisDat,0,4), 1,1,1,size(thisDat,4));
    % thisDat = permute((thisDat-gavg)./gstd,[4 2 3 1]);
    
    thisDat = permute(thisDat,[3 2 4 1]);
    
    for trialI=1:ntrials % leave one out loop
        clc
        fprintf('*** performing loo %d -- %3.2f%% completed ***\n',trialI,(trialI/ntrials)*100);
        
        Btrain   = nan(nch,nconds+1,ntimes);
        Btest    = nan(nch,nconds+1,ntimes);
        
        Restrain = nan(nch,(nconds*ntrials)-nconds,ntimes);
        Restest  = nan(nch,nconds,ntimes);
        for channelI = ch
            
            Ytrain = reshape(permute(thisDat(trains_ind(:,trialI),:,:,channelI),[2 1 3]),[ntimes (ntrials-1)*nconds])';
            Ytest  = permute(thisDat(trialI,:,:,channelI),[3 2 1]);
            
            Btrain(channelI,:,:) = pinv(Xtrain)*Ytrain;
            Btest(channelI,:,:)  = pinv(Xtest)*Ytest;
            
            % check what it looks like
            %Yhat = X*squeeze(Btrain(channelI,:,:)); % Yhat is the model
            Restrain(channelI,:,:)  = Ytrain - Xtrain*squeeze(Btrain(channelI,:,:)); % Yhat is the model
            %Restest(channelI,:,:)   = Ytest - Xtest*squeeze(Btest(channelI,:,:)); % Yhat is the model
        end
        
        for timeI=1:ntimes
            
            SwTrain         = covdiag(Restrain(:,:,timeI)')^0.5;
            
            %SwTest          = covdiag(Restest(:,:,timeI)')^0.5;
            
            % compute LDCs
            wA=(C*Btrain(:,conds,timeI)')/SwTrain;
            wB=(C*Btest(:,conds,timeI)')/SwTrain;
            LDC=sum(wB.*wA,2);%/nch;
            
            RDMs(trialI,timeI,:) = LDC;
        end
    end
    
    if cfg.saveto
        save(cfg.saveto,'RDMs','times','-v7.3');
    end
    
    
catch ME
    fid = fopen('error_report.txt', 'a');
    fprintf(fid,'%s crossnobis error: %s\n%s\n\n', datestr(clock),cfg.dataset, getReport(ME));
    fclose(fid);
end


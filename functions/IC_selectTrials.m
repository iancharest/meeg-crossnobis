function data = IC_selectTrials(data, cond,col,ch,times)
   % convenience function for fieldtrip structures
   % Takes one data structure, the condtion to separate and on
   % which column (in the data.trialinfo) the condition variable can be found
   triali       = find(data.trialinfo(:,col) == cond);
   
   if isempty(triali)
       error('Could not find any trials under condition %d on column %d!', cond, col)
   else
       data           = catcell(3,data.trial(triali));
   end
   data = data(ch,times,:);
end


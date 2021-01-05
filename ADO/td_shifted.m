function vVar = td3(trialnum, choice, ID, subdir)
    % ADO module for dedlay discounting using a shifted hyperbolic model which adds a term "a" to the logistic function 1/(1+exp(e*(Vss - Vll- a)))
    trialnum = floor(trialnum)-1;
    %1 or 0, o is left
    choice = floor(choice);

    fVar = struct('numpractice',4,'numparmest',30,'aHBmin',0,'aHBmax',8,...
                    'aHBgran',80,'aHB',[],'aHBparms',[],'kHBmin',log(0.00001),...
                    'kHBmax',log(0.9),'eHBmin',0.01,'eHBmax',1,'kHBgran',80,...
                    'eHBgran',29,'kHB',[],'eHB',[],'kHBparms',[],'eHBparms',[],... 
                    'nSamp',[],'gtable',[],'nnn',[],'ddisp',[],'filler',[]);
                
    vVar = struct('HBprior',[],'HBpost',[],'index',[],'indexlist',[],...
                    'design',[],'datastore',[],'practicedatastore',[],...
                    'mkHB',[],'meHB',[],'sdkHB',[],'sdeHB',[],'kHB90',[],...
                    'eHB90',[],'f',[],'fyo',[],'maHB',[],'sdaHB',[],'aHB90',[]);

    if trialnum==-1 || trialnum==-2
        
        % Initialize Parameters
        fVar.kHB = [0 linspace(fVar.kHBmin,fVar.kHBmax, fVar.kHBgran+1)];
        fVar.eHB = [0 linspace(fVar.eHBmin,fVar.eHBmax, fVar.eHBgran+1)];
        fVar.aHB = [0 linspace(fVar.aHBmin,fVar.aHBmax, fVar.aHBgran+1)];
        
        [fVar.kHBparms, fVar.eHBparms, fVar.aHBparms] = meshgrid(fVar.kHB(2:end), fVar.eHB(2:end),fVar.aHB(2:end));

        fVar.kHBparms = permute(fVar.kHBparms,[2,1,3]);
        fVar.eHBparms = permute(fVar.eHBparms,[2,1,3]);
        fVar.aHBparms = permute(fVar.aHBparms,[2,1,3]);
        
        fVar.nSamp = 1000;
        
        c = (fVar.eHBgran+1)*(fVar.kHBgran+1)*(fVar.aHBgran+1);
        vVar.HBprior=ones(fVar.kHBgran+1,fVar.eHBgran+1,fVar.aHBgran+1)/c;
        vVar.HBpost = vVar.HBprior;
        
        % initialize design space
        fidx = fopen('design800_jc.txt','rt');
        A = textscan(fidx,'%f %f %f %f %f', 'delimiter', '\n','headerlines',0,'collectoutput',true);

        fVar.gtable = A{1};
        fclose(fidx);
        clear A;
        
        fVar.nnn = length(fVar.gtable);
        vVar.datastore = zeros(fVar.numparmest,7);
        vVar.practicedatastore = zeros(fVar.numpractice,7);
        
        vVar.mkHB = zeros(fVar.numparmest,1);
        vVar.meHB = zeros(fVar.numparmest,1);
        vVar.maHB = zeros(fVar.numparmest,1);
        
        vVar.sdkHB = zeros(fVar.numparmest,1);
        vVar.sdeHB = zeros(fVar.numparmest,1);
        vVar.sdaHB = zeros(fVar.numparmest,1);
        
        vVar.kHB90 = zeros(fVar.numparmest,2);
        vVar.eHB90 = zeros(fVar.numparmest,2);
        vVar.aHB90 = zeros(fVar.numparmest,2);
        
        if ~ischar(ID)
            ID = [num2str(ID),'/'];
        else
            ID = [ID, '\'];
        end
        
        if ~exist([subdir,'/',ID],'dir')
            mkdir([subdir,'/',ID]);
        end
        
        fVar.ddisp = {[  0,    0.143,  0.286,   0.429,  0.571,  0.714,   0.857,...
               1,      1.143,  1.286,   1.428,  1.571,  1.714,   1.857,    2, 2.143,...
               2.286,  2.571,   3,      3.5,    4,       4.5,      5,...
               6,      7,       8.57,   10.714, 12.86,   15,       17.14,...
               21.43,  25.714,  30,     34.286, 38.57,   42.86,    47.14,...
               52.14,  78.214,  104.286,130.36, 156.43,  182.5,    208.57,...
               260.714,312.86,  365,    417.14, 469.286, 521.43],...
               {'now', 'in 1 day', 'in 2 days', 'in 3 days', 'in 4 days', 'in 5 days',...
               'in 6 days','in 7 days', 'in 8 days',  'in 9 days', 'in 10 days',...
               'in 11 days', 'in 12 days','in 13 days', 'in 14 days','in 15 days', 'in 16 days',...
               'in 18 days','in 3 weeks','in 3 weeks and a half','in 4 weeks',...
               'in 4 weeks and a half','in 5 weeks','in 6 weeks','in 7 weeks','in 2 months',...
               'in 2 months and a half', 'in 3 months', 'in 3 months and a half',...
               'in 4 months', 'in 5 months', 'in 6 months', 'in 7 months', 'in 8 months',...
               'in 9 months', 'in 10 months','in 11 months','in 1 year', 'in 1 year and a half',...
               'in 2 years','in 2 years and a half', 'in 3 years', 'in 3 years and a half',...
               'in 4 years','in 5 years','in 6 years','in 7 years','in 8 years','in 9 years',...
               'in 10 years'}};
           
        fVar.filler=[   750, 0.143, 800, 21.43;
                        700, 0.143, 800, 52.14;
                        730, 0.143, 800, 104.286;
                        650, 0.143, 800, 521.43;
                        750, 0.143, 800, 260.714;
                        400, 0.143, 800, 1;
                        500, 0.143, 800, 0.714;
                        300, 0.143, 800, 6;
                        550, 0.143, 800, 0.429;
                        450, 0.143, 800, 2];
        
        vVar.index = randi([0,fVar.nnn-1]);
        vVar.design = fVar.gtable(vVar.index+1,2:end);
        
        mix = rand() + .5;
        if mix>=1
            dtemp = vVar.design(1);
            vVar.design(1) = vVar.design(3);
            vVar.design(3) = dtemp;
            dtemp = vVar.design(2);
            vVar.design(2) = vVar.design(4);
            vVar.design(4) = dtemp;
        end
        obs =[];
        save([subdir,'/',ID,'/data.mat'],'vVar','fVar','obs');
        
        if trialnum ==-2
            fprintf('%1.0f practice & %1.0f estimation trials',fVar.numpractice,fVar.numparmest);

            return
        end
        
    else
        load([subdir,'/',ID,'/data.mat']);
        
        % store the results of the last trial
        if vVar.design(1) <  vVar.design(3)
            tempy = choice;
        else
            tempy = 1- choice;
            temp = vVar.design(1:2);
            vVar.design(1:2) = vVar.design(3:4);
            vVar.design(3:4) = temp;
        end
        
        if trialnum+1 > fVar.numpractice
            currTrial = trialnum -  fVar.numpractice;
            vVar.datastore(currTrial+1,7)=tempy;
            vVar.datastore(currTrial+1,1:6) = [currTrial vVar.index vVar.design];
        else
            vVar.practicedatastore(trialnum+1,7) = tempy;
            vVar.practicedatastore(trialnum+1,1:6) = [trialnum vVar.index vVar.design];
        end
        
        % update parameter estimates
        if trialnum+1 > fVar.numpractice
            vVar.HBprior = vVar.HBpost;
            

            uSSmat = vVar.design(1).*(( exp(fVar.kHBparms).*vVar.design(2)+1).^-1);
            uLLmat = vVar.design(3).*(( exp(fVar.kHBparms).*vVar.design(4)+1).^-1);
            
            
            py1theta = (1+exp(fVar.eHBparms.*(uSSmat - uLLmat-fVar.aHBparms))).^-1; 
            
            if tempy ==1
                temp = py1theta .* vVar.HBprior;
            else
                temp = (1-py1theta) .* vVar.HBprior;
            end
            
            vVar.HBpost = temp / sum(sum(sum(temp)));
        end
        
        if trialnum+1> fVar.numpractice
            fid = fopen([subdir,'/',ID,'/obs.txt'],'a');
            fprintf(fid, '%3.0f %4.0f %6.2f %5.2f %6.2f %5.2f %2.0f\n',vVar.datastore(currTrial+1,1),...
                            vVar.datastore(currTrial+1,2),vVar.datastore(currTrial+1,3),...
                            vVar.datastore(currTrial+1,4),vVar.datastore(currTrial+1,5),...
                            vVar.datastore(currTrial+1,6),vVar.datastore(currTrial+1,7));
            fclose(fid);
        else
            fid = fopen([subdir,'/',ID,'/practobs.txt'],'a');
            fprintf(fid, '%3.0f %4.0f %6.2f %5.2f %6.2f %5.2f %2.0f\n',vVar.practicedatastore(trialnum+1,1),...
                            vVar.practicedatastore(trialnum+1,2),vVar.practicedatastore(trialnum+1,3),...
                            vVar.practicedatastore(trialnum+1,4),vVar.practicedatastore(trialnum+1,5),...
                            vVar.practicedatastore(trialnum+1,6),vVar.practicedatastore(trialnum+1,7));
            fclose(fid);
        end
        CredProb = [0.05 0.95];
        
        if trialnum+1 > fVar.numpractice
            sumHBpost3 = sum(sum(vVar.HBpost,2),1);
            sumHBpost2 = sum(sum(vVar.HBpost,2),3);
            sumHBpost1 = sum(sum(vVar.HBpost,1),3);
            ksumpost = reshape(sumHBpost2,1, fVar.kHBgran+1);
            esumpost = reshape(sumHBpost1,1, fVar.eHBgran+1);
            asumpost = reshape(sumHBpost3,1, fVar.aHBgran+1);
            kHBpost = [0 ksumpost];
            eHBpost = [0 esumpost];
            aHBpost = [0 asumpost];
            
            vVar.mkHB(currTrial+1,1) = sum(kHBpost .* fVar.kHB);
            vVar.meHB(currTrial+1,1) = sum(eHBpost .* fVar.eHB);
            vVar.maHB(currTrial+1,1) = sum(aHBpost .* fVar.aHB);
            
            vVar.sdkHB(currTrial+1,1) = sqrt(sum(((fVar.kHB-vVar.mkHB(currTrial+1,1)).^2).*kHBpost));
            vVar.sdeHB(currTrial+1,1) = sqrt(sum(((fVar.eHB-vVar.meHB(currTrial+1,1)).^2).*eHBpost));
            vVar.sdaHB(currTrial+1,1) = sqrt(sum(((fVar.aHB-vVar.maHB(currTrial+1,1)).^2).*aHBpost));
            
            tempsum = cumsum(kHBpost);
            idx1 = sum(tempsum <= CredProb(1));
            idx2 = length(tempsum) - sum(tempsum >= CredProb(2));
            
            vVar.kHB90(currTrial+1,1) = fVar.kHB(idx1) + (fVar.kHB(idx1+1)-fVar.kHB(idx1))* ...
               (CredProb(1)-tempsum(idx1))/kHBpost(idx1+1); 
            vVar.kHB90(currTrial+1,2) = fVar.kHB(idx2) + (fVar.kHB(idx2+1)-fVar.kHB(idx2))* ...
               (CredProb(2)-tempsum(idx2))/kHBpost(idx2+1);
           
            
            tempsum = cumsum(eHBpost);
            idx1 = sum(tempsum <= CredProb(1));
            idx2 = length(tempsum) - sum(tempsum >= CredProb(2));
            
            vVar.eHB90(currTrial+1,1) = fVar.eHB(idx1) + (fVar.eHB(idx1+1)-fVar.eHB(idx1))* ...
               (CredProb(1)-tempsum(idx1))/eHBpost(idx1+1); 
            vVar.eHB90(currTrial+1,2) = fVar.eHB(idx2) + (fVar.eHB(idx2+1)-fVar.eHB(idx2))* ...
               (CredProb(2)-tempsum(idx2))/eHBpost(idx2+1);
           
            tempsum = cumsum(aHBpost);
            idx1 = sum(tempsum <= CredProb(1));
            idx2 = length(tempsum) - sum(tempsum >= CredProb(2));
            
            vVar.aHB90(currTrial+1,1) = fVar.aHB(idx1) + (fVar.aHB(idx1+1)-fVar.aHB(idx1))* ...
               (CredProb(1)-tempsum(idx1))/aHBpost(idx1+1); 
            vVar.aHB90(currTrial+1,2) = fVar.aHB(idx2) + (fVar.aHB(idx2+1)-fVar.aHB(idx2))* ...
               (CredProb(2)-tempsum(idx2))/aHBpost(idx2+1);
           
            fid = fopen([subdir,'/',ID,'/HBpost_SummaryStats.txt'],'w');
            fprintf(fid,'%25s','Estimation of parameter k\n');
            fprintf(fid, '%7.8s %8.7s %7.7s %10.7s %7.7s\n', 'trial', 'mean','stdev', 'lower90','upper90');
            for i = 1:currTrial
                fprintf(fid,'%7.0f %8.4f %7.4f %10.4f %7.4f\n', i, vVar.mkHB(i,1), vVar.sdkHB(i,1), vVar.kHB90(i,1), vVar.kHB90(i,2)); 
            end
            
            fprintf(fid,'%25s','Estimation of parameter e\n');
            fprintf(fid, '%7.8s %8.7s %7.7s %10.7s %7.7s\n', 'trial', 'mean','stdev', 'lower90','upper90');
            for i = 1:currTrial
                fprintf(fid,'%7.0f %8.4f %7.4f %10.4f %7.4f\n', i, vVar.meHB(i,1), vVar.sdeHB(i,1), vVar.eHB90(i,1), vVar.eHB90(i,2)); 
            end
            
            fprintf(fid,'%25s','Estimation of parameter a\n');
            fprintf(fid, '%7.8s %8.7s %7.7s %10.7s %7.7s\n', 'trial', 'mean','stdev', 'lower90','upper90');
            for i = 1:currTrial
                fprintf(fid,'%7.0f %8.4f %7.4f %10.4f %7.4f\n', i, vVar.maHB(i,1), vVar.sdaHB(i,1), vVar.aHB90(i,1), vVar.aHB90(i,2)); 
            end
            
            fclose(fid);
    
        end
        
    end
 
    if trialnum+1 <= fVar.numpractice
        
        vVar.index = randi([0,fVar.nnn]);
        vVar.design = fVar.gtable(vVar.index+1,2:end);
        
        mix = rand() + .5;
        if mix>=1
            dtemp = vVar.design(1);
            vVar.design(1) = vVar.design(3);
            vVar.design(3) = dtemp;
            dtemp = vVar.design(2);
            vVar.design(2) = vVar.design(4);
            vVar.design(4) = dtemp;
        end
    else
        [paramidx, cnt] = categrnd3d(vVar.HBpost, fVar.nSamp);
        tmp_kHBp = reshape(permute(fVar.kHBparms,[2,1,3]), 1, numel(fVar.kHBparms));
        tmp_eHBp = reshape(permute(fVar.eHBparms,[2,1,3]), 1, numel(fVar.eHBparms));
        tmp_aHBp = reshape(permute(fVar.aHBparms,[2,1,3]), 1, numel(fVar.aHBparms));
        
        tmp_kHBp = tmp_kHBp(paramidx);
        tmp_eHBp = tmp_eHBp(paramidx);
        tmp_aHBp = tmp_aHBp(paramidx);
        numidx = length(cnt);
       
        vSS=fVar.gtable(:,2);  
        tSS=fVar.gtable(:,3);  
        vLL=fVar.gtable(:,4); 
        tLL=fVar.gtable(:,5);  
       
        idx1 = ones(fVar.nnn,1);
        idx2 = ones(1, numidx);
        
        uSS = repmat(vSS,1,numidx).*(tSS.* exp(tmp_kHBp)+1).^-1;
        uLL = repmat(vLL,1,numidx).* (tLL.* exp(tmp_kHBp)+1).^-1;
        
        py1theta = (1+exp((uSS-uLL-repmat(tmp_aHBp,fVar.nnn,1)).*repmat(tmp_eHBp,fVar.nnn,1))).^-1; % sigshift 
        py1theta(py1theta<1e-10)=1e-10;
        py1theta(py1theta>1-1e-10)=1-1e-10;
        
        py1 = (py1theta*cnt)/fVar.nSamp;
        py0 = 1- py1;
        entR = -py1.*log(py1) - py0.*log(py0);
        
        py0theta = 1-py1theta;
        entRcond = ((-py1theta.*log(py1theta)-py0theta.*log(py0theta))*cnt)/fVar.nSamp;
        
        utilmat = entR - entRcond;
        
        [~, loc] = max(utilmat);
        
        vVar.index = loc-1;
        vVar.indexlist = [vVar.indexlist loc-1];
        if (currTrial>2) && ((vVar.indexlist(currTrial+1) == vVar.indexlist(currTrial)) ||...
                (vVar.indexlist(currTrial+1) == vVar.indexlist(currTrial-1)))
            vVar.fyo=1;
            fx = ceil(randi([1, length(fVar.filler)]));
            vVar.f = fVar.filler(fx,:);
            mix = rand() + .5;
            if mix>= 1
                dtemp = vVar.design(1);
                vVar.design(1) = vVar.design(3);
                vVar.design(3) = dtemp;
                dtemp = vVar.design(2);
                vVar.design(2) = vVar.design(4);
                vVar.design(4) = dtemp;
            end
            vVar.design = fVar.gtable(loc,2:5);
            c_arr = [-20 -10 0 10 20];
            vVar.design(1)= vVar.design(1)+c_arr(randi([1 5]));

            if vVar.design(1)<=0
                vVar.design(1) = 10;
            end
            if vVar.design(1) >= vVar.design(3)
                vVar.design(1) = vVar.design(3)-10;
            end
            
            ind_fut = find(fVar.ddisp{1}==vVar.design(4))-1;
            
            if ind_fut == 1
                vVar.design(4) = fVar.ddisp{1}(ind_fut +2);
            elseif ind_fut == length(fVar.ddisp{1})-1
                vVar.design(4) = fVar.ddisp{1}(ind_fut+1);
            else
                ind_fut2 = ind_fut + datasample([1 -1],1);

                vVar.design(4) = fVar.ddisp{1}(ind_fut2+1);
            end
                
        else
            vVar.fyo = 0;
            vVar.f=[];
            vVar.design = fVar.gtable(loc,2:5);
        end
        
        mix = rand() + .5;
        if mix >= 1
            dtemp = vVar.design(1);
            vVar.design(1) = vVar.design(3);
            vVar.design(3) = dtemp;
            dtemp = vVar.design(2);
            vVar.design(2) = vVar.design(4);
            vVar.design(4) = dtemp;
        end   
    end


    if  vVar.fyo ==1
        ind_d1 = find(fVar.ddisp{1}==vVar.f(2))-1;
        ind_d2 = find(fVar.ddisp{1}==vVar.f(4))-1;
        fprintf('$%1.0f %s & $%1.0f %s & \n',vVar.f(1), fVar.ddisp{2}{ind_d1+1},...
                    vVar.f(3),fVar.ddisp{2}{ind_d2+1});
        
        ind_d1 = find(fVar.ddisp{1}==vVar.design(2))-1;
        ind_d2 = find(fVar.ddisp{1}==vVar.design(4))-1;
        fprintf('$%1.0f %s & $%1.0f %s & \n',vVar.design(1), fVar.ddisp{2}{ind_d1+1},...
                    vVar.design(3),fVar.ddisp{2}{ind_d2+1});
    else
        ind_d1 = find(fVar.ddisp{1}==vVar.design(2))-1;
        ind_d2 = find(fVar.ddisp{1}==vVar.design(4))-1;
        fprintf('$%1.0f %s & $%1.0f %s & \n',vVar.design(1), fVar.ddisp{2}{ind_d1+1},...
                    vVar.design(3),fVar.ddisp{2}{ind_d2+1});
    end
    save([subdir,'/',ID,'/data.mat'],'vVar','fVar','obs');
    
    
    
    
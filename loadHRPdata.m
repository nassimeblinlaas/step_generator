function uLINK = loadHRPdata(fname)
%%%% parse openhrp data and make uHRP data
%%%% 2004 Mar.3 s.kajita AIST

%fname = 'HRP2main_full.wrl';
%fname = 'HRP2LRmain.wrl';

fid = fopen(fname);

nest_level = 0;
jointdef = 0;
idx = 1;
while 1
    tline = fgetl(fid);
    if ~ischar(tline), break, end
    
    %% cut comments
    if findstr('#',tline)
      tline = tline(1:findstr('#',tline)-1);
    end

    %% nest level of {}
    nest_level = nest_level + nest_count(tline);
    
    %% find joint definition
    [T1,rest] = strtok(tline);
    
    if strcmp(T1,'DEF')
        [T2,rest] = strtok(rest);
        [T3,rest] = strtok(rest);
        if strcmp(T3,'Joint')
            jointdef = 1; 
            %disp(T2);
            idx = idx+1;
            LINKwk(idx).name = T2;
            if strcmp(T2, 'WAIST')
                LINKwk(idx).mother = 0;
                jid_stack(1) = 0;
                jid_stack(2) = 1;
                LINKwk(idx).a = [0 0 1]';
            else
                LINKwk(idx).mother = jid_stack(nest_level-1);
            end
        elseif ~strcmp(T3,'Segment')
            jointdef = 0;
        end
    elseif jointdef
        %% get joint property
        if strcmp(T1,'jointId')
            jid = str2num(rest)+2;    % jointId
            jid2idx(jid) = idx;       % table from jointId -> idx
            jid_stack(nest_level) = jid;
        elseif strcmp(T1,'mass')
            LINKwk(idx).m = str2num(rest); 
        elseif jointdef & strcmp(T1,'translation')
            LINKwk(idx).b = str2num(rest)';
        elseif strcmp(T1,'jointAxis')
            if findstr('X',rest)
                LINKwk(idx).a = [1 0 0]';
            elseif findstr('Y',rest)
                LINKwk(idx).a = [0 1 0]';
            elseif findstr('Z',rest)
                LINKwk(idx).a = [0 0 1]';
            end
        elseif strcmp(T1,'centerOfMass')
            LINKwk(idx).c = str2num(rest)';
        elseif strcmp(T1,'momentsOfInertia')
            wk = str2num(rest);
            Iwk = reshape(wk,3,3);
            LINKwk(idx).I = (Iwk + Iwk')/2;
        elseif strcmp(T1,'gearRatio')
            LINKwk(idx).gr = str2num(rest);
        elseif strcmp(T1,'rotorInertia')
            LINKwk(idx).Ir = str2num(rest);
        end
    end
end

fclose(fid);

%%% sort link data by jointId
rsize = length(jid2idx);

%global uLINK
uLINK = LINKwk(2);

for n=2:rsize
    idx = jid2idx(n);
    uLINK(n) = LINKwk( idx );
end

%%%%%%% find child and sister %%%%%%

%% initialize
for n=1:rsize
    uLINK(n).sister = 0;
    uLINK(n).child  = 0;
    uLINK(n).q      = 0;
    uLINK(n).dq     = 0;
    uLINK(n).ddq    = 0;
    uLINK(n).vertex = [];
end

uLINK(1).p = [0,0,0]';
uLINK(1).R = eye(3);


%% who am i sequence
for n=1:rsize
    mom = uLINK(n).mother;          % I know you are my mother.
    if mom ~= 0 
        if uLINK(mom).child == 0
            uLINK(mom).child = n;           % Mother, I'm your daughter !!
        else
            elder_sister = uLINK(mom).child;    % I have elder sister!
            while(uLINK(elder_sister).sister ~= 0)
                elder_sister = uLINK(elder_sister).sister;  % I have another elder sister!
            end
            uLINK(elder_sister).sister = n;   % I'm your younger sister!
        end
    end
end

        
    
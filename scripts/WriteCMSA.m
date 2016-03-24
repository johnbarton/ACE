%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright for this implementation:
% 2015 - Simona Cocco
% cocco@lps.ens.fr
%
% Permission is granted for anyone to copy, use, or modify this
% software and accompanying documents for any uncommercial
% purposes, provided this copyright notice is retained, and note is made
% of any changes that have been made. This software and
% documents are distributed without any warranty, express or
% implied. All use is entirely at the user's own risk.
%
%
% Any publications resulting from applications of ACE should cite:
%
% J P Barton, E de Leonardis, A Coucke, and S Cocco (2015)
% ACE: adaptive cluster expansion for maximum entropy
% graphical model inference
% and 
% S Cocco and R Monasson (2011). Adaptive Cluster Expansion for 
% Inferring Boltzmann Machines with Noisy Data. Physical Review 
% Letters 106, 090601.
%
%
% This function takes as input a general alignment of configurations of a
% system (multi-sequence alignment for proteins) and constructs the outputs
% necessary for the ACE algorithm and generative tests of the inferred model.
%
%
% List of inputs:
%
% -filetype (string, default: "numbers")
%   This specifies the format of the input alignment.
%   Three alignment formats are possible:
%   "numbers" - The alignment consists of a matrix with B rows (the sequences)
%       and N columns (the sites). States are identified by numbers.
%   "letters" - The alignment consists of a matrix with B rows (the sequences)
%       and N columns (the sites). States are identified by letters.
%   "fasta" - The alignment is in the standard fasta format.
%
% -filename (string)
%   String with the filename of the input alignment. The three formats above
%   are supported.
%
% -theta (real number >= 0, generally from 0 to 0.3)
%   The reweighting threshold, used to take into account correlated sampling
%   (phylogeny). If zero there is no reweighting. Otherwise it weights the
%   contribution of each sequence to the average frequencies and correlations
%   with a weight that is inversely proportional to the number of sequences
%   in the MSA with Hamming distance smaller than theta * N. A typical value
%   for a MSA is 0.2.
%
% -redmethod (string: "frequency" or "entropy")
%   If "entropy" the reduction of Potts states is based on the single site entropy
%   contribution, otherwise it is based on the frequency of the Potts state.
%
% -redcut (real number >= 0, <=1)
%   Reduction threshold: if 0 only states which are never observed are removed.
%
% -gapred (0 or 1)
%   When gapred=1, gaps in the alignment are replaced with other a.a. according
%   to the frequency of a.a. at the same site in other sequences (without gaps
%   at that site).
%   NOTE: this replacement is stochastic, so the resulting set of correlations
%   will not always be identical when this option is used.
%
% -gaugechoice (string: "least", "cons", "wt", "group"),
%   Choice of the gauge state for writing out the correlations.
%   "least": The least frequent a.a. at each site.
%   "cons": The consensus a.a.
%   "wt": One can give as input a wild-type sequence taken as the gauge;
%       in this case a file wt.dat with the same length of the MSA has to be in
%       the directory.
%   "group": The grouped Potts state is set as the gauge state.
%   For a different gauge choice it is sufficient to give the sequence
%   corresponding to this choice as the wild-type ("wt").
%
%
% Outputs:
%
% -filename.p
%   The frequencies and pairwise correlations required as input for ace.
%
% -filename.cmsa
%   A compressed form of the alignment used by qgt.
%
% -filename.cons
%   The consensus sequence in the reduced alphabet, used by qgt.
%
% -filename.wgt
%   The corresponding weight for each sequence or configuration in the data
%   after the reweighting procedure, used by qgt.
%
% - filename.errj
%   The elementary variance of the inferred parameters in the same format
%   as the correlations (".p" file).
%   Error propagation is needed to calculate error bars in other gauges.
%
% -filename.freq
%   A list of single site frequencies in a vector with N*q rows, with no reduction.
%
% -filename.ind
%   Indices of the selected states in the N*q alphabet.
%
% -filename.indgauge
%   The gauge sequence in the N*q alphabet.
%
% -filename.indgrouped
%   Indices of the regrouped states in the 21-state alignment.
%
% -filename.qgrouped
%   The number of grouped states for each site.
%
% -filename.q
%   The number of effective Potts states for each site, equal to the number
%   of explicitly modeled states if there are no regrouped states, or the
%   number of explicitly modeled states + 1 otherwise.
%
%
% EXAMPLES IN THE PAPER
%
% Erdos-Renyi random graph model ER05 ("examples/ER05-configurations.dat")
% WriteCMSA('numbers','ER05-configurations.dat',0,'frequency',0.05,'group',0)
%
% Lattice protein S_B ("examples/lattice-protein-SB-alignment.faa"):
% NOTE: To reproduce test correlations, change the function letter2number to the
% lattice protein version (see line 714).
% WriteCMSA('letters','lattice-protein-SB-alignment.dat',0,'frequency',0,'least',0)
%
% Trypsin inhibitor PF00014 ("examples/PF00014-alignment.faa"):
% reweighting=0.2, frequency reduction, p0=0, gap substitution, gauge least
% WriteCMSA('fasta','PF00014-alignment.faa',0.2,'frequency',0,'least',1);
% reweighting=0.2, frequency reduction, p0=0.05, gap substitution, gauge cons
% WriteCMSA('fasta','PF00014-alignment.faa',0.2,'frequency',0.05,'cons',1);
%
% HIV protein p7 ("examples/p7-alignment.fasta")
% reweighting=0, entropy reduction, S0=0.9, gap substitution, gauge cons
% WriteCMSA('fasta','p7-alignment.fasta',0,'entropy',0.9,'cons',1);
%
% Cortical neuron recording ()

function WriteCMSA(filetype, filename, theta, redmethod, redcut, gauge, gapred)
    disp('Reading the alignment and reweighting...')
    tic

    if (strcmp(filetype,'fasta'))
        X = fastaread(filename);
        N = size(X(1).Sequence,2);
        M = size(X,1);
        alignc = zeros(M,N);
        for i=1:M
            for j=1:N
                alignc(i,j) = letter2number(X(i).Sequence(j));
            end
        end

    elseif (strcmp(filetype,'letters'))
        X = importdata(filename);
        N = size(X{1},2);
        M = size(X,1);
        alignc = zeros(M,N);
        for i=1:M
            for j=1:N
                alignc(i,j) = letter2number(X{i}(j));
            end
        end

    elseif (strcmp(filetype,'numbers'))
        name = filename;
        alignc = dlmread(name);
        N = size(alignc,2);
        M = size(alignc,1);
    
    end

    q = max(max(alignc)); % maximal number of Potts states

    msa21gap = zeros(M,N*q);
    aa = zeros(1,N);
    for i=1:N;
        for a=1:q; % a.a. are labeled from 1 to q
            msa21gap(:,(i-1)*q+a) = (alignc(:,i)==a);
        end
    end
    
    pgap = sum(msa21gap)/M; % p = frequency of amino acid at site i
    p2cgap = reshape(pgap,q,N);
 
    for i=1:N;
        aa(i) = size(find(p2cgap(:,(i))),1);
    end
    
    align=zeros(M,N);
    if(gapred)
        for i=1:M
            for j=1:N
                if((alignc(i,j)==1)&&(aa(j)>2))
                    % substitution of the gap state with another a.a. with the
                    % probability of the a.a. for the non-gapped sequences
                    % n.b. only if there are at least 2 a.a.+gap with
                    % non zero frequencies at the site
 
                    indreg = find(p2cgap(:,j)>0);
                    indre = indreg(2:size(indreg));
                    cumd = cumsum(p2cgap(indre,j));
                    rd = random('uniform',0,cumd(size(cumd,1)));
                    an = find(rd<cumd, 1 );
                    align(i,j) = indre(an);
                    
                else
                    align(i,j)=alignc(i,j);
                    
                end
            end
        end
        
    else
        align = alignc;
        
    end


    % reweight the alignment with the parameter theta
    % reweighting with C++ mexfile if the alignment is bigger
    % than 200 sequences.
    % % to comment after the first run
    % mex weightCalculator.c

    W = ones(1,M);
    if( theta > 0.0 )
        if ( M<200 )
            W = (1./(1+sum(squareform(pdist(align,'hamm')<theta))));
        else
            align_vec = zeros(1,N*M);
            for i=1:M
                align_vec( 1+(i-1)*N : i*N ) = align( i, : );
            end
            W = weightCalculator(align_vec, theta, M, N);
            clear align_vec
        end
    end
    
    Meff = sum(W);
    w = W'/Meff;

    disp(['Effective number of sequences is B=',num2str(Meff),'.'])
    toc


    % alignment msa21 for N*21 sites
    msa21 = zeros(M,N*q);
    for i=1:N;
        for a=1:q;  %a.a. are labeled from 1 to q
            msa21(:,(i-1)*q+a) = (align(:,i)==a);
        end
    end

    p = sum(msa21)/M; % p = frequency of amino acid at site i
    p2c = reshape(p,q,N);
 
    [pmax,indcon] = max(p2c);
 
    if (strcmp(gauge,'wt'));
        indgauge = load('wt.dat');
    elseif (strcmp(gauge,'cons'));
        [pgauge,indgauge] = max(p2c);
    elseif (strcmp(gauge,'group'));
        indgauge = zeros(1,N);
    else
        indpmin = zeros(1,N);
        for i=1:N
            indnz = find(p2c(:,i)>redcut);
            [pmin(i),indpmin(i)] = min(p2c(indnz,i));
            indgauge(i) = indnz(indpmin(i));
        end
    end

    % pgauge(:,1) = p2c(:indgauge);
    % matrix of probabilities qXN


    disp('Compressing states...');
    tic;

    [qa,qng,ind,qout,indout,cons] = colorreduction(redmethod,redcut,q,p2c,N,indgauge,indcon);

    qtng = zeros(N+1,1);
    qt = zeros(N+1,1);
    qtout = zeros(N+1,1);
    qtout(1) = 0;
 
    for i=2:N+1
        qtng(i) = qtng(i-1)+qng(i-1);
        qt(i) = qt(i-1)+qa(i-1);
        qtout(i) = qtout(i-1)+qout(i-1);
    end
 
    % NN: effective size of the msacut
    % cut the MSA according to state reduction and compute reweighted
    % probabilities
    % msacut: matrix with only selected states

    msacut = zeros(M,qt(N+1));
    for i=1:N
        if (strcmp(gauge,'group'))
            msacut(:,qt(i)+1:qt(i)+qng(i)) = msa21(:,ind(qtng(i)+1:qtng(i)+qng(i))');
        else
            if (qout(i)>0)
                msacut(:,(qt(i)+1:qt(i)+qng(i))) = msa21(:,ind(qtng(i)+1:qtng(i)+qng(i))');
                msacut(:,qt(i)+qa(i)) = sum(msa21(:,indout(qtout(i)+1:qtout(i)+qout(i))),2)>0;
            else
                msacut(:,qt(i)+1:qt(i)+qng(i)) = msa21(:,ind(qtng(i)+1:qtng(i)+qng(i))');
            end
        end
    end

    NN = size(msacut,2);
    toc;
    
    
    % one and two point statistics for ace input file
    disp('Computing one- and two-point correlations...');
    tic;
 
    wm = repmat(w,[1,NN]);
    
    % reweighted probabilities
    pw = sum(msacut.*wm);
    pm = repmat(pw,[M,1]);
    
    % substraction of the mean to the probabilities

    dq = sparse(msacut-pm);
    ww = spdiags(w,0,M,M);
    csp = dq'*ww*dq;
    c = full(csp); % c: reweighted correlation matrix
 
%     % substraction of the mean to the probabilities
%     dq = msacut-pm;
%  
%     % c: reweighted correlation matrix
%     if(theta>0)
%         ww = diag(w);
%         c = dq'*ww*dq;
%     else
%         c = (dq'*dq)./M;
%     end
    
    toc
    
    
    disp('Writing output for ace and qgt');
    tic

    % Write correlations and estimated errors on parameters in sce format
    %if((strcmp(redmethod,'entropy'))
    %firstname=[filename(1:find(filename=='.',1,'last')-1),'_w',num2str(theta)
    
    prename = [filename(1:find(filename=='.',1,'last')-1)];
    if(theta>0)
        prename = [prename,'_w',num2str(theta)];
    end
    
    if(strcmp(redmethod,'entropy'));
        prename = [prename,'_fS',num2str(redcut),'_',gauge];
    else
        prename = [prename,'_po',num2str(redcut),'_',gauge];
    end
    
    if (gapred)
        prename = [prename,'_nogap'];
    end
    
    % write the reweighting vector in the file 'filename_theta_pmin.wgt'
    name = [prename,'.wgt'];
    fid = fopen(name,'w');
    for i=1:M
        fprintf(fid,'%e\n',W(i));
    end
    fclose(fid);

    % write indices for a.a. kept and regrouped
    name = [prename,'.ind'];
    fid = fopen(name,'w');
    for i=1:size(ind,1)
        fprintf(fid,'%d\n',ind(i));
    end
    fclose(fid);

    name = [prename,'.indgrouped'];
    fid = fopen(name,'w');
    for i=1:size(indout,1)
        fprintf(fid,'%d\n',indout(i));
    end
    fclose(fid);

    name = [prename,'.p'];
    %name2=[prename,'.dhdj'];
    name3 = [prename,'.errj'];
    fid = fopen(name,'w');
    %fid2 = fopen(name2,'w');
    fid3 = fopen(name3,'w');

    la = zeros(N,q);
    l = 0;
    for i=1:N
        for a=1:qa(i)
            l = l+1;
            la(i,a) = l;
            fprintf(fid,' %e',pw(l));
            %fprintf(fid2,' %e',(1-pw(l))/((pw(l)+1/Meff)*Meff));
        end
        fprintf(fid,'\n');
        %fprintf(fid2,'\n');
    end
    
    for i=1:N
        pi = pw(qt(i)+1:qt(i)+qa(i))';
        for j=i+1:N
            pj = pw(qt(j)+1:qt(j)+qa(j));
            pp2 = c(qt(i)+1:qt(i)+qa(i),qt(j)+1:qt(j)+qa(j))+pi*pj;
            
            ppic = pi-sum(pp2,2);
            ppicmat = repmat(ppic,[1,qa(j)]);
            ppcj = pj-sum(pp2,1);
            ppcjmat = repmat(ppcj,[qa(i),1]);
 
            ppcc = 1-sum(pi,1)-sum(ppcj,2);
            mone = ones(qa(i),qa(j));
            ppccmat = ppcc.*mone;
 
            reg = 1/(10*Meff).*mone;
            pp2n = pp2+reg;
            ppicmatn = ppicmat+reg;
            ppcjmatn = ppcjmat+reg;
            ppccmatn = ppccmat+reg;
            s1 = (mone-pp2)./pp2n;
            s2 = (mone-ppicmat)./ppicmatn;
            s3 = (mone-ppcjmat)./ppcjmatn;
            s4 = (mone-ppccmat)./ppccmatn;
            errpg = sqrt(s1+s2+s3+s4)./sqrt(Meff);
            errp = s1./Meff;
 
            for a=1:qa(i)
                for b=1:qa(j)
                    fprintf(fid,' %e',max(pp2(a,b),1e-10));
                    %fprintf(fid2,' %e',errp(a,b));
                    fprintf(fid3,' %e',errpg(a,b));
                end
            end
            
            fprintf(fid,'\n');
            %fprintf(fid2,'\n');
            fprintf(fid3,'\n');
            
        end
    end
    
    fclose(fid);
    %fclose(fid2);
    fclose(fid3);


    % write out the MSA in compact format
    name = [prename,'.cmsa'];
    fid = fopen(name,'w');
    
    for i=1:size(msacut,2)
        fprintf(fid,'%d\n',-1);
        for j=1:size(msacut,1)
            if (msacut(j,i)==1)
                fprintf(fid,'%d\n',j);
            end
        end
    end
    fclose(fid);


    % record the consensus index in the reduced states alphabet
    name = [prename,'.cons'];
    name2 = [prename,'.gauge'];
    name3 = [prename,'.indgauge'];
    fid = fopen(name,'w');
    fid3 = fopen(name3,'w');
    for i=1:N
        fprintf(fid,'%d\n',cons(i));
        fprintf(fid3,'%d\n',indgauge(i));
    end
    fclose(fid);
    fclose(fid3);

    % record the number of regruped a.a. for each site
    name = [prename,'.qgrouped'];
    fid = fopen(name,'w');
    for i=1:N
        fprintf(fid,'%d\n',qout(i));
    end
    fclose(fid);

    name = [prename,'.q'];
    fid = fopen(name,'w');
    for i=1:N
        fprintf(fid,'%d\n',qa(i));
    end
    fclose(fid);

    % record the frequencies for each a.a. in the original alignment (no reweighting no reduction)
    name = [prename,'.freq'];
    fid = fopen(name,'w');
    for i=1:N*q
        fprintf(fid,'%2.6f\n',p(i));
    end
    fclose(fid);
    
    toc

%Function which builds the indices of the reduced alphabet
function[qa,qng,ind,qout,indout,cons] = colorreduction(redmethod,redcut,q,p2c,N,indgauge,indmax)
    
    %ind gives the entries of the selected a.a. in the 21*N vector
    %indout = the a.a. grouped in the gauge state
    %build ind on the first sites
    %redmethod='frequency'
    %indmax=indcon;
 
    ind = [];
    indout = [];
    
    % entropy reduction
    if(strcmp(redmethod,'entropy'));
        S_ia = min(-p2c.*log(p2c),p2c~=0);
        S = sum(S_ia);
        [A,sorted] = sort(S_ia,'descend');
        qnonzero = zeros(N,1);
        cons = zeros(N,1);
 
        for ii=1:1;
            qnonzero(ii) = size(find(A(:,ii)),1);
            Sred = 0;
            Sreds = 0;
            k = 0;

            while((k<1)||(Sred/S(ii)<redcut))
                k = k+1;
                Sreds = Sreds+min(-p2c(sorted(k,ii),ii).*log(p2c(sorted(k,ii),ii)),p2c(sorted(k,ii),ii)~=0);
                pout = sum(p2c(sorted(k+1:qnonzero(ii),ii),ii));
                Sredr = min(-pout.*log(pout),pout~=0);
                Sred = Sreds+Sredr;
            end
            
            indt1 = sorted(1:k,ii);
            indoutt = sorted(k+1:qnonzero(ii),ii);
            indt = indt1(indt1~=indgauge(ii));
            cons(ii) = sum(find(indt==indmax(ii)));
            indoutt = indoutt(indoutt~=indgauge(ii));
        end
 
        qa = zeros(1,N);
        qout = zeros(N,1);
        qng = zeros(1,N);
        qng(1) = size(indt,1);
        qout(1) = size(indoutt,1);
        ind = indt(1:qng(1));
        indout = indoutt(1:qout(1));
 
        if(qout(1)>0);
            qa(1) = qng(1)+1;
        else
            qa(1) = qng(1);
        end
 
        %build ind on other sites
        for ii=2:N;
            qnonzero(ii) = size(find(A(:,ii)),1);
            Sred = 0;
            Sreds = 0;
            k = 0;
            
            while((k<1)||(Sred/S(ii)<redcut))
                k = k+1;
                Sreds = Sreds+min(-p2c(sorted(k,ii),ii).*log(p2c(sorted(k,ii),ii)),p2c(sorted(k,ii),ii)~=0);
                pout = sum(p2c(sorted(k+1:qnonzero(ii),ii),ii));
                Sredr = min(-pout.*log(pout),pout~=0);
                Sred = Sreds+Sredr;
            end
 
            indt1 = sorted(1:k,ii);
            indoutt = sorted(k+1:qnonzero(ii),ii);
            indt = indt1(indt1~=indgauge(ii));
            cons(ii) = sum(find(indt==indmax(ii)));
            indoutt = indoutt(indoutt~=indgauge(ii));
            qout(ii) = size(indoutt,1);
            qng(ii) = size(indt,1);
            indt2 = indt(1:qng(ii))+(ii-1)*q;
            
            if(qng(ii)>0)
                ind = [ind;indt2];
            end
                
            indoutt2 = indoutt(1:qout(ii))+(ii-1)*q;

            if(qout(ii)>0);
                qa(ii) = qng(ii)+1;
                indout = [indout;indoutt2];
            else
                qa(ii) = qng(ii);
            end
        
        end
 
    % frequency reduction
    else
        cons = zeros(1,N);
        indt1 = find((p2c(:,1)>redcut));
        indt = indt1(indt1~=indgauge(1));
        cons(1) = sum(find(indt==indmax(1)));
 
        indoutt = find((p2c(:,1))&(p2c(:,1)<=redcut));
        indoutt = indoutt(indoutt~=indgauge(1));
 
        % qa: number of potts states at each site taking into account
        % grouped state (last potts state)
        % qng number of non-grouped potts states at each site
        qa = zeros(1,N);
        qout = zeros(1,N);
        qng = zeros(1,N);
        
        %first site
        qng(1) = size(indt,1);
        qout(1) = size(indoutt,1);
        ind = indt(1:qng(1));
        indout = indoutt(1:qout(1));
 
        if(qout(1)>0);
            qa(1)=qng(1)+1;
        else
            qa(1)=qng(1);
        end
        
        %build ind on other sites
        for ii=2:N
            indt1 = find((p2c(:,ii)>redcut));
            indt = indt1(indt1~=indgauge(ii));
            cons(ii) = sum(find(indt==indmax(ii)));
            qng(ii) = size(indt,1);
            indt2 = indt(1:qng(ii))+(ii-1)*q;
 
            if(qng(ii)>0)
                ind = [ind;indt2];
            end
            
            indoutt = find((p2c(:,ii))&(p2c(:,ii)<=redcut));
            indoutt = indoutt(indoutt~=indgauge(ii));
            qout(ii) = size(indoutt,1);
            indoutt2 = indoutt(1:qout(ii))+(ii-1)*q;

            if(qout(ii)>0);
                qa(ii) = qng(ii)+1;
                indout = [indout;indoutt2];
            else
                qa(ii) = qng(ii);
            end
            
        end
        
    end
    
    if(max(indgauge(:))==0);
        qa = qng;
    end

end


% Auxiliary function

function x=letter2number(a)
    switch(a)
        case '-'
            x=1;
        case 'A'
            x=2;
        case 'C'
            x=3;
        case 'D'
            x=4;
        case 'E' 
            x=5;
        case 'F'
            x=6;
        case 'G' 
            x=7;
        case 'H'
            x=8;
        case 'I' 
            x=9;
        case 'K'
            x=10;
        case 'L' 
            x=11;
        case 'M'
            x=12;
        case 'N' 
            x=13;
        case 'P'
            x=14;
        case 'Q'
            x=15;
        case 'R'
            x=16;
        case 'S' 
            x=17;
        case 'T'
            x=18;
        case 'V'
            x=19;
        case 'W'
            x=20;
        case 'Y'
            x=21;
        otherwise
            x=1;
        end

% For the lattice protein case we have used a different alphabet, which is
% the following, the only difference is the order of the Potts states in the
% input .p file
% Alphabet used in the lattice protein model

%function x=letter2number(a)
%    % full AA alphabet
%    switch(a)
%        case 'C'
%            x=1;
%        case 'M'
%            x=2;
%        case 'F'
%            x=3;
%        case 'I'
%            x=4;
%        case 'L'
%            x=5;
%        case 'V'
%            x=6;
%        case 'W'
%            x=7;
%        case 'Y'
%            x=8;
%        case 'A'
%            x=9;
%        case 'G'
%            x=10;
%        case 'T'
%            x=11;
%        case 'S'
%            x=12;
%        case 'N'
%            x=13;
%        case 'Q'
%            x=14;
%        case 'D'
%            x=15;
%        case 'E'
%            x=16;
%        case 'H'
%            x=17;
%        case 'R'
%            x=18;
%        case 'K'
%            x=19;
%        case 'P'
%            x=20;
%        end
    
    
    end
end
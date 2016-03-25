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
% This function takes as input a file giving the times (or time bins in
% which a variable (neuron) is active) and constructs the outputs necessary
% for the ACE algorithm and generative tests of the inferred model.
%
% NOTE: For data other than neural recordings, use WriteCMSA.m instead.
%
%
% List of inputs:
%
% -filetype (string: "list")
%   This specifies the format of the input data.
%   Two formats are possible:
%   "list" - A list of act times bins in which each variable is active,
%       in the ".cmsa" format
%   "neuro" - A neural recording of activity times and corresponding neurons.
%       The time bin size must be entered in units of the recording time.
%
% -filename (string)
%   String with the filename of the input alignment. The three formats above
%   are supported.
%
% -theta (real number >= 0, generally from 0 to 0.3)
%   The reweighting threshold, used to take into account correlated sampling.
%   If zero there is no reweighting. Otherwise it weights the contribution of
%   each configuration to the average frequencies and correlations with a weight
%   that is inversely proportional to the number of configurations in the data
%   with Hamming distance smaller than theta * N.
%   For neural data a typical choice is 0 (no reweighting).
%
% -bin (real number >=0)
%   The time bin size for neural recordings, in units of the recording time.
%
%
% List of outputs:
%
% -filename.p
%   The frequencies and pairwise correlations required as input for ace.
%
% -filename.cmsa
%   A compressed form of the data used by qgt.
%
% -filename.cons
%   The consensus configuration, used by qgt.
%
% -filename.wgt
%   The corresponding weight for each configuration in the data
%   after the reweighting procedure (if any), used by qgt.
%
% - filename.errj
%   The elementary variance of the inferred parameters in the same format
%   as the correlations (".p" file).
%   Error propagation is needed to calculate error bars in other gauges.
%
% -filename.freq
%   A list of single site frequencies in a vector with N*q rows, with no reduction.
%
% 
% EXAMPLES IN THE PAPER
%
% Cortical neuron recording ("examples/neural-recording.dat")
% Time units in the recording are 0.1ms, and the time bin is 10ms.
% WriteCMSAbin('neuro','neural-recording.dat',0,100);
% Computing the correlations based on the output ".cmsa" file.
% WriteCMSAbin('list','neural-recording.cmsa',0,0);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function WriteCMSAbin(filetype,filename,theta,bin)
    disp('Reading the data and reweighting...')
    tic;

    if((strcmp(filetype,'list')))
        cmsa = load(filename);
        align = spike2msa(cmsa);
    else
        spiketime = load(filename);
        spiketimemin = min(spiketime,[],1);
        spiketimemax = max(spiketime,[],1);

        ttot = size(spiketime,1);
        tbin = floor((spiketime(:,1)-spiketimemin(1))/bin)+1;
        neuro = spiketime(:,2);
        smsa = sparse(tbin,neuro,1);

        %to remove multiple spikes in the same time bin to one
        smsa = min(smsa,1);
        align = full(smsa);
    end

    N = size(align,2);
    disp(['Configuration length is N=',num2str(N),'.'])
    M = size(align,1);
    q = max(max(align)); % maximal number of Potts states

    % reweight the alignment with the parameter theta
    % reweighting with C++ mexfile if the alignment is bigger than 200
    % sequences.
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


    % frequency
    p = sum(align)/M; % p = frequency of amino acid at site i
    NN = size(align,2);
    toc;
    
    
    % one and two point statistics for ace input file
    disp('Computing one- and two-point correlations...');
    tic;
     
    wm = repmat(w,[1,NN]);
    
    % reweighted probabilities
    pw = sum(align.*wm);
    pm = repmat(pw,[M,1]);
    
    % substraction of the mean to the probabilities
    dq = sparse(align-pm);

    if(theta>0);
        ww = spdiags(w,0,M,M);
        csp = dq'*ww*dq;
    else
        csp=dq'*dq./M ;
    end
 
    c = full(csp); % c = reweighted correlation matrix
    toc
  
    
    disp('Writing output for ace and qgt');
    tic

    % Write correlations and estimated errors on parameters in sce format

    prename = [filename(1:find(filename=='.',1,'last')-1)];
    if(theta>0)
        prename = [prename,'_w',num2str(theta)];
    end

    % write the reweighting vector in the file 'filename_theta_pmin.wgt'
    name = [prename,'.wgt'];
    fid = fopen(name,'w');
    for i=1:M
        fprintf(fid,'%e\n',W(i));
    end
    fclose(fid);

    %write the file .p and the  approximate error-bars
    name = [prename,'.p'];
    name2 = [prename,'.errj'];
    
    fid = fopen(name,'w');
    fid2 = fopen(name2,'w');

    la = zeros(N,q);
    l = 0;
    for i=1:N
        fprintf(fid,' %e',pw(i));
        fprintf(fid2,' %e',(1-pw(i))/((pw(i)+1/Meff)*Meff));
        fprintf(fid,'\n');
        fprintf(fid2,'\n');
    end
    
    for i=1:N
        pi = pw(i)';
        for j=i+1:N
            pj = pw(j);
            pp2 = c(i,j)+pi*pj;
            reg = 1/(10*Meff);
            pp2n=pp2+reg;
            ppic=pi-pp2;
            ppjc=pj-pp2;
            ppicn = ppic+reg;
            ppjcn = ppjc+reg;
            ppcc = 1-pi-pj+pp2;
            ppccn=ppcc+reg;
            s1 = (1-pp2)./pp2n;
            s2 = (1-ppic)./ppicn;
            s3 = (1-ppjc)./ppjcn;
            s4 = (1-ppcc)./ppccn;
            errpg = sqrt(s1+s2+s3+s4)./sqrt(Meff);
            errp = s1./Meff;
 
            fprintf(fid,' %e',max(pp2,1e-10));
            fprintf(fid2,' %e',errpg);
            fprintf(fid,'\n');
            fprintf(fid2,'\n');
        end
    end
    
    fclose(fid);
    fclose(fid2);
    

    % write out the MSA in compact format
    name = [prename,'.cmsa'];
    fid = fopen(name,'w');
    
    for i=1:size(align,2)
        fprintf(fid,'%d\n',-1);
        for j=1:size(align,1)
            if (align(j,i)==1)
                fprintf(fid,'%d\n',j);
            end
        end
    end
    fclose(fid);


    % record the consensus index in the reduced states alphabet
     name = [prename,'.cons'];
     fid = fopen(name,'w');

     for i=1:N
         fprintf(fid,'%d\n',0);
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


    % Function to read msa in compact form; -1 separates variables and the data
    % corresponds to configuration in which that variable is present

    function msa = spike2msa(spikes);

        indb = find(spikes>0);
        tbin = spikes(indb);
        ind = find(spikes<0);
        neuro = [];

        N = size(ind,1);

        for i=1:N-1;
            nbi = ones(ind(i+1)-ind(i),1).*i;
            neuro = [neuro;nbi];
        end

        nbi = ones(size(tbin,1)-ind(size(ind,1))+1,1).*N;
        neuro = [neuro;nbi];
        smsa = sparse(tbin,neuro,1);

        %count multiple spikes in the same time bin as one
        smsa = min(smsa,1);
        msa = full(smsa);

    end

end

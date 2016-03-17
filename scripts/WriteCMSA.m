%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Copyright for this implementation:
%2015 - John Barton, Eleonora de Leonardis, R. Monasson, Simona Cocco
%jpbarton@gmail.com,eleonora.deleonardis@gmail.com,monasson@lpt.ens.fr, cocco@lps.ens.fr

%Permission is granted for anyone to copy, use, or modify this
%software and accompanying documents for any uncommercial
%purposes, provided this copyright notice is retained, and note is % made
%of any changes that have been made. This software and
%documents are distributed without any warranty, express or
%implied. All use is entirely at the user's own risk.

%Any publication resulting from applications of ACE should cite:

   %  J Barton, E de Leonardis, R Monasson, and S Cocco (2015)
    % ACE: inference of graphical models describing functional constraints
    % on protein sequences and neural activity patterns
    % and  
    % S Cocco, R Monasson. Adaptive Cluster Expansion for Inferring
    % Boltzmann Machines with Noisy Data. Physical Review Letters 106, 090601
    % (2011)



%This function takes as imput a general alignment of configurations of a system
% (multi-sequence alignment for proteins) and  built up the outputs necessary
% for the ACE algorithm and  the generative test of the model.

% List of inputs:

% - filetype = string (one of the those reported above) it specifies the format
%          of the given alignment.
%   Three alignment formats are possible:
%     "numbers" - the alignment consists of a matrix with B columns
%                 (the sequences) and N rows (the sites). Each a.a. 
%                 if identified by a number
%     "letters" - the alignment consists of a matrix with B columns
%                 (the sequences) and N rows (the sites). Each a.a. 
%                 if identified by a letter
%     "fasta"   - the alignment is in the usual fasta format.
% No specification means "numbers".

% - filename: string with the filename  of the input alignment. The three formats
%              above are supported.

%-theta: numerical value  generally from 0 to 0.3. The reweighting threshold to take into account correlated sampling of configuration/sequences.
%If it is zero there is no reweighting.
%Otherwise it counts each sequence in the averages frequencies and correlations with a weigth inversly proportional to  the number of sequences  in the MSA at an Hamming distance smaller than theta. 
% Typical value for a MSA =0.2

% -redmethod: string: 'frequency' or 'entropy', if entropy the reduction of Potts symbols is based on the 
%site entropy contribution otherwise it is based on the frequency of the Potts symbol.

% redcut: numerical value from 0 to 1,
% reduction threshold if = 0  only symbols which are never present on a site are cut from the  Potts symbols. 

%-gapred: binary value 0,1
%When gapred=1 substitution of gap with a.a. using the statistic of non-gapped sequences
% no reduction if gapred=0

%gaugechoice:string: 'least', 'cons', 'wt'
%The gauge can be choosen to: 'least': the least frequent a.a. on each site, 'cons': the consensus a.a. or 
%'wt' one can  give as input a whild-type sequence taken as the gauge.  
%in this case a file wt.dat with the same length of the MSA has to be in
%the directory
%for a different gauge choice it is sufficient to give the sequence
%corresponding to this choice as wt

%Outputs:
% - filename.cmsa = it is the alignment in a compact format useful for the Gentest.
%-1 sepatate Potts states in the reduced binary representation of the alignment and
%numbers give sequences in which that Potts state is present
% - filename.cons = it is the consensus sequence in the reduced alphabet (useful for Gentest).
% - filename.errap = it is the approximate statistical error on the inferred
%parameters in the same format as the .p.
%-filename.freq = single site frequencies in a vector with N*q columns
%-filename.ind = indexes of the selected symbols in the N*q alphabet
%-filename.indgauge =the gauge sequence in the N*q alphabet
% -filename.p = frequencies and correlations needed for the ace program 
%-filename.indgrouped= indexes of the regrouped symbols on the 21-symbol alignement
% -filename.qgrouped the number of the grouped symbols for each site
%-filename.q the number of sffetive Potts  symbols for each site (= nb. selected symbols if   
%there are not regrouped symbols and nb. selected+1 if there are regrouped
%symbols

% Examples with the input alignment file  PF00014-2014.faa :
%  reweighting=0.2,frequency reduction, pmin=0.05, gapsubstitution, gauge least 
%WriteCMSA('fasta','PF00014-2014.faa',0.2,'frequency',0.05,'least',1);
% reweighting=0.2,frequency reduction, pmin=0, gapsubstitution, gauge least
%WriteCMSA('fasta','PF00014-2014.faa',0.2,'frequency',0,'least',1);
%Example with a lattice-protein alignment file align_B.faa
%WriteCMSA('letters','align_B.faa',0,'frequency',0,'least',0)




function WriteCMSA(filetype, filename, theta, redmethod, redcut,gauge,gapred)
disp('Read alignment and reweighting...')
    tic






    if (strcmp(filetype,'fasta')) 
     X = fastaread(filename);
     N = size(X(1).Sequence,2);
     M = size(X,1);
     alignc = zeros(M,N);
    
     for i=1:M
      for j=1:N
        alignc(i,j)=letter2number(X(i).Sequence(j));
        
      end
     end

    elseif (strcmp(filetype,'letters'))
      
     X=importdata(filename);
     N=size(X{1},2);
     M=size(X,1);
     alignc = zeros(M,N);
    
     for i=1:M
       
        for j=1:N
                
          alignc(i,j)=letter2number(X{i}(j));
        end
     end
    
    
    

    elseif (strcmp(filetype,'numbers'))
     name=filename;
     alignc=dlmread(name);
     N=size(alignc,2);
     M=size(alignc,1);

    end
%q=maximal number of Potts states
q=max(max(alignc));


msa21gap=zeros(M,N*q);
aa=zeros(1,N);
for i=1:N;
  
    for a=1:q;
        
%  a.a. are labeled from 1 to q
    
        msa21gap(:,(i-1)*q+a)=(alignc(:,i)==a);
    end
    
end

 pgap=sum(msa21gap)/M; % p = frequence of amino acid  on site i
 p2cgap=reshape(pgap,q,N);
 
 
 
for i=1:N;
   aa(i)=size(find(p2cgap(:,(i))),1);
    
    
    
end
align=zeros(M,N);   
if(gapred)
    for i=1:M
    
    for j=1:N
       
        if((alignc(i,j)==1)&&(aa(j)>2))
            %substitution of the gap symbol with another a.a. with the
            %probability of the a.a. on the non-gapped sequences
            %n.b. only if there are at least 2 a.a.+gap with 
            % non zero frequencies on the site
            
            indreg=find(p2cgap(:,j)>0);
            indre=indreg(2:size(indreg));
            cumd=cumsum(p2cgap(indre,j));
            rd=random('uniform',0,cumd(size(cumd,1)));
            an=find(rd<cumd, 1 );
            align(i,j)=indre(an);
        else
            align(i,j)=alignc(i,j);
            
        
        end
        
    end
    end
else
    align=alignc;
end
% reweigth the alignment with the parameter theta 
% reweighting with C++ mexfile if the alignment is bigger than 200
% sequences.
%  % to comment after the first run
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
Meff=sum(W);
 w=W'/Meff;



toc


%alignment msa21 on N*21 sites

msa21=zeros(M,N*q);
for i=1:N;
   
    
    for a=1:q;
        
%  a.a. are labeled from 1 to q
    
        msa21(:,(i-1)*q+a)=(align(:,i)==a);
    end
    
end

 
 p=sum(msa21)/M; % p = frequence of amino acid  on site i
 p2c=reshape(p,q,N);
 
 [pmax,indcon]=max(p2c);
 
  if   (strcmp(gauge,'wt'));
              indgauge=load('wt.dat');

  
  
  
  
  elseif  (strcmp(gauge,'cons'));
      [pgauge,indgauge]=max(p2c);
      
      
   else 
      
         indpmin=zeros(1,N); 
 
     for i=1:N
     indnz=find(p2c(:,i)>redcut);
     [pmin(i),indpmin(i)]=min(p2c(indnz,i));
     
     indgauge(i)=indnz(indpmin(i));
    end

   end	
            
              %pgauge(:,1)=p2c(:indgauge);
              
 % matrix pf probabilities qXN
 
 disp('Compression...');
    tic ;
[qa,qng,ind,qout,indout,cons]=colorreduction(redmethod,redcut,q,p2c,N,indgauge,indcon);




 qtng=zeros(N+1,1);

 qt=zeros(N+1,1);
 
 qtout=zeros(N+1,1);
 qtout(1)=0;
 
 
 for i=2:N+1
 qtng(i)=qtng(i-1)+qng(i-1);
 qt(i)=qt(i-1)+qa(i-1);
 qtout(i)=qtout(i-1)+qout(i-1);
 end
  
 
 %NN: effective size of the msacut
 % cut the MSA according to color reduction and compute reweighted
	% probabilities
%msacut: matrix with only selected symbols


  msacut=zeros(M,qt(N+1)); 
  for i=1:N
      if (qout(i)>0) 
msacut(:,(qt(i)+1:qt(i)+qng(i)))=msa21(:,ind(qtng(i)+1:qtng(i)+qng(i))');
msacut(:,qt(i)+qa(i))=sum(msa21(:,indout(qtout(i)+1:qtout(i)+qout(i))),2)>0;
      else
   
         msacut(:,qt(i)+1:qt(i)+qng(i))=msa21(:,ind(qtng(i)+1:qtng(i)+qng(i))');
      end
  end
NN=size(msacut,2);
    toc;
 % one and two point statistic for input file to the inverse problem 
    disp('One and 2 Point Statistics...');
    tic;
   

 
    wm=repmat(w,[1,NN]);
    %reweighted probabilities
    pw=sum(msacut.*wm);
 
     pm=repmat(pw,[M,1]);
    
    % substraction of the mean to the probabilities
   
    dq=msacut-pm;
    
    % c :reweighted correlation matrix
    ww=diag(w);
    c=dq'*ww*dq;
    
   
  
    
toc
disp('Writing outputs for ACE and Gentest');
    tic

%write magnetisations and correlation and  estilated errors on parameters in sce format

%if((strcmp(redmethod,'entropy'))
%firstname=[filename(1:find(filename=='.',1,'last')-1),'_w',num2str(theta)
prename=[filename(1:find(filename=='.',1,'last')-1)];
if(theta>0) 
prename=[prename,'_w',num2str(theta)];
end
if(strcmp(redmethod,'entropy'));
prename=[prename,'_fS',num2str(redcut),'_',gauge];
else
prename=[prename,'_po',num2str(redcut),'_',gauge];
end
if (gapred)
prename=[prename,'_nogap'];
end
%write the reweighting vector in the file 'filename_theta_pmin.wgt' 
name=[prename,'.wgt'];
fid = fopen(name,'w');

for i=1:M
    fprintf(fid,'%e\n',W(i));
end
fclose(fid);

%write indexes on a.a. kept and regrouped

name=[prename,'.ind'];

fid = fopen(name,'w');

for i=1:size(ind,1)
    fprintf(fid,'%d\n',ind(i));
end
fclose(fid);

name=[prename,'.indgrouped'];

fid = fopen(name,'w');

for i=1:size(indout,1)
    fprintf(fid,'%d\n',indout(i));
end
 fclose(fid);


name=[prename,'.p'];

name2=[prename,'.errapp'];

fid = fopen(name,'w');
fid2 = fopen(name2,'w');



la=zeros(N,q);
l=0;
for i=1:N   
    for a=1:qa(i)
        l=l+1;
        la(i,a)=l;
        %fprintf(fid,' %e %i',pw(lif (pw(l)>0)),1);
        fprintf(fid,' %e',pw(l));
        fprintf(fid2,' %e',sqrt((1-pw(l))/((pw(l)+1/Meff)*Meff))+ sqrt((1-pmax(i))/((pmax(i)+1/Meff)*Meff)));
    end
    fprintf(fid,'\n');
    fprintf(fid2,'\n');
end
for i=1:N
pi=pw(qt(i)+1:qt(i)+qa(i))';
    for j=i+1:N
 
       pj=pw(qt(j)+1:qt(j)+qa(j)) ;
       pp2=c(qt(i)+1:qt(i)+qa(i),qt(j)+1:qt(j)+qa(j))+pi*pj;
       
       ppic=pi-sum(pp2,2);
       ppicmat=repmat(ppic,[1,qa(j)]);
       
       
       ppcj=pj-sum(pp2,1);
        ppcjmat=repmat(ppcj,[qa(i),1]);
       ppcc=1-sum(pj)-sum(ppcj);
       mone=ones(qa(i),qa(j));
       ppccmat=ppcc.*mone;
       
       
        reg=1/(10*Meff).*mone;
        pp2n=pp2+reg;
        ppicmatn=ppicmat+reg;
        ppcjmatn=ppcjmat+reg;
        ppccmatn=ppccmat+reg;
        errp=(sqrt((mone-pp2)./pp2n)+sqrt((mone-ppicmat)./ppicmatn)+sqrt((mone-ppcjmat)./ppcjmatn)+sqrt((mone-ppccmat)./ppccmatn))./sqrt(Meff);
       
        
        
        
        
        for a=1:qa(i)
            for b=1:qa(j)
               % Pij=c(la(i,a),la(j,b))+pw(la(i,a))*pw(la(j,b));
                 
                 fprintf(fid,' %e',max(pp2(a,b),1e-10));
               
               fprintf(fid2,' %e',errp(a,b));
                
            end
        end
        fprintf(fid,'\n');
        fprintf(fid2,'\n');
    end
end
fclose(fid);
fclose(fid2);









% write down the MSA in compact format

    name=[prename,'.cmsa'];

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



% write down the consensus index in the reduced symbols alfabet 

    name=[prename,'.cons'];
    name2=[prename,'.gauge'];
    name3=[prename,'.indgauge'];
    
fid = fopen(name,'w');
fid3 = fopen(name3,'w');
for i=1:N
    fprintf(fid,'%d\n',cons(i));
   fprintf(fid3,'%d\n',indgauge(i));
end
fclose(fid);




% write down the number of regruped a.a. for each site 
    name=[prename,'.qgrouped'];


fid = fopen(name,'w');

for i=1:N
    fprintf(fid,'%d\n',qout(i));
   
end
fclose(fid);

    name=[prename,'.q'];

fid = fopen(name,'w');

for i=1:N
    fprintf(fid,'%d\n',qa(i));
end
fclose(fid);


%write down the frequences for each a.a. in the original alignement (no reweighting no reduction)

name=[prename,'.freq'];
fid = fopen(name,'w');
for i=1:N*q
    fprintf(fid,'%2.6f\n',p(i));
   
end
fclose(fid);

toc


%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
function[qa,qng,ind,qout,indout,cons]=colorreduction(redmethod,redcut,q,p2c,N,indgauge,indmax)
    %ind give the entries of the selected a.a. in the 21*N vector
 %indout=the a.a. grouped in the gauge state 
  %built up ind on the first sites
 
    
 if(strcmp(redmethod,'entropy'));
  % entropy computation
  S_ia=min(-p2c.*log(p2c),p2c~=0);
  S=sum(S_ia);
	
    
  [A,sorted]=sort(S_ia,'descend');
  qnonzero=zeros(N,1);
  cons=zeros(N,1);
    
  for ii=1:1;
    qnonzero(ii)=size(find(A(:,ii)),1);
    Sred=0;
    Sreds=0;
        
    k=0;
    while((k<1)||(Sred/S(ii)<redcut))
      k=k+1;   
      Sreds=Sreds+min(-p2c(sorted(k,ii),ii).*log(p2c(sorted(k,ii),ii)),p2c(sorted(k,ii),ii)~=0);
	  pout=sum(p2c(sorted(k+1:qnonzero(ii),ii),ii));
      Sredr=min(-pout.*log(pout),pout~=0);
      Sred=Sreds+Sredr;
    end
    indt1=sorted(1:k,ii);
    indoutt=sorted(k+1:qnonzero(ii),ii);
    indt=indt1(indt1~=indgauge(ii));
    cons(ii)=sum(find(indt==indmax(ii)));
    indoutt=indoutt(indoutt~=indgauge(ii));
  end

 
  %ind give the entries of the selected a.a. in the 21*N vector
  %indout=the a.a. grouped in the gauge state 
  %built up ind on the first sites
 

 
 
  
 
  qa=zeros(1,N);
  qout=zeros(N,1);
 
  qng=zeros(1,N);
 
  qng(1)=size(indt,1);
  qout(1)=size(indoutt,1);
 
  ind=indt(1:qng(1));
  indout=indoutt(1:qout(1));
 
 
 
  if(qout(1)>0);
     qa(1)=qng(1)+1;
  else
     qa(1)=qng(1);
  end
  %built up ind on other sites
  for ii=2:N;
    
      qnonzero(ii)=size(find(A(:,ii)),1);
      Sred=0;
      Sreds=0;
      k=0;
	  while((k<1)||(Sred/S(ii)<redcut))
        k=k+1;   
        Sreds=Sreds+min(-p2c(sorted(k,ii),ii).*log(p2c(sorted(k,ii),ii)),p2c(sorted(k,ii),ii)~=0);
	    pout=sum(p2c(sorted(k+1:qnonzero(ii),ii),ii));
        Sredr=min(-pout.*log(pout),pout~=0);
        Sred=Sreds+Sredr;
      end
      indt1=sorted(1:k,ii);
      indoutt=sorted(k+1:qnonzero(ii),ii);
      indt=indt1(indt1~=indgauge(ii));
      cons(ii)=sum(find(indt==indmax(ii)));
      indoutt=indoutt(indoutt~=indgauge(ii));
      qout(ii)=size(indoutt,1);
      qng(ii)=size(indt,1);
      indt2=indt(1:qng(ii))+(ii-1)*q;
      if(qng(ii)>0)
          ind=cat(1,ind,indt2);
      end
      indoutt2=indoutt(1:qout(ii))+(ii-1)*q;

      if(qout(ii)>0);
        qa(ii)=qng(ii)+1;
        indout=cat(1,indout,indoutt2);
      else
        qa(ii)=qng(ii);
      end
 
 
  end
 
 
 
 else %frequency reduction
  cons=zeros(1,N);

  indt1=find((p2c(:,1)>redcut));
  indt=indt1(indt1~=indgauge(1));
  cons(1)=sum(find(indt==indmax(1)));
 
 

  indoutt=find((p2c(:,1))&(p2c(:,1)<=redcut));
  indoutt=indoutt(indoutt~=indgauge(1));
 
 
  
 
 
  
  %qa: number of potts states on each site taking into account 
  %grouped symbol (last potts symbol)
  %qng number of non-grouped potts states on each site
  qa=zeros(1,N);
  qout=zeros(1,N);
  qng=zeros(1,N);
  %first site
  qng(1)=size(indt,1);
  qout(1)=size(indoutt,1);
  ind=indt(1:qng(1));
  indout=indoutt(1:qout(1));
 
 
 
  if(qout(1)>0);
     qa(1)=qng(1)+1;
  else
     qa(1)=qng(1);
  end
  %built up ind on other sites
  for ii=2:N
 
 
   indt1=find((p2c(:,ii)>redcut));
   indt=indt1(indt1~=indgauge(ii));
   cons(ii)=sum(find(indt==indmax(ii)));
   qng(ii)=size(indt,1);
   indt2=indt(1:qng(ii))+(ii-1)*q;
   if(qng(ii)>0)
    ind=cat(1,ind,indt2);
   end
   indoutt=find((p2c(:,ii))&(p2c(:,ii)<=redcut));
   indoutt=indoutt(indoutt~=indgauge(ii));
   qout(ii)=size(indoutt,1);
   indoutt2=indoutt(1:qout(ii))+(ii-1)*q;

   if(qout(ii)>0);
     qa(ii)=qng(ii)+1;
     indout=cat(1,indout,indoutt2);
   else
     qa(ii)=qng(ii);
   end
 
  end
 end
end

% function x=letter2number(a)
% switch(a)
% 
%     case '-'
%          x=1;
%     case 'A'    
%         x=2;    
%     case 'C'    
%         x=3;
%     case 'D'
%         x=4;
%     case 'E'  
%         x=5;
%     case 'F'
%         x=6;
%     case 'G'  
%         x=7;
%     case 'H'
%         x=8;
%     case 'I'  
%         x=9;
%     case 'K'
%         x=10;
%     case 'L'  
%         x=11;
%     case 'M'
%         x=12;
%     case 'N'  
%         x=13;
%     case 'P'
%         x=14;
%     case 'Q'
%         x=15;
%     case 'R'
%         x=16;
%     case 'S'  
%         x=17;
%     case 'T'
%         x=18;
%     case 'V'
%         x=19;
%     case 'W'
%         x=20;
%     case 'Y'
%         x=21;
%     otherwise
%         x=1;
% end
% end
% Alphabet used in the lattice Protein Model

function x=letter2number(a)
    switch(a)
        % full AA alphabet
        case 'C'
             x=1;
        case 'M'    
            x=2;    
        case 'F'    
            x=3;
        case 'I'
            x=4;
        case 'L'  
            x=5;
        case 'V'
            x=6;
        case 'W'  
            x=7;
        case 'Y'
            x=8;
        case 'A'  
            x=9;
        case 'G'
            x=10;
        case 'T'  
            x=11;
        case 'S'
            x=12;
        case 'N'  
            x=13;
        case 'Q'
            x=14;
        case 'D'
            x=15;
        case 'E'
            x=16;
        case 'H'  
            x=17;
        case 'R'
            x=18;
        case 'K'
            x=19;
        case 'P'
            x=20;
        
       
    end
end
end
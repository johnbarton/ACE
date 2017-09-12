import pickle
import numpy as np


#####################################################################################
#                                                                                   #
# Example usage from the Python command line (in the main ACE directory)            #
#                                                                                   #
# Import all functions                                                              #
# > from scripts.ACEtools import *                                                  #
#                                                                                   #
# Compute correlations from MSA sequences                                           #
# > WriteCMSA('fasta', 'p7-alignment.fasta', 0, 'entropy', 0.9, 'cons', 1, 'p7')    #
#                                                                                   #
# Read in sequences and tags from MSA                                               #
# > msa, tag = getmsa('examples/p7-alignment.fasta')                                #
#                                                                                   #
# Read in fields (h) and couplings (J)                                              #
# > h, J = getj('examples/p7-out-learn')                                            #
#                                                                                   #
# Read in the map from sequence to Potts state                                      #
# > smap = getsmap('p7')                                                            #
#                                                                                   #
# Compute energy of the first sequence in the MSA                                   #
# > E = getE(h, J, smap, msa[0])                                                    #
#                                                                                   #
# Convert couplings from consensus to zero sum gauge (e.g. for contact prediction)  #
# > J_ZS = consensus2zerosum(h, J)                                                  #
#                                                                                   #
#####################################################################################



# File extension conventions

jext='.j'       # Couplings
pext='.p'       # Standard (one- and two-point) correlations
pnext='.pn'     # Probability of n mutations
paaext='.paa'   # Probability of all AA at each protein site
repext='.rep'   # Nucleotide -> AA sequence report
sceext='.sce'   # Report from a run of the selective cluster expansion



# WriteCMSA port

def WriteCMSA(filetype='fasta', filename='input', theta=0, redmethod='frequency', redcut=0, gauge='cons', gapred=1, out=''):
    """
    A Python version of the WriteCMSA Matlab script. This script calls on a number of auxiliary functions 
    and outputs additional files useful for computing energies, etc. See above for usage examples. Input is
    the same as in the Matlab script.
    
    Input:
    
    - filetype      (string: "fasta", "binary", or "binaryT")
                    If "fasta", process a multiple sequence alignment file in FASTA format. 
                    If "binary", process a binary file, where all entries are zeros or ones, 
                    (e.g. binarized neural recording), where each LINE represents one measurement of the system.
                    If "binaryT", process a "binary" file, but where each COLUMN represents a recording.
    
    - filename      File name of the input alignment
    
    - theta         (real number >= 0, generally from 0 to 0.3)
                    The reweighting threshold, used to take into account correlated sampling
                    (phylogeny). If zero there is no reweighting. Otherwise it weights the
                    contribution of each sequence to the average frequencies and correlations
                    with a weight that is inversely proportional to the number of sequences
                    in the MSA with Hamming distance smaller than theta * N. A typical value
                    for a MSA is 0.2.
                    
    - redmethod     (string: "frequency" or "entropy")
                    If "entropy" the reduction of Potts states is based on the single site entropy
                    contribution, otherwise it is based on the frequency of the Potts state.
                    
    - redcut        (real number >= 0, <=1)
                    Reduction threshold: if 0 only states which are never observed are removed.
                    If reduction is according to frequency, this is the minimum frequency that an AA
                    must have in order to be included explicitly as a Potts state. If reduction is
                    according to entropy, the number of states at each site will be chosen in order
                    to capture at least this fraction of the total single site entropy.

    - gauge         (string: "least", "cons", "group", or a file name),
                    Choice of the gauge state for writing out the correlations.
                    "least": The least frequent (non-grouped) AA at each site.
                    "cons": The consensus AA.
                    "group": The grouped Potts state is set as the gauge state.
                    The gauge can also be manually specified by placing here the name of a file 
                    containing a string with the gauge AA at each site.
    
    - gapred        (0 or 1)
                    When gapred=1, gaps in the alignment are replaced with other AA according
                    to the frequency of AA at the same site in other sequences (without gaps
                    at that site).
                    NOTE: By default this replacement is stochastic, so the resulting set of
                    correlations will not always be identical when this option is used. To use a
                    deterministic replacement (equivalent to the stochastic version in the limit
                    of a large number of replacements), set fillConvention='smooth' in the options below.
    """

    # Set options for passing to getaaq/getbin
    
    reweight    = (theta>0)
    byEntropy   = False
    byFrequency = False
    
    if redmethod=='frequency':
        byEntropy   = False
        byFrequency = True
    elif redmethod=='entropy':
        byEntropy   = True
        byFrequency = False

    if out=='': out = '.'.join(filename.split('.')[:-1])+'-output'
    
    
    if filetype=='fasta':
        getaaq(filein=filename, out=out, saveS=False, removeSingular=True, removeConsGap=False, fillGaps=(gapred==1), fillConvention='smooth', qcut=redcut, byEntropy=byEntropy, byFrequency=byFrequency, gauge=gauge, reweight=reweight, threshold=1.0-theta, pseudocount=0, gaplim=1.0, xlim=1.0, useX=False, savep3=False)

    elif filetype=='binary':
        getbin(filein=filename, out=out, removeSingular=True, reweight=reweight, threshold=1.0-theta, transpose=False, savep3=False)

    elif filetype=='binaryT':
        getbin(filein=filename, out=out, removeSingular=True, reweight=reweight, threshold=1.0-theta, transpose=True, savep3=False)

    else:
        print('Filetype "%s" is not in the current list of options.' % filetype)


def getbin(filein, out='', removeSingular=True, loadWeight=False, reweight=False, threshold=1.0, transpose=False, savep3=False, **kwargs):
    """
    This function reads in sequence data from a binary file, then outputs correlations and supplementary information.
    
    Input:
    - filein            The file for the binary sequence data
    - out               File prefix for writing out correlations, etc
    - removeSingular    If true, remove sites with no variation when writing out correlations
    - loadWeight        Sequence weights can be read in from an input file, specified here
    - reweight          If true, reweight sequences
    - threshold         Similarity fraction for reweighting sequences
    - transpose         If true, transpose the data before processing
    - savep3            If true, write out three-point correlations
    """

    # Read in binary sequences and transpose if necessary

    msa = np.loadtxt(filein)

    if transpose: msa = msa.T
    
    N   = len(msa[0])

    # Sequence reweighting

    B      = float(len(msa))
    Beff   = B
    weight = np.array([1.0 for i in range(len(msa))])
    count  = []

    if loadWeight:
        weight = np.loadtxt(loadWeight, float)
        Beff   = np.sum(weight)
    if reweight: Beff, weight = seqreweight(msa, tag, threshold=threshold)

    # Remove sites with no variation

    p1      = np.sum(weight * msa.T, axis=1) / Beff
    nonsing = (p1>0) * (p1<1)

    if removeSingular:
        msa = msa[:,nonsing]
        N = len(msa[0])

    # Compute correlations

    p12 = np.einsum('i,ij,ik->jk', weight, msa, msa) / Beff

    # Three-point correlations

    if savep3 and out:
        f    = open(out+'.p3', 'w')
        p123 = np.einsum('i,ij,ik,il->jkl', weight, msa, msa, msa) / Beff
        for i in range(N):
            for j in range(i+1,N):
                for k in range(j+1,N):
                    f.write('%.6e\n' % p123[i,j,k])
        f.close()

    # Compute probability of n "mutations" (e.g. active neurons)
    
    nmut = np.sum(msa>0,axis=1)
    pn   = np.zeros(N+1)

    for i in range(len(nmut)):
        padd          = np.zeros(N+1)
        padd[nmut[i]] = weight[i] / Beff
        pn           += padd

    # If there is an output file, then print to file, else return info.
    
    print('Beff = %lf' % Beff)
    
    if out:
        fileout = out
        
        # Write one- and two-point correlations
        
        f = open(fileout+pext,'w')
        
        for i in range(N):
            f.write('%.8e\n' % p12[i,i])
        
        for i in range(N):
            for j in range(i+1,N):
                f.write('%.8e\n' % p12[i,j])

        f.close()

        # Write P(n) "mutations"
        
        fn = open(fileout+pnext,'w')
        for i in pn:
            fn.write('%.8e\t' % i)
        fn.write('\n')
        fn.close()

        # Write report (higher level information)

        consensus  = [str(int(p12[i,i]>=0.5)) for i in range(N)]
        states     = 1 + nonsing

        printReport(fileout, B, Beff, N, consensus, states, nonsing, 1.0-np.array([p12[i,i] for i in range(N)]), 0)
        
        # Write supplementary CMSA files (Matlab equivalent)
        
        printCMSAbin(msa, out=fileout+'.cmsa')      # CMSA

        f = open(fileout+'.cons', 'w')              # Consensus (all zeros b/c of reordering)
        for i in range(N): f.write('%d\n' % 0)
        f.close()

        f = open(fileout+'.wgt', 'w')               # Weight for each sequence
        for i in weight: f.write('%.6e\n' % i)
        f.close()


def getaaq(filein, out='', saveS=False, removeSingular=True, removeConsGap=0.0, fillGaps=False, fillConvention='smooth', qcut=21, byEntropy=False, byFrequency=False, gauge='cons', convert=False, loadWeight=False, reweight=False, threshold=1.0, pseudocount=0, gaplim=1.0, xlim=1.0, useX=False, savep3=False, **kwargs):
    """
    This function reads in sequence data from a .fasta file, then outputs correlations and supplementary information.
    
    Input:
    - filein            The fasta file for the sequence data
    - out               File prefix for writing out correlations, etc
    - saveS             Record single site entropies (T/F)
    - removeSingular    If true, remove sites with no variation when writing out correlations
    - removeConsGap     If true, remove sites with gap frequency > x (float)
    - fillGaps          Replace gaps with random amino acids, two different conventions (see below)
    - fillConvention    'noisy'  : replacement amino acids are chosen at random from the observed distribution
                        'smooth' : gaps are replaced with ambiguous amino acids, which are themselves replaced by a vector mixture representing the single site AA distribution
    - qcut              Cutoff on the number of states at each site, based on number, entropy fraction, or frequency
    - byEntropy         Reduce number of states based on entropy fraction (T/F)
    - byFrequency       Reduce number of states based on frequency
    - gauge             Sets the gauge state at each site, options below
                        'cons'  : the consensus (most frequently observed) AA
                        'least' : the least frequently observed AA
                        'group' : the regrouped state, chosen according to the reduction conventions above
                        else    : an input sequence ('wild-type' in the WriteCMSA Matlab script), the string should point to the sequence file
    - convert           If true, translate an input DNA/RNA sequence into an amino acid sequence
    - loadWeight        Sequence weights can be read in from an input file, specified here
    - reweight          If true, reweight sequences
    - threshold         Similarity fraction for reweighting sequences
    - pseudocount       Adjust output correlations using a pseudocount of x (float, 0 <= x <= 1)
    - gaplim            For quality control, exclude sequences where the total fraction of gaps is > x
    - xlim              For quality control, exclude sequences where the total fraction of ambiguous amino acids is > x
    - useX              Treat amibiguous amino acids ("X") as independent states (T/F)
    - savep3            If true, write out three-point correlations
    """

    # Read in sequences, verify all have the same length
    
    msa, tag = getmsa(filein, convert=convert, noArrow=True)
    
    lengths    = [len(s) for s in msa]
    lmin, lmax = np.min(lengths), np.max(lengths)
    
    assert lmin==lmax, "Error: Sequences are of different lengths. Realign the MSA to ensure equal lengths for all sequences."
    
    # Compute insertion/deletion frequency

    gapfreq = np.sum(msa=='-',1)/float(len(msa))
    
    # Switch case, convert to arrays
    
    for i in range(len(msa)):
        for j in range(len(msa[i])): msa[i][j] = msa[i][j].upper()
    
    msa = np.array(msa)
    tag = np.array(tag)
    N   = len(msa[0])
    
    # Remove sequences with too many ambiguous amino acids or gaps

    nogapseq = [i for i in range(len(msa[0]))]
    msa, tag = cleanAlignment(msa, tag, nogapseq=nogapseq, gapcut=np.ceil(gaplim * len(nogapseq)), xcut=np.ceil(xlim * len(nogapseq)))
    
    # Sequence reweighting

    B      = float(len(msa))
    Beff   = B
    weight = np.array([1.0 for i in range(len(msa))])
    count  = []

    if loadWeight:
        weight = np.loadtxt(loadWeight, float)
        Beff   = np.sum(weight)
    if reweight: Beff, weight = seqreweight(msa, tag, threshold=threshold)

    AA = np.array(['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V','*','-'])
    if useX: AA = np.insert(AA, len(AA), 'X')
    count = np.array([[np.sum((msa[:,i]==aa) * weight) for aa in AA] for i in range(N)], float)

    allowed = [1 for i in range(N)]
    if removeConsGap:
        allowed = count[:,-1]<removeConsGap*Beff
        msa     = np.array(msa[:, allowed])
        tag     = np.array(tag)
        N       = len(msa[0])
        
        # Reevaluate number of gaps/ambiguous amino acids after shortening
        
        B       = float(len(msa))
        Beff    = B
        weight  = np.array([1.0 for i in range(len(msa))])
        if reweight: Beff, weight = seqreweight(msa, tag, threshold=threshold)
        count   = np.array([[np.sum((msa[:,i]==aa) * weight) for aa in AA] for i in range(N)], float)

    # Set number of states, get consensus and maps from AA to Potts configuration
    
    if fillGaps:
        for s in msa:
            for i in range(len(s)):
                if (s[i]=='-'):
                    if fillConvention=='noisy': s[i] = AA[np.random.choice(range(len(AA)),1,count[i]/np.sum(count[i]))[0]]
                    else:                       s[i] = 'X'
        count = np.array([[np.sum((msa[:,i]==aa) * weight) for aa in AA] for i in range(N)], float)

    states, consensus, seqmap, vecmap = getstates(count, qcut, byEntropy=byEntropy, byFrequency=byFrequency, gauge=gauge, useX=useX)
    
    # Remove sites with no variation
    
    nonzero = states>1
    if byFrequency: nonzero = states>2
    if removeSingular:
        nstates = states[nonzero]
        nseqmap = seqmap[nonzero]
        vecmap  = vecmap[nonzero]
        msa     = msa[:,nonzero]
    else:
        nstates = states
        nseqmap = seqmap
    N = len(msa[0])

#    # Save configuration, NOTE: MUST COMMENT OUT REMOVAL OF FULLY CONSERVED SITES TO GET FULL LENGTH BINARY SEQUENCE
#
#    pconfig = np.array([[seqmap[i][seq[i]] for i in range(N)] for seq in msa])
#    #msa1p1s = get1p1s(msa, tag, acc_badlist=acc_blacklist, acc_dupelist=acc_dupelist)
#    p1      = [np.array([np.sum((pconfig[:,i]==q) * weight) for q in range(1,states[i])]) for i in range(N)]
#    nx      = [np.sum((pconfig[:,i]==-1) * weight) for i in range(N)]
#    for i in range(N):  p1[i] /= (Beff - nx[i])
#
#    f = open('binary.dat','w')
#    for s in msa:
#        for i in range(len(s)):
#            binval = seqmap[i][s[i]]
#            if binval==-1:
#                p = np.array([1.-np.sum(p1[i])]+list(p1[i]),float)
#                f.write('%d ' % np.random.choice(range(len(p1[i])+1), p=p))
#            else: f.write('%d ' % binval)
##            if binval==-1:
##                p = np.array(list(p1[i])+[1.-np.sum(p1[i])],float)
##                f.write('%d ' % np.random.choice(range(len(p1[i])+1), p=p))
##            elif binval>0: f.write('%d ' % (binval-1))
##            else:          f.write('%d ' % len(p1[i]))
#        f.write('\n')
#    f.close()
#    f = open('weight.dat', 'w')
#    for i in weight: f.write('%.6e\n' % i)
#    f.close()
#    return 0

    # Convert MSA sequences to Potts configurations + compute one-body correlations, define Potts vector configurations

    pconfig = np.array([[nseqmap[i][seq[i]] for i in range(N)] for seq in msa])
    p1      = [np.array([np.sum((pconfig[:,i]==q) * weight) for q in range(1,nstates[i])]) for i in range(N)]
    nx      = [np.sum((pconfig[:,i]==-1) * weight) for i in range(N)]
    
    for i in range(N):
        p1[i]         /= (Beff - nx[i])
        vecmap[i][-1]  = np.array([j for j in p1[i]])

    # Record entropy (optional)
    
    if saveS:
        f     = open(out+'-S.dat','w')
        count = 0
        for i in range(len(nonzero)):
            if nonzero[i]:
                S      = -(1. - np.sum(p1[count])) * np.log(1. - np.sum(p1[count]))
                S     -= np.sum([p1[count][j] * np.log(p1[count][j]) for j in range(len(p1[count]))])
                count += 1
                f.write('%lf\n' % S)
            else:
                f.write('%lf\n' % 0)
        f.close()

    # Compute two-point correlations

    pcolvector = [np.array([vecmap[i][pconf[i]] for pconf in pconfig]) for i in range(N)]
    p2         = [(np.sum((pcolvector[i][:,:,np.newaxis] * pcolvector[j][:,np.newaxis,:]).T * weight, 2).T).flatten() / Beff for i, j in pairs(N)]

    # Pseudocount

    if pseudocount>0:
        for i in range(N):
            for qi in range(len(p1[i])): p1[i][qi] = ((1. - pseudocount) * p1[i][qi]) + (pseudocount / float(nstates[i]))
        for i, j in pairs(N):
            idx = index(i,j,N)
            for qi in range(len(p1[i])):
                for qj in range(len(p1[j])):
                    sidx          = sindex(qi,qj,len(p1[i]),len(p1[j]))
                    p2[idx][sidx] = ((1. - pseudocount) * p2[idx][sidx]) + (pseudocount / float(nstates[i] * nstates[j]))

    # Three-point correlations

    if savep3 and out:
        f      = open(out+'.p3', 'w')
        idxset = [[i,j,k] for i in range(N) for j in range(i+1,N) for k in range(j+1,N)]
        p3     = [(np.sum((pcolvector[i][:,:,np.newaxis,np.newaxis] * pcolvector[j][:,np.newaxis,:,np.newaxis] * pcolvector[k][:,np.newaxis,np.newaxis,:]).T * weight, -1).T).flatten() / Beff for i, j, k in idxset]
        for i in range(len(p3)):
            if len(p3[i])==0: f.write('%.6e\n' % 0)
            for p in p3[i]:   f.write('%.6e\n' % p)
        f.close()

    # Compute probability of n mutations
    
    nmut = np.sum(pconfig>0,1)
    pn   = np.zeros(N+1)

    for i in range(len(nmut)):
        padd          = np.zeros(N+1)
        padd[nmut[i]] = 1.0

        # Count contribution of imputed amino acids to P(n)
        for j in range(len(pconfig[i])):
            if pconfig[i][j]==-1:
                ptot = np.sum(p1[j])
                padd = ((1. - ptot) * padd) + (ptot * np.array([0] + [padd[k-1] for k in range(1,N+1)],float))
                
        padd *= weight[i] / Beff
        pn   += padd

    # If there is an output file, then print to file, else return info.
    
    print('Beff = %lf' % Beff)
    
    if out:
        fileout = out
        
        # Write one- and two-point correlations
        
        f = open(fileout+pext,'w')
        
        for i in p1:
            if len(i)==0: f.write('%.8e' % 0)
            else:         f.write('%.8e' % i[0])
            for j in i[1:]: f.write('\t%.8e' % j)
            f.write('\n')
        
        for i in p2:
            if len(i)==0: f.write('%.8e' % 0)
            else:         f.write('%.8e' % i[0])
            for j in i[1:]: f.write('\t%.8e' % j)
            f.write('\n')
        f.close()

        # Write P(n) mutations
        
        fn = open(fileout+pnext,'w')
        for i in pn:
            fn.write('%.8e\t' % i)
        fn.write('\n')
        fn.close()

        # Write AA to Potts configuration map
        
        with open(fileout+paaext,'wb') as faa: pickle.dump(seqmap, faa)

        ## Write out accession numbers
        #
        #f = open(fileout+'-accessions.dat', 'w')
        #for i in tag: f.write('%s\n' % i[1:])
        #f.close()

        # Write report (higher level information)

        newnonzero = np.array([k for k in allowed])
        count      = 0
        for i in range(len(allowed)):
            if allowed[i]:
                if not nonzero[count]: newnonzero[i] = 0
                count += 1

        printReport(fileout, B, Beff, N, consensus, states, newnonzero, 1.0-np.array([np.sum(p) for p in p1]), 0)
        
        # Write supplementary CMSA files (Matlab equivalent)
        
        printCMSA(msa, seqmap, out=fileout+'.cmsa') # CMSA

        f = open(fileout+'.cons', 'w')              # Consensus (all zeros b/c of reordering)
        for i in range(N): f.write('%d\n' % 0)
        f.close()

        f = open(fileout+'.wgt', 'w')               # Weight for each sequence
        for i in weight: f.write('%.6e\n' % i)
        f.close()
        
    else: return p1, p2


def cleanAlignment(msa, tag, nogapseq=[], gapcut=False, xcut=False):
    """
    Remove sequences with many gaps/insertions (possible alignment errors) or ambiguous amino acids.
    """

    deltot = 0
    delgap = 0
    delx   = 0
    lenall = len(msa)

    if gapcut:
        #gapcount = np.sum(msa=='-',1)
        gapcount = np.array([np.sum([seq[k]=='-' for k in nogapseq]) for seq in msa])
        gapmean  = np.mean(gapcount)
        selected = np.fabs(gapcount-gapmean)<gapcut
        
        msa     = msa[selected]
        tag     = tag[selected]
        delgap  = lenall-len(msa)
        deltot += delgap
        
    if xcut:
        numx     = np.sum(msa=='X',1)
        selected = numx<xcut
        
        msa     = msa[selected]
        tag     = tag[selected]
        delx    = lenall-delgap-len(msa)
        deltot += delx

    Xcount = 0
    for i in range(len(msa)):
        if 'B' in msa[i] or 'Z' in msa[i] or 'J' in msa[i]:
            for j in range(len(msa[i])):
                if msa[i][j]=='B' or msa[i][j]=='Z' or msa[i][j]=='J':
                    msa[i][j]='X'
                    Xcount  +=1
    
    print('Removed %d of %d sequences (%lf) (%d gap, %d ambiguous), set %d inconclusive amino acids (B, Z, J) to ambiguous (X)' % (deltot,lenall,float(deltot)/float(lenall),delgap,delx,Xcount))

    return msa, tag


def seqreweight(msa, msatag, threshold=1.0):
    """
    Return weighting for a set of sequences based on similarity.
    """

    # Sequence reweighting (similarity)
    
    Beff=1.0
    weight=np.array([1.0 for i in range(len(msa))])

    thresh=int((1.0 - threshold) * len(msa[0]))
    for i in range(len(msa)):
        for j in range(i+1,len(msa)):
            if hamming(msa[i],msa[j])<thresh:
                weight[i]+=1.0
                weight[j]+=1.0
    weight=1.0/weight
    Beff=weight.sum(0)
        
    return Beff, weight


def getstates(count, qcut, byEntropy=False, byFrequency=False, gauge='cons', useX=False, use21=False):
    """
    Determine allowed number of states at each site, and return a map from sequence to Potts configuration.
    """

    AA = np.array(['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V','*','-'])
    
    if useX: AA = np.insert(AA, len(AA), 'X')

    # Determine the maximum allowed states at each site

    seqmap    = []
    vecmap    = []
    states    = []
    consensus = []
    gaugeseq  = []
    
    if gauge not in ['cons', 'least', 'group']: gaugeseq = fopen(gauge)
    
    # Use all 21 possible AA (+ gap) at each site
    
    if use21:
        defsmap = {'X' : -1}
        defvec  = np.zeros(len(AA))
        for i in range(len(AA)): defsmap[AA[i]] = i+1
        defvmap     = {}
        defvmap[0]  = np.array([v for v in defvec])
        defvmap[-1] = np.array([v for v in defvec])
        for i in range(1,len(AA)+1):
            defvec[i-1] = 1.0
            defvmap[i]  = np.array([v for v in defvec])
            defvec[i-1] = 0.0
        for c in count:
            seqmap.append(defsmap)
            vecmap.append(defvmap)
            states.append(len(AA))
            consensus.append('X')

    # Compress the number of states

    else:
        for c in count:
        
            # Treat 'X' as missing data, rather than a specific state (if not in AA)
            tempsmap = {}
            if 'X' not in AA: tempsmap['X'] = -1
            tempvmap = {}
            
            # Sort AAs and counts according to frequency, descending order
            AAsort   = AA[np.argsort(c)][::-1]
            csort    =  c[np.argsort(c)][::-1]
            freq     = csort/np.sum(csort)
            gaugeidx = 0
            nstates  = 0
            
            # REDUCE BY ENTROPY
            if byEntropy:
            
                # Get total entropy
                frnz = freq[freq>0]
                Stot = np.sum([-p * np.log(p) for p in frnz])
                Sfrc = 0.0
                ptot = 0.0

                # Determine number of states based on entropy fraction
                if Stot==0: nstates=1
                else:
                    for i in range(len(frnz)-1):
                        if ptot > 0.0: Sfrc += (1. - ptot) * np.log(1. - ptot)
                        
                        ptot += frnz[i]
                        Sfrc -= frnz[i] * np.log(frnz[i]) + (1. - ptot) * np.log(1. - ptot)
                        
                        if Sfrc > qcut * Stot:
                            nstates=i+2
                            break
                        elif i==len(frnz)-2: nstates=i+2

                # Sanity check
                if nstates<=0: print(freq, nstates)
        
            # REDUCE BY FREQUENCY
            elif byFrequency:
                nstates = 1+np.sum(freq>qcut)
            
            # SIMPLE REDUCE, number of states (min 1, max qcut)
            else:
                nstates = np.min([qcut, np.sum(c>0)])
                nstates = np.max([   1,     nstates])

            # Get gauge state
            if gauge=='cons':                           gaugeidx = 0
            elif gauge=='least' and nstates<len(freq):  gaugeidx = nstates-2
            elif gauge=='least' and nstates>=len(freq): gaugeidx = nstates-1
            elif gauge=='group':                        gaugeidx = nstates-1
            else:                                       gaugeidx = AAsort.index(gaugeseq[i])
            if gaugeidx>=nstates: gaugeidx = nstates-1  # map to grouped state if 'wt' state is grouped
        
            # Map from AA to state
            vec = np.zeros(nstates-1)
            ct  = 1
            
            for i in range(nstates):
                map = ct
                
                # Map gauge state to the zero vector
                if i==gaugeidx:
                    map = 0
                    tempvmap[map] = np.array([v for v in vec])
                
                # Map other states to correct index
                else:
                    vec[ct-1]     = 1.0
                    tempvmap[map] = np.array([v for v in vec])
                    vec[ct-1]     = 0.0
                    ct           += 1
                
                # Map AA to state
                tempsmap[AAsort[i]] = map
                if i==nstates-1:
                    if nstates>1 and gaugeidx!=nstates-1: vec[-1] = 1.0
                    for j in range(nstates, len(AA)):
                        tempsmap[AAsort[j]] = map
                        tempvmap[map]       = np.array([v for v in vec])

            # Append maps to list
            seqmap.append(tempsmap)
            vecmap.append(tempvmap)
            states.append(nstates)
            consensus.append(AAsort[0])

    # Return the result

    return np.array(states), consensus, np.array(seqmap), np.array(vecmap)


def printReport(fileout, B, Beff, N, groundstate, variation, nonzero, concentration, indelfreq):
    """
    Write report (higher level information) from getaa.
    """
        
    fr=open(fileout+repext,'w')
    fr.write('number of samples: %d\t%d\n' % (B,Beff))
    fr.write('final length: %d\n' % N)
    
    fr.write('most probable sequence:\n')
    fr.write(''.join(groundstate))
    #for i in groundstate:
    #    fr.write('%d\t' % i)
    fr.write('\n')
    
    fr.write('observed variation:\n')
    for i in variation:
        fr.write('%d\t' % i)
    fr.write('\n')
    
    fr.write('deleted sites (singular, p=0 or p=1):\n')
    for i in range(len(nonzero)):
        if not nonzero[i]:
            fr.write('%d\t' % i)
    fr.write('\n')
    
    fr.write('observed concentration:\n')
    for i in concentration:
        fr.write('%.4e\t' % i)
    fr.write('\n')
    
    fr.write('average concentration: %.4e\n' % np.mean(concentration))
    fr.write('median concentration: %.4e\n' % np.median(concentration))
    fr.write('insertion/deletion frequency: %.4e\n' % indelfreq)
    
    fr.close()


def printCMSA(msa, smap, out=''):
    """
    Given an input MSA and sequence map, write the corresponding cmsa to 'out'.
    """

    # Get amino acids at each site
    
    AA   = np.array(['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V','*','-'])

    N    = len(msa[0])
    q    = [np.max([smap[i][aa] for aa in AA]) for i in range(N)]
    aaid = [[[] for x in range(q[i])]          for i in range(N)]
    
    for i in range(N):
        for j in range(len(msa)):
            if (msa[j][i] in AA):
                id = smap[i][msa[j][i]]
                if id>0: aaid[i][id-1].append(j)

    # Write out cmsa

    f = open(out,'w')

    for i in range(N):
        for j in range(q[i]):
            f.write('-1\n')
            for k in aaid[i][j]: f.write('%d\n' % k)

    f.close()


def printCMSAbin(msa, out=''):
    """
    Given input binary sequences, write the corresponding cmsa to 'out'.
    """

    # Get amino acids at each site
    
    N    = len(msa[0])
    q    = [1 for i in range(N)]
    aaid = [[[] for x in range(q[i])] for i in range(N)]
    
    for i in range(N):
        for j in range(len(msa)):
            id = int(msa[j][i])
            if id>0: aaid[i][id-1].append(j)

    # Write out cmsa

    f = open(out,'w')

    for i in range(N):
        for j in range(q[i]):
            f.write('-1\n')
            for k in aaid[i][j]: f.write('%d\n' % k)

    f.close()


# Reading in data


def getmsa(ref, convert=False, noArrow=False, noq=True):
    """
    Take an input FASTA file and return the multiple sequence alignment, along with corresponding tags.
    """

    # Read in format for FASTA file
    
    msatemp=[i.split() for i in open(ref).readlines()]
    msatemp=[i for i in msatemp if len(i)>0]
    msa=[]
    tag=[]
    
    for i in msatemp:
        if i[0][0]=='>':
            msa.append('')
            if noArrow: tag.append(i[0][1:])
            else:       tag.append(i[0])
        else: msa[-1]+=i[0]

    msa=np.array([list(i) for i in msa])

    if convert:
        aa=[[codon2aa(i[j]+i[j+1]+i[j+2],noq=noq) for j in range(0,len(msa[0]),3)] for i in msa]
        aa=np.array([list(i) for i in aa])
        return aa, tag

    else: return msa, tag


def getsmap(ref):
    """
    Return a map from AA sequence to Potts configuration from a paa file.
    """
    
    with open(ref+paaext,'rb') as faa: smap=pickle.loads(faa.read())

    return smap


def getj(ref, seqref='', **kwargs):
    """
    Get couplings and fields from .j files.
    """

    return getvec(ref, seqref, ext=jext, **kwargs)


def getp1(ref, seqref=''):
    """
    Return one-point correlations from a pfile. (DEFAULT TO NEW FORMAT)
    """
    
    return getvec1(ref, seqref, ext=pext)


def getp2(ref, seqref=''):
    """
    Return two-point correlations from a pfile. (DEFAULT TO NEW FORMAT)
    """

    v1, v2 = getvec(ref, seqref, ext=pext)

    return v2


def getp1p2(ref, seqref=''):
    """
    Return one-point correlations from a pfile, including deleted sites. (DEFAULT TO NEW FORMAT)
    """

    return getvec(ref, seqref, ext=pext)


def getvec(ref, seqref='', ext='', **kwargs):
    """
    Get couplings/correlations from files.
    """

    # Read in data
    
    dat = fopen(ref+ext)
    
    # Extract one- and two-point variables

    N  = int(s2l(dat))
    v1 = np.array([np.array(i,float) for i in dat[:N]])
    v2 = np.array([np.array(i,float) for i in dat[N:]])
    
    # If there is a sequence reference, use this to determine deleted and vmin
    
    if seqref:
    
        # Get deleted sites and map from full index (including deleted sites) to reduced index
    
        deleted = getdeleted(seqref)
        newN    = N+len(deleted)
        idx     = []
        count   = 0
        
        for i in range(newN):
            if i in deleted: count+=1
            idx.append(i-count)
        
        # Fill in deleted entries
    
        vmin  = 0.0
        newv1 = []
        newv2 = []
        count = 0
        
        if ext=='.j': vmin = np.min([np.min(i) for i in v1])
        qstate = np.max([len(i) for i in v1])>1
        
        if qstate:
            vmin = np.array([vmin])
            for i in range(newN):
                if i in deleted: newv1.append(vmin)
                else:            newv1.append(v1[idx[i]])
            for i, j in pairs(newN):
                if i in deleted and j in deleted: newv2.append(np.zeros(1))
                elif i in deleted:                newv2.append(np.zeros(len(v1[idx[j]])))
                elif j in deleted:                newv2.append(np.zeros(len(v1[idx[i]])))
                else:                             newv2.append(v2[index(idx[i],idx[j],N)])
        
        else:
            vmin = np.array([vmin])
            for i in range(newN):
                if i in deleted: newv1.append(vmin)
                else:            newv1.append(v1[idx[i]])
            for i, j in pairs(newN):
                if i in deleted or j in deleted: newv2.append(np.zeros(1))
                else:                            newv2.append(v2[index(idx[i],idx[j],N)])
            
        v1=np.array(newv1)
        v2=np.array(newv2)

    return v1, v2



# Computations


def getE(h, J, smap, seq, passframeshift=False, passGaps=False, passX=False):
    """
    Get the energy of an input amino acid sequence, using an input sequence map to convert from sequence to Potts configuration.
    """
    
    # Sanity check
    
    assert len(h)==len(seq), 'Length mismatch between sequence and couplings! len(h) = %d, len(seq) = %d' % (len(h),len(seq))
    
    # Assign frameshift/gap/ambiguous codons to wild-type if specified
    
    if passframeshift:
        for i in range(len(smap)): smap[i]['#'] = 0
    
    if passGaps:
        for i in range(len(smap)): smap[i]['-'] = 0

    if passX:
        for i in range(len(smap)): smap[i]['X'] = 0
    
    # Get Potts configuration
    
    N       = len(seq)
    pconfig = [smap[i][seq[i]] for i in range(N)]
    nz      = [i for i in range(N) if pconfig[i]!=0]

    # Compute energy

    E = 0.0

    for i in range(len(nz)):
        E -= h[nz[i]][pconfig[nz[i]]-1]
        for j in range(i+1,len(nz)):
            E -= J[index(nz[i],nz[j],N)][sindex(pconfig[nz[i]]-1,pconfig[nz[j]]-1,len(h[nz[i]]),len(h[nz[j]]))]

    return E



# Accessory functions


def consensus2zerosum(h, J):
    """
    Return a set of couplings in the zero sum gauge, given input in the consensus gauge.
    """
    
    N    = len(h)
    q    = np.array([len(i)          for i in h])
    invq = np.array([1./(float(i)+1) for i in q],float)
    newJ = []

    for i, j in pairs(N):
        idx   = index(i,j,N)
        tempJ = []
        
        # Get values for converting to new gauge
        
        suma = np.zeros(q[i])
        sumb = np.zeros(q[j])
        for a in range(q[i]):
            for b in range(q[j]):
                suma[a] += J[idx][sindex(a,b,q[i],q[j])]
                sumb[b] += J[idx][sindex(a,b,q[i],q[j])]
        sumall = np.sum(suma)

        # Compute zero sum gauge couplings and store in the new vector

        for a in range(q[i]):
            for b in range(q[j]):
                sab = sindex(a,b,q[i],q[j])
                K   = J[idx][sab] - (suma[a] * invq[j]) - (sumb[b] * invq[i]) + (sumall * invq[i] * invq[j])
                tempJ.append(K)
                
            # Add couplings for gauge state (site j)
            K = (sumall * invq[i] * invq[j]) - (suma[a] * invq[j])
            tempJ.append(K)

        # Add couplings for gauge state (site i)
        for b in range(q[j]):
            K = (sumall * invq[i] * invq[j]) - (sumb[b] * invq[i])
            tempJ.append(K)

        # Add double gauge state coupling
        K = sumall * invq[i] * invq[j]
        tempJ.append(K)

        newJ.append(np.array(tempJ))

    # Return new vector J in zero sum gauge

    return newJ


def index(i,j,N):
    """
    Returns the location of {i,j} in the set {{0,1},...,{N-2,N-1}}.
    """
    
    x=[i,j,N]
    x.sort()
    
    i, j, N = x[0], x[1], x[2]
    return int(i * (N - 2) - (i * (i - 1)) / 2 - 1 + j)


def sindex(i,j,qi,qj):
    """
    Returns the location of {i,j} in the set {{0,0},{0,1},...,{0,qj-1},...{qi-1,qj-1}}.
    """
    
    return int((i * qj) + j)


def pairs(N):
    """
    Returns the ordered set of pairs of elements in {0, ..., N-1} for iterating.
    """

    pind=[]

    for i in range(N):
        for j in range(i+1,N): pind.append((i,j))

    return pind


def s2l(x):
    """
    Return the size (N) of a N*(N+1)/2 element vector, given its length.
    """

    return ((np.sqrt(1 + 8 * len(x)) - 1) / 2)


def fopen(ref):
    """
    Open a file and return a list of lines.
    """

    return [i.split() for i in open(ref).readlines()]


def getconsensus(ref):
    """
    Return the consensus amino acid sequence from a report file.
    """
    
    d=fopen(ref+'.rep')[3]
    
    return d[0]


def getdeleted(ref):
    """
    Return the deleted positions (p=0 or p=1) from a report file.
    """

    return [int(i) for i in fopen(ref+'.rep')[7]]


def hamming(seq1,seq2):
    """
    Return the Hamming distance between two sequences (of any kind).
    """

    d=np.sum(np.array(seq1)!=np.array(seq2))

    return d


def codon2aa(c, noq=False):
    """
    Returns the amino acid character corresponding to the input codon.
    """
    
    # If all nucleotides are missing, return gap
    
    if c[0]=='-' and c[1]=='-' and c[2]=='-': return '-'
    
    # Else if some nucleotides are missing, return '?'
    
    elif c[0]=='-' or c[1]=='-' or c[2]=='-':
        if noq: return '-'
        else:   return '?'
    
    # If the first or second nucleotide is ambiguous, AA cannot be determined, return 'X'
    
    elif c[0] in ['W', 'S', 'M', 'K', 'R', 'Y'] or c[1] in ['W', 'S', 'M', 'K', 'R', 'Y']: return 'X'
    
    # Else go to tree
    
    elif c[0]=='T':
        if c[1]=='T':
            if    c[2] in ['T', 'C', 'Y']: return 'F'
            elif  c[2] in ['A', 'G', 'R']: return 'L'
            else:                          return 'X'
        elif c[1]=='C':                    return 'S'
        elif c[1]=='A':
            if    c[2] in ['T', 'C', 'Y']: return 'Y'
            elif  c[2] in ['A', 'G', 'R']: return '*'
            else:                          return 'X'
        elif c[1]=='G':
            if    c[2] in ['T', 'C', 'Y']: return 'C'
            elif  c[2]=='A':               return '*'
            elif  c[2]=='G':               return 'W'
            else:                          return 'X'
        else:                              return 'X'
        
    elif c[0]=='C':
        if   c[1]=='T':                    return 'L'
        elif c[1]=='C':                    return 'P'
        elif c[1]=='A':
            if    c[2] in ['T', 'C', 'Y']: return 'H'
            elif  c[2] in ['A', 'G', 'R']: return 'Q'
            else:                          return 'X'
        elif c[1]=='G':                    return 'R'
        else:                              return 'X'
        
    elif c[0]=='A':
        if c[1]=='T':
            if    c[2] in ['T', 'C', 'Y']: return 'I'
            elif  c[2] in ['A', 'M', 'W']: return 'I'
            elif  c[2]=='G':               return 'M'
            else:                          return 'X'
        elif c[1]=='C':                    return 'T'
        elif c[1]=='A':
            if    c[2] in ['T', 'C', 'Y']: return 'N'
            elif  c[2] in ['A', 'G', 'R']: return 'K'
            else:                          return 'X'
        elif c[1]=='G':
            if    c[2] in ['T', 'C', 'Y']: return 'S'
            elif  c[2] in ['A', 'G', 'R']: return 'R'
            else:                          return 'X'
        else:                              return 'X'
        
    elif c[0]=='G':
        if   c[1]=='T':                    return 'V'
        elif c[1]=='C':                    return 'A'
        elif c[1]=='A':
            if    c[2] in ['T', 'C', 'Y']: return 'D'
            elif  c[2] in ['A', 'G', 'R']: return 'E'
            else:                          return 'X'
        elif c[1]=='G':                    return 'G'
        else:                              return 'X'

    else:                                  return 'X'

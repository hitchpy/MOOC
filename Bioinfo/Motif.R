### Bioinformatic Ch3 Motif finding

neighbors = function(pa, d)
{
  if(d == 0)return(pa)
  if(nchar(pa) == 1)return(c('A','T','C','G'))
  res = NULL
  subs = substr(pa, 2, nchar(pa))
  sn = neighbors(subs, d)
  for(tt in sn)
  {
    ## HammingDistance
    dis = sum(strsplit(tt, split = '')[[1]] != strsplit(subs, split ='')[[1]])
    if(dis < d){res = c(res, paste0('A',tt), paste0('T',tt),paste0('C',tt),paste0('G',tt))}
    else{res = c(res, paste0(substr(pa, 1,1), tt))}
  }
  res
}

motif_enu = function(dna, k, d)
{
  res = NULL
  ll = list()
  for(i in 1: length(dna))
  {
    kmers = sapply(1:(nchar(dna[i])-k+1), function(x){substr(dna[i], x, x+(k-1))})
    ll[[i]] = as.vector(sapply(kmers, neighbors, d))
  }
  pa = unique(unlist(ll))
  for(pp in pa)
  {
    pre = sapply(1:length(dna), function(i){pp %in% ll[[i]]})
    if(all(pre))res = c(res, pp)
  }
  res
}



##### Median string

numtostr = function(num) ## Generate all the possible num patterns
{
  if(num == 1)return(c('A', 'T', 'C', 'G'))
  subs = numtostr(num - 1)
  as.vector(sapply(c('A', 'T', 'C', 'G'), function(h){paste0(h, subs)}))
}

## Find the minimum sum of dis between pattern and all strings in dna
disps = function(pa, dna)
{
  aa = strsplit(pa, split ='')[[1]]
  k = nchar(pa)
  dis = 0
  for(tt in dna)
  {
    kmers = sapply(1:(nchar(tt)-k+1), function(x){substr(tt, x, x+(k-1))})
    bb = strsplit(kmers, split = '')
    hd = min(sapply(bb, function(x)sum( x != aa)))
    dis = dis + hd
  }
  dis
}

medianstring = function(dna, k)
{
  res = NULL
  dis = 100000
  for(kk in numtostr(k))
  {
    temp = disps(kk, dna)
    if(temp < dis){dis = temp; res = kk}
  }
  res
}

## Given a text, find the kmer that has the highest profile
## Should use the log prob and addition instead of multiplication
## For numeric reasons
profile_kmer = function(pa, k, profile)
{
  row.names(profile) = c('A', 'C', 'G', 'T')
  kmers = sapply(1:(nchar(pa)-k+1), function(x){substr(pa, x, x+(k-1))})
  helper <- function(kmer){
    res  = strsplit(kmer, split = '')[[1]]
    prod(sapply(1:k, function(i)profile[res[i], i]))
  }
  scores = sapply(kmers, helper)
  kmers[which.max(scores)]
}

## Greedy motif search algorithm, sequencailly build up 
## Profile matrix and update the best motif.
gmsearch = function(dna, k, t, pseudo = TRUE)
{
  ## Better to work with splited version
  dna2 = strsplit(dna, split = '')
  bestmotif = sapply(dna2, '[', 1:k)
  kmers = sapply(1:(length(dna2[[1]]) - k + 1), function(x)dna2[[1]][x: (x+k -1)])
  for(i in 1: ncol(kmers))
  {
    temp = matrix('', k, t)
    temp[,1] = kmers[, i]
    for(j in 2:t)
    {
      if(pseudo) ## Laplace rule of succession (All add one, a small perturbation)
      {
        profile = apply(temp[, 1:(j -1), drop = FALSE], 1, 
                        function(arow)c(sum(arow == 'A')+1, sum(arow == 'C')+1,sum(arow == 'G')+1,sum(arow == 'T')+1)/(j +3))
      }else{
        profile = apply(temp[, 1:(j -1), drop = FALSE], 1, 
                        function(arow)c(sum(arow == 'A'), sum(arow == 'C'),sum(arow == 'G'),sum(arow == 'T'))/(j - 1))
      }
      temp[, j] = strsplit(profile_kmer(dna[j], k, profile), split = '')[[1]]
    }
    score = function(mat)
    {
      mm = apply(mat, 1, function(arow)names(sort(table(arow), decreasing = TRUE))[1])
      sum(sapply(1:nrow(mat), function(i)sum(mat[i, ] != mm[i])))
    }
    if(score(temp) < score(bestmotif))bestmotif = temp
  }
  apply(t(bestmotif), 1, paste, collapse = '')
}




### Randomized Motif search algorithm, random initialization then update until 
### Score stop increasing 

random_motif_search = function(dna, k, t)
{
  len = nchar(dna[1]) -k + 1
  hd = ceiling(runif(t, 0, len)) # motif's starting position
  bestmotif = sapply(1:t, function(i)substr(dna[i], hd[i], (hd[i]+k -1)))
  score = function(mat) # Helper function to calculate the score
  {
    mm = apply(mat, 2, function(arow)names(sort(table(arow), decreasing = TRUE))[1])
    sum(sapply(1:ncol(mat), function(i)sum(mat[, i] != mm[i])))
  }
  while(TRUE)
  {
    temp = do.call(rbind,strsplit(bestmotif, split = ''))
    profile = apply(temp, 2, 
                    function(acol)c(sum(acol == 'A')+1, sum(acol == 'C')+1,
                                    sum(acol == 'G')+1,sum(acol == 'T')+1)/(t +4))
    motifs = sapply(dna, profile_kmer, k, profile)
    names(motifs) = c(1:t)
    temp2 = do.call(rbind,strsplit(motifs, split = ''))
    if(score(temp2) < score(temp))
    {bestmotif = motifs}else{return(list(score = score(temp2), motif = bestmotif))}
  }
}


##### Gibbs sampling version with a fixed number of iteration

gibbs_motif_search = function(dna, k, t, N)
{
  bs = 0
  len = nchar(dna[1]) -k + 1
  hd = ceiling(runif(t, 0, len))
  bestmotif = sapply(1:t, function(i)substr(dna[i], hd[i], (hd[i]+k -1)))
  score = function(mat) # Helper function to calculate the score
  {
    mm = apply(mat, 2, function(arow)names(sort(table(arow), decreasing = TRUE))[1])
    sum(sapply(1:ncol(mat), function(i)sum(mat[, i] != mm[i])))
  }
  helper = function(profile, text, k, len) # Helper to do random update in each step
  {
    text = strsplit(text, split = '')[[1]]
    res = matrix(0, nrow = len, ncol = k)
    for(i in 1:len){res[i,] = text[i:(i+k -1)]}
    scores = apply(res, 1, function(arow)prod(sapply(1:k, function(x)profile[arow[x], x])))
    res[sample(1:len, 1, prob = scores),]
  }
  idx = ceiling(runif(N, 0, t)) ## All the random rows we need
  for(i in 1:N) # Main Loop
  {
    temp = do.call(rbind,strsplit(bestmotif, split = ''))
    profile = apply(temp[-(idx[i]),], 2, 
                    function(acol)c(sum(acol == 'A')+1, sum(acol == 'C')+1,
                                    sum(acol == 'G')+1,sum(acol == 'T')+1)/(t +4))
    row.names(profile) = c('A', 'C', 'G', 'T')
    temp2 = temp
    temp2[idx[i], ] = helper(profile, dna[idx[i]], k, len)
    if(score(temp2)< score(temp)){
      bestmotif[idx[i]] = paste(temp2[idx[i], ], collapse = '')
      bs = score(temp2)
    }
  }
  list(score = bs, motif = bestmotif)
}









setwd('~/Documents/MOOC//Bioinfo')
condons = read.table('RNA_codon_table_1.txt', 
                     stringsAsFactors = FALSE, header = FALSE, fill = TRUE)

astr = 'ATGGCCATGGCCCCCAGAACTGAGATCAATAGTACCCGTATTAACGGGTGA'
dat = readLines('dataset_96_5.txt')
tran = function(astr)
{
  astr = strsplit(astr, split='')[[1]]
  res = NULL
  for(i in 1:floor(length(astr)/3))
  {
    temp = paste0(astr[((i-1)*3 +1):((i-1)*3 +3)], collapse = '')
    res = paste0(res, condons[condons[,1] == temp, 2])
  }
  res
}

pcode = function(astr, peptite)
{
  res = NULL
  len2 = nchar(peptite)
  astr = strsplit(astr, split='')[[1]]
  len = length(astr)
  astr[astr == 'T'] = 'U'
  rstr = rep('', length(astr))
  rstr[astr =='C'] = 'G'
  rstr[astr =='U'] = 'A'
  rstr[astr =='A'] = 'U'
  rstr[astr =='G'] = 'C'
  rstr = rev(rstr)
  for(i in 1:2)
  {
    if(i == 1) # Original string
    {
      for(j in 1:3) # over 3 starting points
      {
        temp = paste0(astr[j:len], collapse = '')
        aa = gregexpr(peptite, tran(temp))[[1]]
        if(aa[1] == -1){}
        else{
          for(k in 1:length(aa))
          {
            res = c(res, substr(temp, (aa[k]-1)*3 +1, (aa[k]-1)*3 + 3*len2))
          }
        }
      }
    }
    else # compliment string
    {
      for(j in 1:3) # over 3 starting points
      {
        temp = paste0(rstr[j:len], collapse = '')
        aa = gregexpr(peptite, tran(temp))[[1]]
        if(aa[1] == -1){}
        else{
          for(k in 1:length(aa))
          {
            ss = substr(temp, (aa[k]-1)*3 +1, (aa[k]-1)*3 + 3*len2)
            ss = rev(strsplit(ss, split = '')[[1]])
            sf = rep('', length(ss))
            sf[ss == 'A'] = 'U'
            sf[ss == 'U'] = 'A'
            sf[ss == 'G'] = 'C'
            sf[ss == 'C'] = 'G'
            sf = paste0(sf, collapse = '')
            res = c(res, sf)
          }
        }
      }
    }
  }
  gsub('U', 'T',res)
}

mass = read.table('integer_mass_table.txt', stringsAsFactors = FALSE)
### Function to generate the whole mass spectrum for a given peptite
cyclicspectrum = function(peptite)
{
  if(nchar(peptite) == 1)return(c(0, mass[mass[,1] == peptite,2]))
  pp = strsplit(peptite, split ='')[[1]]
  len = length(pp)
  prefixmass = rep(0, len)
  prefixmass[1] = mass[mass[,1] == pp[1],2]
  for(i in 2:length(prefixmass))
  {
    prefixmass[i] = prefixmass[i-1] + mass[mass[,1] == pp[i], 2]
  } ## Compute the prefixmass table
  pmass = prefixmass[len] ## Total weight
  res = NULL
  for(i in 1:(len -1))
    for(j in (i+1):len)
    {
      res = c(res, prefixmass[j] - prefixmass[i])
      if( j<= len) res = c(res, pmass - (prefixmass[j] - prefixmass[i]))
    }
  c(0, sort(res), pmass)
}

spec = scan('dataset_100_5.txt', nlines = 1)
nodup = mass[!duplicated(mass[,2]),]

# Generate the linear spectrum instead of cyclic one.
linearspec = function(peptite)
{
  if(peptite == '')return(0)
  if(nchar(peptite) == 1)return(c(0, mass[mass[,1] == peptite,2]))
  pp = strsplit(peptite, split ='')[[1]]
  len = length(pp)
  prefixmass = rep(0, len)
  prefixmass[1] = mass[mass[,1] == pp[1],2]
  for(i in 2:length(prefixmass))
  {
    prefixmass[i] = prefixmass[i-1] + mass[mass[,1] == pp[i], 2]
  } ## Compute the prefixmass table
  pmass = prefixmass[len] ## Total weight
  res = NULL
  for(i in 1:(len -1))
    for(j in (i+1):len)
    {
      res = c(res, prefixmass[j] - prefixmass[i])
    }
  c(0, sort(c(res, prefixmass)) )
}

consisten = function(peptite, spec) ## Helper function for cps function
{
  if(peptite == '')return(FALSE)
  if(nchar(peptite) == 1)return(if(mass[mass[,1] == peptite,2] %in% spec)TRUE else FALSE)
  lin = linearspec(peptite)
  if(!all(lin %in% spec))return(FALSE)
  t1 = table(lin)
  t2 = table(spec)
  for(na in names(t1))
    if(t1[na] > t2[na])
      return(FALSE)
  TRUE
}

### given a full spectrum, output all possible Cyclic peptides that can generate the spectrum.
cps = function(spec)
{
  pep = list(" " = 0)
  len = length(spec)
  pmass = spec[len]
  lx = 0
  res = NULL
  while(length(pep)>0)
  {
    for(i in 1:length(pep))
    {
      for(j in 1:18)
      {
        pep[paste0(names(pep)[i], nodup[j,1])] = pep[[i]] + nodup[j,2]
      }
    }
    if(lx >0)for(i in 1:lx)pep[[1]] = NULL
    for(nas in names(pep))
    {
      nn = substr(nas, 2, nchar(nas))
      if(pep[[nas]] == pmass)
      {
        if(all(cyclicspectrum(nn) == spec))res = c(res, nn)
        pep[[nas]] = NULL
      }
      else if(!consisten(nn, spec)) pep[[nas]] = NULL
    }
    lx = length(pep)
  }
  res
}

### Just a hack, can merge into cps to process the output.
process_out = function(ss)
{
  res = NULL
  for(s in ss)
  {
    s = strsplit(s, split ='')[[1]]
    temp = as.character(nodup[s[1]== nodup[,1],2])
    for(j in 2:length(s))
    {
      temp = paste0(temp,"-",as.character(nodup[s[j] == nodup[,1], 2]))
    }
    res = c(res, temp)
  }
  res
}

## Given a peptide, get its full spectrum, 
## output the score between experiment spec and full spec of the given pep
cycScore = function(pep, spec)
{
  full = cyclicspectrum(pep)
  temp1 = table(full) ; temp2 = table(spec)
  n1 = names(temp1) ; n2 = names(temp2)
  res = 0
  for(nn in n2) # Loop over the given spec, to count the matches
  {
    if(nn %in% n1) res = res + min(temp2[nn], temp1[nn])
  }
  res
}

linScore = function(pep, spec)
{
  if(pep == '') return(0)
  full = linearspec(pep)
  temp1 = table(full) ; temp2 = table(spec)
  n1 = names(temp1) ; n2 = names(temp2)
  res = 0
  for(nn in n2) # Loop over the given spec, to count the matches
  {
    if(nn %in% n1) res = res + min(temp2[nn], temp1[nn])
  }
  res
}

## Helper function for leadercps, to trim the leaderboard
trim = function(lb, spec, N)
{
  if(length(lb) == 0) return(NULL) ## Final case when all the mass past the parentmass.
  len = length(lb)
  scores = rep(0, len)
  nn = names(lb)
  nn = substr(nn, 2, nchar(nn[1]))
  #for(i in 1:len)
  #{
  #  scores[i] = linScore(substr(nn[i], 2, nchar(nn[i])), spec)
  #}
  scores = sapply(nn, linScore, spec)
  idx = order(scores, decreasing = TRUE)
  lb = lb[idx] ; scores = scores[idx]
  if(N < len){
  pre = scores[N]
  }else{pre = scores[len]}
  lb[scores >= pre]
}

### Cyclic peptide sequencing when the spec is not theoretical but experimental one
### Using leaderboard method with N as threshold
leadercps = function(spec, N)
{
  leaderboard = list(' ' = 0)
  leaderpeptite = list(' ', 0)
  len = length(spec)
  pmass = spec[len]
  while(length(leaderboard) > 0)
  {
    temp = list()
    for(i in 1:length(leaderboard)) ## Expand the leaderboard
    {
      for(j in 1:18)
      {
        temp[paste0(names(leaderboard)[i], nodup[j,1])] = leaderboard[[i]] + nodup[j,2]
      }
    }
    leaderboard = temp
    for(nas in names(leaderboard))
    {
      nn = substr(nas, 2, nchar(nas))
      if(leaderboard[[nas]] == pmass){
        ddd = linScore(nn,spec)
        if(ddd > leaderpeptite[[2]])
          leaderpeptite = list(nas ,ddd)
      }else if(leaderboard[[nas]] > pmass)leaderboard[[nas]] = NULL
    }
    leaderboard = trim(leaderboard, spec, N)
    print(nchar(names(leaderboard)[1]))
  }
  leaderpeptite
}


## Compute the convolution of a given spectrum
convolution = function(spec)
{
  len = length(spec)
  res = NULL
  for(i in 2:len)
  {
    temp = rep(0, i -1)
    for(j in 1:(i-1))
    {
      temp[j] = spec[i] - spec[j]
    }
    res = c(res, temp)
  }
  res[res != 0]
}










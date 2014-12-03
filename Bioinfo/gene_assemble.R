
# Given text, break up to kmers
string_composition = function(k, text)
{
  tmp = strsplit(text, split = '')[[1]]
  len = length(tmp) -k +1
  res = rep('', len)
  for(i in 1:len)
  {
    res[i] = paste(tmp[i: (i+k -1)], collapse = '')
  }
  sort(res)
}

# Given Genome path, construct sequence
genome_path = function(strs)
{
  k = nchar(strs[1])
  len = k + length(strs) -1
  res = rep('', len)
  res[1:k] = strsplit(strs[1], split = '')[[1]]
  j = 2
  for(i in (k+1):len)
  {
    res[i] = substr(strs[j], k, k)
    j = j+1
  }
  paste(res, collapse = '')
}

dat = c('ACCGA',
        'CCGAA',
        'CGAAG',
        'GAAGC',
        'AAGCT')

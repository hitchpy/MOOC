
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


####  Overlap graph, given a collation of k mers, return overlap graph as adjacency list
overlap = function(patterns)
{
  res = NULL
  patterns = sort(patterns)
  k = nchar(patterns[1]) # k mers' K
  prefix = substr(patterns,1, k-1) # Suffix and prefix for matching 
  suffix = substr(patterns, 2, k)
  len = length(patterns) 
  for(i in 1: len)
  {
    if(suffix[i] %in% prefix)
    {
      idx = which(suffix[i] == prefix)
      res = c(res, paste(patterns[i], '->', patterns[idx]))
    }
  }
  res
}

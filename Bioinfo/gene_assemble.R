
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

## DeBruijn Graph, the k mers are the edges and k-1 mers are the nodes
## Output the graph in adjacency list format
DeBruijn = function(k, dat)
{
  # Base on the definition, first construct All unique k-1 mers, then
  # Based on all the edges, add number to the corresponding adjacency matrix.
  edges = string_composition(k, dat)
  nodes = sort(unique(as.vector(sapply(edges, function(ss){c(substr(ss, 1, k-1), substr(ss, 2, k))})))) # Create k-1 mers
  len = length(edges) # For loop over and adding edges
  res = matrix(0, nrow = length(nodes), ncol = length(nodes)) # Adjacent matrix
  rownames(res) = colnames(res) = nodes
  for(i in 1:len)  # Fill in the counts
  {
    pre = substr(edges[i], 1, k -1)
    suf = substr(edges[i], 2, k)
    res[pre, suf] = res[pre, suf] + 1
  }
  fil = apply(res, 1, sum)
  res = res[fil >0, ]
  ## Based on the adjacency matrix, output result
  len = nrow(res)
  temp = rownames(res)
  sapply(1:len, function(i)paste(temp[i], '->', paste(rep(nodes, times = res[i, ]), collapse = ',')))
}













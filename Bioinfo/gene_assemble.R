
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
  nodes = sort(unique(as.vector(sapply(edges, 
                                       function(ss){c(substr(ss, 1, k-1), 
                                                      substr(ss, 2, k))})))) # Create k-1 mers
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


### Similar to the above function, but the starting point is the provided edges(kmers)
### So very similar operation!
DeBruijn_patterns = function(patterns)
{
  k = nchar(patterns[1])
  edges = patterns
  nodes = sort(unique(as.vector(sapply(edges, 
                                       function(ss){c(substr(ss, 1, k-1), 
                                                      substr(ss, 2, k))})))) # Create k-1 mers
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



## Eulerian algorithm for finding the eulerian cycle in a Eulerian directed graph.
## Input elist: a E directed graph as adjacency list
## Output the path that traverse the graph using all the edges exactly once.
dat = c(
'0 -> 3',
'1 -> 0',
'2 -> 1,6',
'3 -> 2',
'4 -> 2',
'5 -> 4',
'6 -> 5,8',
'7 -> 9',
'8 -> 7',
'9 -> 6')
EulerianCycle = function(elist)
{
  ## Convert input to list format, with node name as list name, outlink as elements
  temp = strsplit(elist, split = ' -> ')
  len = length(temp)
  xx = lapply(1:len, function(i)as.numeric(strsplit(temp[[i]][2], split = ',')[[1]]))
  names(xx) = sapply(1:len, function(i)temp[[i]][1])
  ## Construct Eulerian cycle
  res = as.numeric(names(xx)[1]) # for storing the output nodes' order. Init to first node
  while(length(xx) >0 )
  {
    ll = length(res)
    cur_node = as.character(res[ll])
    if(!cur_node %in% names(xx)) # Got stuck in a node with no outlink
    {
      # Rearrange the res, as if start from new node and re-traverse the previous edges.
      position = which(as.character(res) %in% names(xx))[1] # It's guarantee there will be at least one.
      res = c(res[c(position:(ll-1), 1:(position - 1))], res[position])
    }else
    {
      next_node = xx[[cur_node]][1] # Take the first out node in the output link.
      if(length(xx[[cur_node]]) == 1){xx[[cur_node]] = NULL}# Remove the node that use up all the outlinks
      else{xx[[cur_node]] = xx[[cur_node]][-1]} # Else just use the first available outlink
      res = c(res, next_node)
    }
  }
  paste(res, collapse = '->') # output the path
}




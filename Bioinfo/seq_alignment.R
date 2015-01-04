### UCSD Bioinfo MOOC 
### Yu Pei 12-28-2014

## Given the possible coins and the money we need to get the change,
## Output the minimun number of coins for money amount of change.
DPChange = function(money, coins)
{
  len = length(coins)
  res = rep(0, money+1)
  for(i in 2:(money+1)) # i equals the change plus 1, since the first pos is for 0
  {
    res[i] = 9999 ## Arbitrary large number
    for(k in 1:len)
    {
      if(i-1 >= coins[k])
      {
        if(res[i - coins[k]]+1 < res[i])
        {
          res[i] = res[i - coins[k]]+1
        }
      }
    }
  }
  res[money+1]
}


### Dynamic programming  algorithm for finding the length 
### of a longest path in the Manhattan Tourist problem.
### Parameters: n: number of rows, m: number of columns
### Down: weight for the down route as matrix
### Right: edge weight for the horizontal route
### Output: the length of the longest path as a number
ManhattanTourist = function(n, m, Down, Right)
{
  res = matrix(0, nrow = n+1, ncol = m+1)
  ## Initialize the sides of the rectangle
  for(i in 2:(n+1)){res[i, 1] = res[i-1, 1] + Down[i-1, 1]}
  for(i in 2:(m+1)){res[1, i] = res[1, i-1] + Right[1, i-1]}
  
  ## Fill in the remaining nodes in the graph.
  for(i in 2:(n+1))
  {
    for(j in 2:(m+1))
    {
      res[i, j] = max(res[i-1, j] + Down[i-1, j], res[i, j-1] + Right[i, j-1])
    }
  }
  res[n+1, m+1]
}


### More to come
## Longest commom subsequence problem of two sequences, with implicit 
## topological ordering, using backtrack matrix to output the common characters

## Input: two strings v and w
## Output: the backtrack matrix that contain the longest route from source to sink.

LCSBacktrack = function(v, w)
{
  vlen = nchar(v)
  wlen = nchar(w)
  res = matrix(0, nrow = vlen+1, ncol = wlen+1)
  backtrack = matrix(0, nrow = vlen, ncol = wlen) # return value
  for(i in 1:(vlen +1)){res[i, 1] = 0}
  for(i in 1:(wlen +1)){res[1, i] = 0}# initilize the first row and col
  # fill in the remaining
  for(i in 2:(vlen+1))
  {
    for(j in 2:(wlen+1))
    {
      temp = res[i-1, j-1]
      if(substr(v,i-1, i-1) == substr(w, j-1, j-1))
      {
        temp = temp + 1 ## a match
      }
      res[i, j] = max(res[i-1, j] , res[i, j-1], temp)
      ## based on res values, 
      if(res[i,j] == res[i-1, j]){backtrack[i-1, j-1] = 1}
      if(res[i,j] == res[i, j-1]){backtrack[i-1, j-1] = 2} 
      
       ## The sudo code is wrong in this part!(here is the correct version)
      if(substr(v,i-1, i-1) == substr(w, j-1, j-1)){backtrack[i-1, j-1] = 3}
    }
  }
  backtrack
}

OutputLCS = function(backtrack, v, i, j)
{
  if(i ==0 || j ==0){return()}
  if(backtrack[i,j] == 1){OutputLCS(backtrack,v,i-1,j)}
  else if(backtrack[i,j] == 2){OutputLCS(backtrack,v,i,j-1)}
  else{OutputLCS(backtrack,v,i-1,j-1);cat(substr(v,i,i))}
}


## Helper function for LCS for arbitrary DAG
## Ordering a DAG, given an adjacency list, the order is not unique
## Input graph: a graph as adjacency list
## Output: nodes names in topological order

TopologicalOrdering = function(graph)
{
  res = NULL ## Initilize the return ordering
  candidates = NULL 
  for(nn in names(graph))
  {
    ## No inlink
    if(all(!sapply(graph, '%in%', x = nn))){candidates = c(candidates, nn)}
  }
  while(length(candidates)>0)
  {
    cur = candidates[1]
    res = c(res, candidates[1])
    candidates = candidates[-1]
    temp = graph[[cur]]
    graph[[cur]] = NULL
    for(xx in temp)
    {
      if(all(!sapply(graph, '%in%', x = xx))){candidates = c(candidates, xx)}
    }
  }
  res
}

## Longest path in a DAG problem
## Input: raw input from Stepic.org
## s1: source node number, s2: sink node number
## glist: the remaining input graph as linked list 
## Output: The length of the longest path in number and the path as link nodes.

Reachable = function(graph, root)
{
  if(is.null(graph[[root]])){return(root)}
  temp = NULL
  for(x in graph[[root]])
  {
    temp = c(temp, x, Reachable(graph, x))
  }
  temp
}

LPDAG = function(s1, s2, glist)
{
  temp = strsplit(glist, split = ':')
  l1 = list() # Initilize two list, one for the edges
  l2 = list() # one for the weight in each edge
  for(i in 1:length(temp))
  {
    tmp = strsplit(temp[[i]], split = '->')[[1]]
    l1[[tmp[1]]] = c(l1[[tmp[1]]], tmp[2])
    tn = names(l2[[tmp[1]]])
    l2[[tmp[1]]] = c(l2[[tmp[1]]], as.numeric(temp[[i]][2]))
    names(l2[[tmp[1]]]) = c(tn, tmp[2])
  }
  reachable = c(s1, unique(Reachable(l1, s1)))
  l1[!names(l1)%in% reachable] = NULL
  l2[!names(l2)%in% reachable] = NULL
  torder = TopologicalOrdering(l1) ## Guaranteed s1 is the first one?
  res = rep(0, length(torder))
  names(res) = torder
  res2 = rep('', length(torder))
  names(res2) = torder
  for(i in 1:length(torder))
  {
    pre = which(sapply(l1, '%in%', x = names(res)[i]))
    if(length(pre)>0)
    {
      pre = names(l1)[pre]
      res[i] = max(res[pre] + sapply(pre, function(ii)l2[[ii]][names(res)[i]]))
      res2[i] = pre[which.max(res[pre] + sapply(pre, function(ii)l2[[ii]][names(res)[i]]))]
    }
  }
  output = s2
  cur = s2
  while(res2[cur] != s1)
  {
    output = paste0(res2[cur], '->', output)
    cur = res2[cur]
  }
  output = paste0(s1, '->', output)
  list(res[s2], output)
}








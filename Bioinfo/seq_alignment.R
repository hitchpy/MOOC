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














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


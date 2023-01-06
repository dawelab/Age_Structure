### Pseudocode for LIS implementation in alignment chaining
## LIS algorithm
- Input: Array A of integers, length n 
  * In the alignment context, array A stores the starting coordinate of each query alignment in minimap2 output
- Create array L and set L[0] = alignedlen[0] * math.log10(mapqual[0] + 0.001) - abs(query[0]-0)/100
  * In the alignment context, array L stores all scores for each alignment i, of which index is equivalent to that of Array A. 
- Create array prev and set prev[j] = NULL for all j.
  * In the alignment context, array prev stores the index of the best processor based on score calculation
- for i = 1, . . . , n do
- for j = 0, . . . , i-1 do
  * Loop through alignment i and the alignment j (all the predecessors of i based on indexing)
- if A[i] > A[j] and L[i] < L[j] + alignedlen[i] * math.log10(mapqual[i] + 0.001) - abs(query[i]-query[j])/100:
  * In the alignment context, alignment scores are dynamically computed for all the predecessors, 
  the best predecessor is the one with the highest score (if starting coordinate of alignment i is larger that of alignment j)
- set L[i] = L[j]+ alignedlen[i] * math.log10(mapqual[i] + 0.001) - abs(query[i]-query[j])/100
  * The alignment score of i is the accumulative score of i and its predecessors 
- set prev[i] = j.
- end if
- end for
- end for

## Reconstruct-LIS through traceback
- Find the index of largest value of L and set variable x=index(max(L))
  * In the alignment context, we are searching for the end alignment with the biggest accumulative score and its index
- Create array X to store index: X=[x]
- while x != prev[x] and x != 0:
- set x = prev[x]
  * In the alignment context, we are finding the chain of index of all the best predecessors
- X.append(x)
- return X[::-1]
  * Return the chain of alignment indexes in correct orientation 
- end while

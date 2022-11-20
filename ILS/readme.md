## LIS algorithm
- Input: Array A of integers, length n 
  * In the alignment context, array A is the number of query alignments in minimap2 output
- Create array L and set L[1] = alignedlen[0] * math.log10(mapqual[0] + 0.001) - abs(query[0]-0)/100
  * In the alignment context, array L stores all scores for each alignment, of which index is equivalent to that of Array A. 
- Create array prev and set prev[j] = NULL for all j.
  * In the alignment context, array prev stores the index of each alignment
- for i = 1, . . . , n-1 do
- for j = 0, . . . , i-1 do
- if A[i] > A[j] and L[i] < L[j] + alignedlen[i] * math.log10(mapqual[i] + 0.001) - abs(query[i]-query[j])/100:
- set L[i] = L[j]+ alignedlen[i] * math.log10(mapqual[i] + 0.001) - abs(query[i]-query[j])/100
- set prev[i] = j.
- end if
- end for
- end for

## Reconstruct-LIS through traceback
- Find the index of largest value of L and set variable x=index(max(L))
- Create array X to store index: X=[x]
- while x != prev[x] and x != 0:
- set x = prev[x]
- X.append(x)
- return syn_index[::-1]
- end while

def ILS(query,alignedlen,mapqual):
    n = len(query)
    # Declare the list (array) for LIS and
    # initialize LIS values for all indexes
    score = [0]*n
    #The chaining score of anchor i in the first round was calculated with following equation, where len(i) and q(i) are respectively 
    #the length and the mapping quality of anchor i. Gap (i,j) is the distance between anchors i and j, 
    # and x and y are the start and end coordinates.
    score[0] = alignedlen[0] * math.log10(mapqual[0] + 0.001) - abs(query[0]-0)/100
    prev = [0]*n
    for i in range(0, n):
        prev[i] = i

    # Compute optimized LIS values in bottom up manner
    for i in range (1 , n):
        for j in range(0 , i):
            if query[i] > query[j] and score[i] < score[j] + alignedlen[i] * math.log10(mapqual[i] + 0.001) - abs(query[i]-query[j])/100:
                score[i] = score[j]+ alignedlen[i] * math.log10(mapqual[i] + 0.001) - abs(query[i]-query[j])/100
                #aligned segment length * adjusted mapping quality value - adjusted gap penalty
                prev[i] = j
                #print([i,score[i]])
    # traceback
    #find the index of largest value
    idx=np.argmax(score)
    # Create variables to store index
    syn_index=[idx]
    while idx != prev[idx] and idx != 0:
        idx = prev[idx]
        syn_index.append(idx)
    #len(syn_index)
    return syn_index[::-1]

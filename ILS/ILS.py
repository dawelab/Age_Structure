def ILS(A,alignedlen,mapqual):
    n = len(A)
    # Declare the list (array) for LIS and
    # initialize LIS values for all indexes
    L = [0]*n
    #The chaining L of anchor i in the first round was calculated with following equation, where len(i) and q(i) are respectively 
    # the length and the mapping quality of anchor i. Gap (i,j) is the distance between anchors i and j, 
    # and x and y are the start and end coordinates.
    L[0] = alignedlen[0] * math.log10(mapqual[0] + 0.001) - abs(A[0]-0)/100
    
    prev = [0]*n
    for i in range(0, n):
        prev[i] = i

    # Compute optimized LIS values in bottom up manner
    for i in range (1 , n):
        for j in range(0 , i):
            if A[i] > A[j] and L[i] < L[j] + alignedlen[i] * math.log10(mapqual[i] + 0.001) - abs(A[i]-A[j])/100:
                L[i] = L[j]+ alignedlen[i] * math.log10(mapqual[i] + 0.001) - abs(A[i]-A[j])/100
                #aligned segment length * adjusted mapping quality value - adjusted gap penalty
                prev[i] = j
                #print([i,L[i]])
    # traceback
    #find the index of largest value
    x=np.argmax(L)
    # Create variables to store index
    X=[x]
    while x != prev[x] and x != 0:
        x = prev[x]
        X.append(x)
    return X[::-1]

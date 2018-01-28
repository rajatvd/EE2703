def trapForLoop(f,start,ends,h):
    """
    Integrates the function 'f' from 'start' to all the values in 'ends'
    using the Trapezoidal method with step size 'h'. 
    
    Returns a vector of values of the integral for each value in 'ends'.
    
    Assumes 'ends' is sorted in ascending order.
    
    Note: If 'h' is too large when compared to the resolution of values 
    in 'ends', the closest value which is a multiple of 'h' will be used.
    
    (using a for loop)
    """
    
    # convert to numpy array if not
    ends = array(ends)
    
    # initialize variables
    f0 = f(start)
    x=start
    s=0 # cumulative sum
    integrals=[]
    val=0
    
    # reuse cumulative sum for calculating multiple integral values
    for end in ends:
        
        # sum until each end value and store result
        while x<=end:
            val = f(x)
            s += val
            x += h
        
        integrals.append((s-0.5*(f0+val))*h)
        
    return integrals
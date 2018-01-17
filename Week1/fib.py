def fib():
    n=1
    nold=1
    new=0
    for k in range(3,11,1):
        new=n+nold
        nold=n
        n=new
        print k,new
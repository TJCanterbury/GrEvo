
time_func = function(t1 = 10, n=100, r= 100, pals=10){
    print("seconds:")
    t = t1*r*(n^2 + n)
    t = t / pals
    print(t)
    print("minutes:")
    print(t / 60)
    print("hours:")
    print(t / (60*60))
}
def KodurFrr(fc,D,L,C,K,is_circular=True,aggregate=True):
    if is_circular:
        if aggregate:f=0.07
        else:f=0.08
    else:
        if aggregate:f=0.06
        else:f=0.07
    R=f*(fc+20)/(K*L-1000)*D**2*(D/C)**0.5
    return R

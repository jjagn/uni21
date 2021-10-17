CF = 5000
i = 10.3 # IN %!!!
n = 8

FV_Annuity = CF/(i/100)*((1+i/100)**n-1)
FV_Annuity_Due = FV_Annuity * (1+i/100) + CF

print(FV_Annuity_Due)
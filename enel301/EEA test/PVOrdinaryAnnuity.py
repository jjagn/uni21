CF = 50000
i = 3 # IN %!
n = 9

PV_Annuity = CF / (i/100) * (1- 1/ (1 + i/100) ** n)
print(PV_Annuity)



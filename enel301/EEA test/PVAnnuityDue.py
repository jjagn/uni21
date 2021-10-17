CF = 50000
i = 3 # IN %!
n = 10

PV_Annuity_Due = (CF / (i/100) * (1- 1/ (1 + i/100) ** n)) * (1+i/100)

print(PV_Annuity_Due)



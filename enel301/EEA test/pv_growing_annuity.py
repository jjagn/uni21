CF = 675000  # cash flow p.a.
i = 18       # interest rate IN %!!!
g = 13     # rate of growth of annuity in %!!!
n = 15     # number of years annuity lasts for

PV = CF/(i/100 - g/100) * (1 - ((1 + g/100)/(1 + i/100)) ** n)

print(PV)
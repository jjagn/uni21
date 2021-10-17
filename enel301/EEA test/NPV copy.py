cashflows = [-93030, -18667, -18667, -18667, -18667, -18667, -18667, (-18667+5499)]
cashflows_enumerated = enumerate(cashflows)

r = 12                                  # discount rate IN %!
NPV = 0

n = len(cashflows) - 1

annuity_factor = (1 - 1/(1 + r/100) ** n) / (r/100)

for year, yearly_cashflow in cashflows_enumerated:
    NPV += yearly_cashflow/(1+r/100) ** year

    print("Year: {}, Cashflow: {}, Rolling NPV: {}".format(year, yearly_cashflow, NPV))

print("Final NPV: {}".format(NPV))

eac = NPV / annuity_factor

print("EAC: {}".format(eac))


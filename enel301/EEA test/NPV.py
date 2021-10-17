cashflows = [-12126, 1315.15, 1315.15, 1315.15, 1315.15, 1315.15, 1315.15, 1315.15, 1315.15, 1315.15, 1315.15, 1315.15, 1315.15, 1315.15, 1315.15, 1315.15, 1315.15, 1315.15, 1315.15, 1315.15, 1315.15]
cashflows_enumerated = enumerate(cashflows)

r = 4                                  # discount rate IN %!
NPV = 0

n = len(cashflows) - 1

annuity_factor = (1 - 1/(1 + r/100) ** n) / (r/100)

for year, yearly_cashflow in cashflows_enumerated:
    yearly_cashflow_adj_inflation = (yearly_cashflow*(1+0.0196)**year)
    # yearly_cashflow_adj_inflation = yearly_cashflow
    pv_cashflow = yearly_cashflow_adj_inflation/(1+r/100) ** year

    NPV += pv_cashflow

    print("Year: {}, Cashflow: {}, Cashflow adjusted for inflation: {}, PV of cashflow: {}, Rolling NPV: {}".format(year, yearly_cashflow,yearly_cashflow_adj_inflation, pv_cashflow, NPV))

print("Final NPV: {}".format(NPV))

eac = NPV / annuity_factor

print("EAC: {}".format(eac))


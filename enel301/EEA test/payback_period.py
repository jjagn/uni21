cashflows = [-12126, 1315.15, 1315.15, 1315.15, 1315.15, 1315.15, 1315.15, 1315.15, 1315.15, 1315.15, 1315.15, 1315.15, 1315.15, 1315.15, 1315.15]
cashflows_enumerated = enumerate(cashflows)

cumulative_cashflow = 0

for year, cashflow in cashflows_enumerated:
    cumulative_cashflow += cashflow
    if cumulative_cashflow >= 0:
        # print(year, cashflow, cumulative_cashflow)
        PB = year - 1 - ((cumulative_cashflow-cashflow)/cashflow)
        break

print(PB) 
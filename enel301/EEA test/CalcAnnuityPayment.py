PV = 23500  # current value of loan
i = 8.4     # IN %!!!
n = 7       # no. of payments

CF = (PV * i/100) / ((1 - (1 / (1 + i/100) ** n)) * (1 + i/100))

print(CF)

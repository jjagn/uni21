# FUTURE VALUE WITH MULTIPLE CASH FLOWS

# input values
values = [13227, 15611, 18970, 19114]
print(values)
interest_rate = 8 #IN %!

print("cashflows: {}".format(values))
print("interest rate = {}%\n".format(interest_rate))

# loop
total_future_value = 0
i = 1
for value in values:
    print("year {}".format(i))
    current_value = value * (1 + interest_rate / 100)**(len(values)-i)
    print("cashflow = {}".format(value))
    print("adjusted for interest = {}".format(current_value))
    print("\n")
    total_future_value += current_value
    i += 1

print(total_future_value)

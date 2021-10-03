# PRESENT VALUE WITH MULTIPLE CASH FLOWS

# input values
values = [1225, 1350, 1500, 1600, 1600]
discount_rate = 8 #IN %!

print("cashflows: {}".format(values))
print("discount rate = {}%\n".format(discount_rate))

# loop
present_value = 0
i = 1
for value in values:
    print("year: {}".format(i))
    current_value = value/(1 + discount_rate / 100)**i
    present_value += current_value
    i += 1
    print("cashflow = {}".format(value))
    print("adjusted for interest = {}\n".format(current_value))

print(present_value)

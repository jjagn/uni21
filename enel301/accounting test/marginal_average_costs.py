variable_cost = [20.500, 10.500, 8.000]
fixed_cost = [100000, 350000, 600000]
units = 100000
s_p = 85

for i in range(0, len(variable_cost)):
    v_c = variable_cost[i]
    f_c = fixed_cost[i]
    print(i+1)
    i = i + 1
    marginal_cost = (f_c + v_c * (units + 1)) - (f_c + v_c * units) 
    print("marginal cost = ${}".format(marginal_cost))

    average_cost = (f_c/units + v_c)
    print("average cost = ${}".format(average_cost))

    breakeven_units = f_c / (s_p-v_c) 
    print("breakeven quantity = {}".format(breakeven_units))

    print("\n")
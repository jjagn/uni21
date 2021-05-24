import math

# required revolutions in 20 year life
required_life = 506.88

f0 = 13.6 # calculation factor from bearing specs
C = 229 # basic dynamic load rating C
C0 = 216 # basic static load rating C0
e = 0.29 # load ratio factor


X = 0.56 # factor for axial load
Y = 1.5 # factor for radial load

c = 0.2 # distance of pulley from bearing

Fa = 18 # axial force in kN
Fr = math.sqrt(4100.45 ** 2 + (7411.24 + 1543.32 * c) ** 2) / 1000 # radial force in kN

calc_factor = f0 * Fa / C0

print(calc_factor)

# print(Fr)

P = X*Fr + Y*Fa # equivalent dynamic load

print("C must be > {:.2f}".format(8*P))

rated_life = (C/P) ** 3

if rated_life > required_life:
    print("bearing passes, rated for: {:.2f} million cycles".format(rated_life))
else:
    print("bearing fails life expectancy, only rated for: {:.2f} million cycles of 506.88".format(rated_life))


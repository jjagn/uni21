import math

# for SKF bearing designation 22324 CCJA/W33VA406

# required revolutions in 20 year life
required_life = 506.88

c = 0.2 # distance of pulley from bearing

# loads
Fa = 18 # axial force in kN
Fr = math.sqrt(4100.45 ** 2 + (7411.24 + 1543.32 * c) ** 2) / 1000 # radial force in kN

# parameters from bearing datasheet
C = 652 # basic dynamic load rating C
C0 = 675 # basic static load rating C0
B =  58 # bearing width [mm]
d =  120 # bearing bore diameter [mm]

# from calc factor and tables
e = 0.26 # load ratio factor
Y0 = 2.5
Y1 = 2.6
Y2 = 3.9

axial_load_allowed = 0.003 * B * d
print("allowed axial load = {:.2f}kN".format(axial_load_allowed))

P0 = Fr + Y0 * Fa

if Fa/Fr <= e:
    P = Fr + Y1 * Fa
else:
    P = 0.67 * Fr + Y2 * Fa # equivalent dynamic load

print("P0 = {:.2f}kN".format(P0))
print("P = {:.2f}kN".format(P))

rated_life = (C/P) ** (10/3)

if rated_life > required_life:
    print("bearing passes, rated for: {:.2f} million cycles".format(rated_life))
else:
    print("bearing fails life expectancy, only rated for: {:.2f} million cycles of 506.88".format(rated_life))

print('C/P = {:.2f}'.format(C/P))


# Converting dadi units to physical units
# constants
mu = 1.86*10^(-8)
gen.small = 2
gen.big = 5
L.genom = 390206788

# dadi params
theta = 6183
T = 3.68
T2 = 1.91

# effective sequence params
unfil.snps = 33523695
sfs.snps = 323336.46


# Step 1. Calculate effective L
L.effective = L.genom * (sfs.snps/unfil.snps)

# Step 2. Calculate Nref (theta = 4Nref*mu*L)
nref = theta / (4*mu*L.effective)
cat("Nref is", nref, "\n")

# Divergence time
Tgen = 2 * nref * T
Tyears.1 = Tgen * gen.small
Tyears.2 = Tgen * gen.big
cat("Divergence time between", Tyears.1, "and", Tyears.2, "years\n")

# T2
System = <NAME>

# Job control
LCHARG = .TRUE.
LAECHG = .TRUE.
NPAR   = 4

# Dipole correction
LDIPOL = .TRUE.
IDIPOL = 3
DIPOL  = 3.94128 2.2755 7.19335

# Electronic minimization
PREC     = Accurate
GGA      = BF 
LUSE_VDW = .TRUE. 
Zab_VDW  = -1.8867
LBEEFENS = .TRUE.
ENCUT    = 500
ALGO     = VeryFast
VOSKOWN  = 1
ISPIN    = 2

# DOS related values
ISMEAR = 0
SIGMA  = 0.05

# Ionic relaxation
IBRION = 2
EDIFF  = 1E-5
EDIFFG = -0.01
NELM   = 250
NSW    = 0 

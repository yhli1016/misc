include(defs.m4)dnl
System = NAME

# Job control
LCHARG = .FALSE.
NPAR   = 4

# Dipole correction
LDIPOL = .TRUE.
IDIPOL = 3
DIPOL  = 3.94128 2.2755 7.19335

# Electronic minimization
PREC     = Normal 
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

# NEB
ICHAIN = 0
IMAGES = NIMAGE
SPRING = -5.0
LCLIMB = .TRUE.
IBRION = 3
POTIM  = 0
IOPT   = 3
EDIFF  = 1E-4
EDIFFG = -0.05
NELM   = 250
NSW    = 100

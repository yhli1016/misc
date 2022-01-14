System = <NAME>

# Job control
ISTART = 1        (read WAVECAR if it exists)
ICHARG = 0        (charge density from wave function)
LWAVE  = .TRUE.   (save wave functions for next run)
LCHARG = .FALSE.  (do not save charge density)
NPAR   = 4

# System
ISPIN = 2

# Electronic minimization
PREC     = Normal 
ENCUT    = 500
VOSKOWN  = 1
EDIFF    = 1E-4   (may reduce to 1E-5 for final run)
NELMIN   = 8      (opt. requires better forces)
NELM     = 60     (change ALGO if convergence is difficult)
ALGO     = Fast   (Davidson+RMM-DIIS)

# Ionic relaxation
EDIFFG = -0.05    (reduce abs. to -0.02 for final run)
NSW    = 100      (rely on your insights for convergence)
IBRION = 3        (disable VASP optimizers)
POTIM  = 0        (disable VASP optimizers)

# NEB
ICHAIN = 0
IMAGES = <NIMAGE>
SPRING = -5.0    (increase abs. to -10.0 for strong bonds)
LCLIMB = .FALSE. (disable climbing imagefor initial iterations)
IOPT   = 7       (1-quasi_newton/2-cg/3-quick_min/7-fire)

# DOS related values
ISMEAR = 0
SIGMA  = 0.05

# Dipole correction
LDIPOL = .TRUE.
IDIPOL = 3
DIPOL  = 3.94128 2.2755 7.19335

# vdW correction
GGA      = BF 
LUSE_VDW = .TRUE. 
Zab_VDW  = -1.8867
LBEEFENS = .TRUE.

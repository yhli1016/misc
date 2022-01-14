System = <NAME>

# Job control
ISTART = 1        (read WAVECAR if it exists)
ICHARG = 0        (charge density from wave function)
LWAVE  = .TRUE.   (save wave functions for next run)
LCHARG = .TRUE.   (save charge density)
LAECHG = .TRUE.   (save core-electron charge)
NPAR   = 4

# System
ISPIN = 2

# Electronic minimization
PREC     = Normal 
ENCUT    = 500
VOSKOWN  = 1
EDIFF    = 1E-4
NELMIN   = 2      (no need for larger value for scf)
NELM     = 60     (change ALGO if convergence is difficult)
ALGO     = Fast   (Davidson+RMM-DIIS)

# Ionic relaxation
NSW    = 0      (fix all atoms)
IBRION = -1     (fix all atoms)

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

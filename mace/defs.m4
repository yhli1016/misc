define(`PREFIX', `test')dnl
define(`NTYP', `1')dnl
define(`NAT', `4')dnl
define(`ECUTWFC', `60.5')dnl
define(`ECUTRHO', `242.0')dnl
define(`TEST', `# TEST')dnl
define(`TEST2', `# TEST2')dnl
define(`DP1', `
    tefield          = .TRUE.
    dipfield         = .TRUE.
')dnl
define(`DP2', `
    edir             = 3
    emaxpos          = 0.90
    eopreg           = 0.05
    eamp             = 0.0
')dnl
define(`POS', `
ATOMIC_SPECIES
   P  30.97   P.pbe-mt_fhi.UPF
CELL_PARAMETERS {alat}
   3.343472076   0.000000000   0.000000000
   0.000000000   4.675443939   0.000000000
   0.000000000   0.000000000  20.000000000
ATOMIC_POSITIONS {crystal}
   P   0.800035018   0.092073099   0.567952767
   P   0.800156942   0.258150821   0.460485652
   P   0.300083492   0.575667132   0.474444252
   P   0.299987230   0.775400398   0.576765798
')dnl
define(`KPT', `K_POINTS {automatic}
 8 6 1 0 0 0
')dnl

changequote([,])

# CONTROL
define([PREFIX], [''])

# SYSTEM
define([NTYP], )
define([NAT], )
define([ECUTWFC], )

define([NBND], )
define([NBNDQ], )
define([NBND_FI], )
define([NBND_PATH], )

# ATOMIC_POSITIONS
define([POS], [
])


# K_POINTS
define([KPT_SCF], [
K_POINTS {automatic}
])

define([KPT], include(wfn.out))
define([KPTQ], include(wfnq.out))
define([KPT_FI], include(wfn_fi.out))

define([KPT_PATH], [
K_POINTS {crystal_b}
])

# 02-wfn
define([NK1], )
define([NK2], )
define([NK3], )
define([DK1], )
define([DK2], )
define([DK3], )
define([XC_MIN], )
define([XC_MAX], )

# 03-wfnq
define([DK1Q], )
define([DK2Q], )
define([DK3Q], )

# 05-wfn_fi
define([NK1_FI], )
define([NK2_FI], )
define([NK3_FI], )

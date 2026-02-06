### constants.jl
#
# Authors: Sofia Calgaro, Toby Dixon
# 
###

"""
    Constants

Module defining physical and experimental constants used throughout ZeroNuFit.jl.

# Constants
- `Qbb`: Q-value for ⁷⁶Ge double-beta decay (2039.06 keV)
- `N_A`: Avogadro's number (6.022E23 mol⁻¹)
- `me_keV`: Electron mass in keV (510.99895 keV, PDG value)
- `m_76`: Molar mass of ⁷⁶Ge (75.92E-3 kg/mol)
- `sig_units`: Signal units scaling factor (1e-27)
- `gamma_2113_keV`: Gamma line energy at 2113.513 keV
- `gamma_2123_keV`: Gamma line energy at 2123.513 keV
- `gamma_2098_keV`: Gamma line energy at 2098.511 keV
- `gamma_2108_keV`: Gamma line energy at 2108.511 keV
- `gamma_2199_keV`: Gamma line energy at 2199.100 keV
- `gamma_2209_keV`: Gamma line energy at 2209.100 keV
"""
module Constants

# define global constants
const Qbb = 2039.06 # keV
const N_A = 6.022E23
const me_keV = 510.99895000 # PDG
const m_76 = 75.92E-3 # kg/mol
const sig_units = 1e-27 # signal is in units of this
const gamma_2113_keV = 2113.513
const gamma_2123_keV = 2123.513
const gamma_2098_keV = 2098.511
const gamma_2108_keV = 2108.511
const gamma_2199_keV = 2199.100
const gamma_2209_keV = 2209.100

end

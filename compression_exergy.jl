# pressure drop in a gas pipe line
# assumptions:
# the velocity range is 60 to 80 ft/s
# for CO2, the maxiumum velocity is 50 ft/s
# minimum gas velocity is 10 to 15 ft/s
# source: http://petrowiki.org/Pipeline_design_consideration_and_standards#Gas_line_sizing
include("CoolProp.jl")

# gas flow rate must be specified
Q_g= 100 # m^3/s at injection condition
gas_type="CO2"
T_pipe=25+273.15 # K pipe line temperature
p_inj=
if gas_type=="CO2"
  v_g=50*0.3048 # ft/s to m/s
elseif gas_type="N2"
  v_g=80*0.3048 # ft/s to m/s
end

rho_gas=CoolProp.PropsSI("DMOLAR", "T", temperature_in_K,
  "P", pressure_in_Pa, gas_type) # gas density mol/m3
s_gas=CoolProp.PropsSI("SMOLAR", "T", temperature_in_K,
  "P", pressure_in_Pa, gas_type) # molar density J/(mol.K)
h_gas=CoolProp.PropsSI("HMOLAR", "T", temperature_in_K,
  "P", pressure_in_Pa, gas_type) # gas enthaply J/mol
mu_gas=CoolProp.PropsSI("VISCOSITY", "T", temperature_in_K,
  "P", pressure_in_Pa, gas_type) # gas viscosity Pa.s
Z_gas=CoolProp.PropsSI("Z", "T", temperature_in_K,
  "P", pressure_in_Pa, gas_type) # gas compressibility factor

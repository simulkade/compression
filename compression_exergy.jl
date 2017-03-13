# pressure drop in a gas pipe line
# assumptions:
# the velocity range is 60 to 80 ft/s
# for CO2, the maxiumum velocity is 50 ft/s
# minimum gas velocity is 10 to 15 ft/s
# %% source: http://petrowiki.org/Pipeline_design_consideration_and_standards#Gas_line_sizing
using Plots, Roots, CoolProp, DataFrames
# include("CoolProp.jl")

# %% constants
R=8.314 # J/(mol.K) gas constants
p_ref=101325 # 1atm in Pa
T_ref=21.0+273.15 # reference temperature
ex_steel_pipe=60e6 # J/kg steel Szargut
rho_steel=8000 # kg/m3 average value (a bit high)
t_life_pipe=30*365*24*3600 # [s] project life time
results_out=Dict("flow rate"=>zeros(0), "pipe diamter"=>zeros(0),
            "transport comp stages"=>zeros(0), "transport exergy"=>zeros(0),
            "steel exergy"=>zeros(0), "injection exergy"=>zeros(0),
            "preessure drop bar"=>zeros(0), "pipe length"=>zeros(0),
            "pipe thickness"=>zeros(0), "transport exergy th"=>zeros(0),
            "injection exergy th"=>zeros(0))
# %% unit conversion
m_to_mile=0.000621371
K_to_degR=1.8
pa_to_psi=0.000145038
m3_to_ft3=35.3147
# sm3s_to_mmscfd=m3_to_ft3/1e6*24.0*3600.0
sm3s_to_mmscfd= 3.05119008
m_to_inch=39.3701

# %% reservoir condition
p_res=280e5 # Pa
T_res= 90+273.15 # K
# gas_type="CO2"
gas_type="N2"
# a mixture can be defined as:
#gas_type="CO2[0.7]&N2[0.3]"

# pipe line specifications
Q_g= 109260/(24*3600) # m^3/s at reservoir condition
L_pipe=200000.0 # m
n_segments=200
T_pipe=35+273.15 # K pipe line temperature
p_end=80e5 # pressure at the end of pipe line
# p_sat_gas=CoolProp.PropsSI("P","T",T_pipe,"Q",0.0,gas_type)
C_gas_res=CoolProp.PropsSI("DMOLAR", "T", T_res,
  "P", p_res, gas_type) # gas density mol/m3
Q_g_molar=Q_g*C_gas_res # mol/s molar flow rate

# calculate the pipe size based on the maximum velocity
# at the final segment of the pipe
if gas_type=="CO2"
  v_g=50*0.3048 # ft/s to m/s
elseif gas_type=="N2" || gas_type=="CH4"
  v_g=80*0.3048 # ft/s to m/s
else
  v_g=80*0.3048 # ft/s to m/s
end


C_gas_end=CoolProp.PropsSI("DMOLAR", "T", T_pipe,
  "P", p_end, gas_type) # gas density mol/m3
Q_end=Q_g_molar/C_gas_end # m^3/s flow at the end segment
d=(4*Q_end/(π*v_g))^0.5 # m pipe diameter
d_inch=d*m_to_inch

# p_sat_gas=CoolProp.PropsSI("P","T",T_pipe,"Q",0.0,gas_type)
C_gas_res=CoolProp.PropsSI("DMOLAR", "T", T_res,
  "P", p_res, gas_type) # gas density mol/m3
Q_g_molar=Q_g*C_gas_res # mol/s molar flow rate

C_gas_end=CoolProp.PropsSI("DMOLAR", "T", T_pipe,
  "P", p_end, gas_type) # gas density mol/m3
Q_end=Q_g_molar/C_gas_end # m^3/s flow at the end segment
C_gas_std=CoolProp.PropsSI("DMOLAR", "T", T_ref,
  "P", p_ref, gas_type) # gas density mol/m3
Q_g_std=Q_g_molar/C_gas_std
d=(4*Q_end/(π*v_g))^0.5 # m pipe diameter
d_inch=d*m_to_inch
e_pipe=0.9 # pipe efficiency
# p_sat_gas=CoolProp.PropsSI("P","T",T_pipe,"Q",0.0,gas_type)
C_gas_res=CoolProp.PropsSI("DMOLAR", "T", T_res,
  "P", p_res, gas_type) # gas density mol/m3
Q_g_molar=Q_g*C_gas_res # mol/s molar flow rate

C_gas_end=CoolProp.PropsSI("DMOLAR", "T", T_pipe,
  "P", p_end, gas_type) # gas density mol/m3
Q_end=Q_g_molar/C_gas_end # m^3/s flow at the end segment
d=(4*Q_end/(π*v_g))^0.5 # m pipe diameter
d_inch=d*m_to_inch
d_inch_min=d_inch-d_inch%2 # make it a standard pipe diameter

# start the calculations here
# workflow:
# choose a few pipe diameters
# design the pipe thickness: X60 grade,
# yield strength: Sy=413 MPa (60000 psi)
# t=P [psi]* d_pipe (inch, outside)/(2*F*E*Sy)
# F=0.72 (average), E=0.6-1.0 weld joint factor
# standard nominal sizes: 10:2:36 (most of the time equal to the OD)
# calculate the pressure drop and pipe material exergy
Sy=413e6 # Pa yield strength
e_weld=0.8
F_pipe=0.72
for j in 1:5
  d_inch+=2.0 # standard pipe sizes
  p=zeros(n_segments+1)
  p[end]=p_end
  for i in n_segments:-1:1
    sg_gas=CoolProp.PropsSI("DMASS", "T", T_ref,"P", p_ref, gas_type)/
      CoolProp.PropsSI("DMASS", "T", T_ref,"P", p_ref, "Air")
    Z_gas=CoolProp.PropsSI("Z", "T", T_pipe,"P", p[i+1], gas_type) # gas compressibility factor
    C_gas=CoolProp.PropsSI("DMOLAR", "T", T_pipe, "P", p[i+1], gas_type) # gas density mol/m3
    p[i]=sqrt((p[i+1]*pa_to_psi)^2+sg_gas^0.961*Z_gas*L_pipe/n_segments*m_to_mile*T_pipe*K_to_degR
    *(Q_g_std*sm3s_to_mmscfd/(0.028*e_pipe*d_inch^2.53))^(1.0/0.51))/pa_to_psi
  end

  # plot(linspace(0, L_pipe, n_segments+1)/1000, p/1e5,
  #      xlabel="Pipe length [km]", ylabel="P [bar]")

  t_pipe=p[1]*d_inch/m_to_inch/(2*F_pipe*e_weld*Sy) # m pipe thickness
  println("the pipe thickness is $t_pipe m")
  v_steel=π/(4*m_to_inch^2)*((d_inch+t_pipe)^2-d_inch^2)*L_pipe
  m_steel=v_steel*rho_steel
  ex_steel=m_steel*ex_steel_pipe/t_life_pipe # J/s for pipe material
  println("the pipe exergy requirement is $ex_steel J/s")
  dp_pipe=p[1]-p[end]

  # compression calculations
  # one compression stage
  p_0=1e5 # CO2 is delivered at this pressure
  T_comp=40.0+273.15 # K gas enters the compressor
  comp_ratio = 3.5 # compression ratio
  n_stage=Int(round(log(p[1]/p_0)/log(comp_ratio)))
  p_comp=zeros(n_stage+1)
  p_comp[1]=p_0
  p_comp[end]=p[1]
  for i in 2:n_stage
    p_comp[i]=p_comp[i-1]*comp_ratio
  end

  # Compression loop, with inter-cooling stages
  w_min_transport=0.0
  for i in 1:n_stage
    S_in=CoolProp.PropsSI("SMOLAR", "T", T_comp, "P", p_comp[i], gas_type) # molar entropy J/(mol.K)
    H_in=CoolProp.PropsSI("HMOLAR", "T", T_comp, "P", p_comp[i], gas_type) # molar enthalpy J/mol
    # find the output T for isentropic compressor
    #S_in-CoolProp.PropsSI("SMOLAR", "T", T_comp, "P", p_out, gas_type)
    T_out=CoolProp.PropsSI("T", "SMOLAR", S_in, "P", p_comp[i+1], gas_type)
    H_out=CoolProp.PropsSI("HMOLAR", "T", T_out, "P", p_comp[i+1], gas_type) # molar enthalpy J/mol
    w_min_transport+=H_out-H_in # J/mol
  end
  eta_comp=0.7
  eta_driver=0.9
  eta_pp=0.4
  w_real_transport=Q_g_molar*w_min_transport/(eta_comp*eta_driver*eta_pp)
  println("Compression exergy for gass transmission with a flow of $Q_g m3/s in a $d_inch inch pipe line is $w_real_transport W")
  # compression loop for CO2 injection
  p_wellhead=p[end]
  n_stage_well=Int(round(log(p_res/p_wellhead)/log(comp_ratio)))
  p_comp_well=zeros(n_stage_well+1)
  p_comp_well[1]=p_wellhead
  p_comp_well[end]=p_res
  for i in 2:n_stage_well
    p_comp_well[i]=p_comp_well[i-1]*comp_ratio
  end

  # Compression loop, with inter-cooling stages, for CO2 injection
  w_min_well=0.0
  for i in 1:n_stage_well
    S_in=CoolProp.PropsSI("SMOLAR", "T", T_comp, "P", p_comp_well[i], gas_type) # molar entropy J/(mol.K)
    H_in=CoolProp.PropsSI("HMOLAR", "T", T_comp, "P", p_comp_well[i], gas_type) # molar enthalpy J/mol
    # find the output T for isentropic compressor
    #S_in-CoolProp.PropsSI("SMOLAR", "T", T_comp, "P", p_out, gas_type)
    T_out=CoolProp.PropsSI("T", "SMOLAR", S_in, "P", p_comp_well[i+1], gas_type)
    H_out=CoolProp.PropsSI("HMOLAR", "T", T_out, "P", p_comp_well[i+1], gas_type) # molar enthalpy J/mol
    w_min_well+=H_out-H_in # J/mol
  end
  eta_comp=0.7
  eta_driver=0.9
  eta_pp=0.4
  w_real_well=Q_g_molar*w_min_well/(eta_comp*eta_driver*eta_pp)
  println("Compression exergy for gas injection with a flow of $Q_g m3/s and a reservoir pressure of $p_res Pa is $w_real_well W")
  push!(results_out["flow rate"], Q_g)
  push!(results_out["pipe diamter"], d_inch)
  push!(results_out["transport comp stages"], n_stage)
  push!(results_out["transport exergy"], w_real_transport)
  push!(results_out["steel exergy"], ex_steel)
  push!(results_out["injection exergy"], w_real_well)
  push!(results_out["preessure drop bar"], dp_pipe/1e5)
  push!(results_out["pipe length"], L_pipe)
  push!(results_out["pipe thickness"], t_pipe)
  push!(results_out["transport exergy th"], Q_g_molar*w_min_transport)
  push!(results_out["injection exergy th"], Q_g_molar*w_min_well)
end
d_results=DataFrame(results_out)
writetable("result_$gas_type.csv", d_results)

# writedlm("result_$gas_type.csv", results_out, ",")
# rho_gas=CoolProp.PropsSI("DMOLAR", "T", temperature_in_K,
#   "P", pressure_in_Pa, gas_type) # gas density mol/m3
# s_gas=CoolProp.PropsSI("SMOLAR", "T", temperature_in_K,
#   "P", pressure_in_Pa, gas_type) # molar density J/(mol.K)
# h_gas=CoolProp.PropsSI("HMOLAR", "T", temperature_in_K,
#   "P", pressure_in_Pa, gas_type) # gas enthaply J/mol
# mu_gas=CoolProp.PropsSI("VISCOSITY", "T", temperature_in_K,
#   "P", pressure_in_Pa, gas_type) # gas viscosity Pa.s
# Z_gas=CoolProp.PropsSI("Z", "T", temperature_in_K,
#   "P", pressure_in_Pa, gas_type) # gas compressibility factor
#
# p_sat_gas=CoolProp.PropsSI("P","T",T_in_K,"Q",0.0,gas_type)

import math


ETHANOL_FLOW_RATE = 55 # flow sensor Flow rate in LPM
GOX_FLOW_RATE = 250 # Gas regulator  flow in LPM
DENSITY_ETHANOL = 789  # kg/m^3
DIAMETER_PIPE = 20 / 1000  # 20 mm to m
PIPE_LENGTH_GOX =  2.046 # in meters
PIPE_LENGTH_ETHANOL = 2.030 #in meters
TEMPERATURE = 300  # K (Determined by Brian's Simulation)
FEED_PRESSURE = 600000 # Pa

# Fittings
K_ELBOWS = 0.9  # 0.3 - 0.9
K_T_SECTION = 0.6  # 0.2 - 0.6
ELBOWS_ETHANOL = 6
ELBOWS_GOX = 6
T_SECTIONS_ETHANOL = 1
T_SECTIONS_GOX = 1


# Valves
SOLENOID_VALVE_CV = 5.84
CHECK_VALVE_CV = 2.2
BALL_VALVE_CV = 1.5
SOLENOID_VALVES_ETHANOL = 1
SOLENOID_VALVES_GOX = 1
CHECK_VALVES_ETHANOL = 1
CHECK_VALVES_GOX = 1
BALL_VALVES_ETHANOL = 1
BALL_VALVES_GOX = 1

def convert_flows_to_m3_s(ethanol_lpm, gox_lpm):
    """
    Convert ethanol and GOX flow rates from LPM to m³/s.
    
    Parameters:
        ethanol_lpm (float): Ethanol flow in LPM
        gox_lpm (float): Gaseous oxygen flow in LPM

    Returns:
        tuple: (ethanol_flow_m3_s, gox_flow_m3_s)
    """
    def lpm_to_m3_per_s(lpm):
        return lpm / 1000 / 60

    ethanol_flow = lpm_to_m3_per_s(ethanol_lpm)
    gox_flow = lpm_to_m3_per_s(gox_lpm)

    return ethanol_flow, gox_flow

def pressure_loss_in_valves(Q, Cv, G):
    """
    Calculate the pressure loss in a valve.

    Parameters
    ----------
    Q : float
        Flow rate (m^3/s).
    Cv : float
        Valve flow coefficient.
    G : float
        Specific gravity of the fluid.

    Returns
    -------
    delta_P : float
        Pressure loss in the valve (Pa).

    """
    delta_P = (Q / Cv)**2 * G
    return delta_P


def calculate_gas_density(pressure: float = FEED_PRESSURE, temperature: float = TEMPERATURE) -> float:
    """
    Calculate the density of a gas.

    Parameters
    ----------
    temperature : float
        Temperature of the gas (K).
    pressure : float
        Pressure of the gas (Pa).

    Returns
    -------
    rho : float
        Density of the gas (kg/m^3).

    """
    R_o2 = 259.8  # J/kg K (Oxygen)
    rho = pressure / (R_o2 * temperature)
    return rho

def calculate_mass_flow_rate(volumetric_flow_m3_s, density_kg_m3):
    """
    Calculate mass flow rate (kg/s) from volumetric flow rate (m³/s) and density (kg/m³).

    Parameters:
        volumetric_flow_m3_s : float
            Flow rate in m³/s
        density_kg_m3 : float
            Density of the fluid in kg/m³

    Returns:
        float : Mass flow rate in kg/s
    """
    return volumetric_flow_m3_s * density_kg_m3

def calculate_specific_gravity(fluid_type: str, rho: float) -> float:
    """
    Calculate the specific gravity of the fluid.

    Returns
    -------
    G : float
        Specific gravity of the fluid.

    """
    RHO_AIR = 1.225
    RHO_WATER = 1000

    if fluid_type == 'liquid':
        G = rho / RHO_WATER
    elif fluid_type == 'gas':
        G = rho / RHO_AIR

    return G


def calculate_volumetric_flow_rate(m_dot: float, rho: float) -> float:
    """
    Calculate the volumetric flow rate.

    Parameters
    ----------
    m_dot : float
        Mass flow rate (kg/s).
    rho : float
        Density of the fluid (kg/m**3).

    Returns
    -------
    Q : float
        Volumetric flow rate (GPM).

    """
    Q = (m_dot / rho) * 15850

    return Q


def calculate_velocities(Q, D: float = DIAMETER_PIPE) -> float:
    """
    Calculate the velocity of the fluid in the pipe.

    Parameters
    ----------
    Q : float
        Flow rate (m^3/s).
    D : float
        Diameter of the pipe (m).

    Returns
    -------
    v : float
        Velocity of the fluid in the pipe (m/s).

    """
    A = (math.pi * D**2) / 4
    v = Q / A
    return v


def calculate_pipe_pressure_loss_darcy_weisbach(v, D, L, rho, f: float = 0.024) -> float:
    """
    Calculate the pressure loss in a pipe using the Darcy-Weisbach equation.

    Parameters
    ----------
    v : float
        Velocity of the fluid in the pipe (m/s).
    D : float
        Diameter of the pipe (m).
    L : float
        Length of the pipe (m).
    f : float
        Darcy-Weisbach friction factor.

    Returns
    -------
    delta_P : float
        Pressure loss in the pipe (Pa).

    """
    delta_P = f * rho * (L / D) * (v**2) / 2

    return delta_P


def calculate_fitting_pressure_loss(v, K, rho) -> float:
    """
    Calculate the pressure loss in a fitting.

    Parameters
    ----------
    v : float
        Velocity of the fluid in the pipe (m/s).
    K : float
        Loss coefficient of the fitting.
    rho : float
        Density of the fluid (kg/m**3).

    Returns
    -------
    delta_P : float
        Pressure loss in the fitting (Pa).

    """
    delta_P = K * rho * (v**2) / 2

    return delta_P


ethanol_flow_m3_s, gox_flow_m3_s = convert_flows_to_m3_s(ETHANOL_FLOW_RATE, GOX_FLOW_RATE)

density_gox = calculate_gas_density()

MASS_FLOW_RATE_ETHANOL = calculate_mass_flow_rate(ethanol_flow_m3_s, DENSITY_ETHANOL)
MASS_FLOW_RATE_GOX = calculate_mass_flow_rate(gox_flow_m3_s, density_gox)

volumetric_flow_rate_ethanol = calculate_volumetric_flow_rate(
    MASS_FLOW_RATE_ETHANOL, DENSITY_ETHANOL)
volumetric_flow_rate_gox = calculate_volumetric_flow_rate(
    MASS_FLOW_RATE_GOX, calculate_gas_density())

volumetric_flow_rate_ethanol_m3s = volumetric_flow_rate_ethanol * 0.00006309019640344
volumetric_flow_rate_gox_m3s = volumetric_flow_rate_gox * 0.00006309019640344


ethanol_velocity = calculate_velocities(volumetric_flow_rate_ethanol_m3s)
gox_velocity = calculate_velocities(volumetric_flow_rate_gox_m3s)

ethanol_line_pressure_loss = calculate_pipe_pressure_loss_darcy_weisbach(
    ethanol_velocity, DIAMETER_PIPE, PIPE_LENGTH_ETHANOL, DENSITY_ETHANOL)
gox_line_pressure_loss = calculate_pipe_pressure_loss_darcy_weisbach(
    gox_velocity, DIAMETER_PIPE, PIPE_LENGTH_GOX, calculate_gas_density())

ethanol_fitting_pressure_loss_elbows = calculate_fitting_pressure_loss(
    ethanol_velocity, K_ELBOWS, DENSITY_ETHANOL) * ELBOWS_ETHANOL
gox_fitting_pressure_loss_elbows = calculate_fitting_pressure_loss(
    gox_velocity, K_ELBOWS, calculate_gas_density()) * ELBOWS_GOX
ethanol_fitting_pressure_loss_t_sections = calculate_fitting_pressure_loss(
    ethanol_velocity, K_T_SECTION, DENSITY_ETHANOL) * T_SECTIONS_ETHANOL
gox_fitting_pressure_loss_t_sections = calculate_fitting_pressure_loss(
    gox_velocity, K_T_SECTION, calculate_gas_density()) * T_SECTIONS_GOX
ethanol_fitting_losses = ethanol_fitting_pressure_loss_elbows + \
    ethanol_fitting_pressure_loss_t_sections
gox_fitting_losses = gox_fitting_pressure_loss_elbows + \
    gox_fitting_pressure_loss_t_sections

specific_gravity_ethanol = calculate_specific_gravity(
    'liquid', DENSITY_ETHANOL)
specific_gravity_gox = calculate_specific_gravity(
    'gas', calculate_gas_density())

ethanol_valve_pressure_loss_solenoid = pressure_loss_in_valves(
    volumetric_flow_rate_ethanol_m3s, SOLENOID_VALVE_CV, specific_gravity_ethanol) * SOLENOID_VALVES_ETHANOL
gox_valve_pressure_loss_solenoid = pressure_loss_in_valves(
    volumetric_flow_rate_gox_m3s, SOLENOID_VALVE_CV, specific_gravity_gox) * SOLENOID_VALVES_GOX
ethanol_valve_pressure_loss_check = pressure_loss_in_valves(
    volumetric_flow_rate_ethanol_m3s, CHECK_VALVE_CV, specific_gravity_ethanol) * CHECK_VALVES_ETHANOL
gox_valve_pressure_loss_check = pressure_loss_in_valves(
    volumetric_flow_rate_gox_m3s, CHECK_VALVE_CV, specific_gravity_gox) * CHECK_VALVES_GOX
ethanol_valve_pressure_loss_ball = pressure_loss_in_valves(
    volumetric_flow_rate_ethanol_m3s, BALL_VALVE_CV, specific_gravity_ethanol) * BALL_VALVES_ETHANOL
gox_valve_pressure_loss_ball = pressure_loss_in_valves(
    volumetric_flow_rate_gox_m3s, BALL_VALVE_CV, specific_gravity_gox) * BALL_VALVES_GOX
ethanol_valve_losses = ethanol_valve_pressure_loss_solenoid + \
    ethanol_valve_pressure_loss_check + ethanol_valve_pressure_loss_ball
gox_valve_losses = gox_valve_pressure_loss_solenoid + \
    gox_valve_pressure_loss_check + gox_valve_pressure_loss_ball


total_ethanol_loss = ethanol_line_pressure_loss + \
    ethanol_fitting_losses + ethanol_valve_losses
total_gox_loss = gox_line_pressure_loss + gox_fitting_losses + gox_valve_losses

print("--- Pressure Losses ---\n")

print("INPUTS:")
print(f"Mass Flow Rate Ethanol: {ETHANOL_FLOW_RATE}  LPM")
print(f"Mass Flow Rate GOX: {GOX_FLOW_RATE}  LPM")
print(f"Pipe Diameter: {DIAMETER_PIPE * 100} mm -> {DIAMETER_PIPE} m")
print(f"Pipe Length GOX: {PIPE_LENGTH_GOX} m")
print(f"Pipe Length ETHANOL: {PIPE_LENGTH_ETHANOL} m")
print(f"Temperature: {TEMPERATURE} K")
print(f"Feed Pressure: {FEED_PRESSURE} Pa -> {FEED_PRESSURE / 100000} bar \n")

print("--- RESULTS ---\n")

print(
    f"Volumetric Flow Rate Ethanol: {volumetric_flow_rate_ethanol} GPM -> {volumetric_flow_rate_ethanol * 0.00006309019640344} m^3/s")
print(
    f"Volumetric Flow Rate GOX: {volumetric_flow_rate_gox} GPM -> {volumetric_flow_rate_gox * 0.00006309019640344} m^3/s\n")

print(f"Velocity Ethanol: {ethanol_velocity} m/s")
print(f"Velocity GOX: {gox_velocity} m/s\n")

print(
    f"Pipe Pressure Loss Ethanol: {ethanol_line_pressure_loss} Pa -> {ethanol_line_pressure_loss / 1000} kPa -> {ethanol_line_pressure_loss / 100000} bar")
print(
    f"Pipe Pressure Loss GOX: {gox_line_pressure_loss} Pa -> {gox_line_pressure_loss / 1000} kPa -> {gox_line_pressure_loss / 100000} bar\n")

print(
    f"Fitting Pressure Loss Ethanol: {ethanol_fitting_losses} Pa -> {ethanol_fitting_losses / 1000} kPa -> {ethanol_fitting_losses / 100000} bar")
print(
    f"Fitting Pressure Loss GOX: {gox_fitting_losses} Pa -> {gox_fitting_losses / 1000} kPa -> {gox_fitting_losses / 100000} bar\n")

print(
    f"Valve Pressure Loss Ethanol: {ethanol_valve_losses} Pa -> {ethanol_valve_losses / 1000} kPa -> {ethanol_valve_losses / 100000} bar")
print(
    f"Valve Pressure Loss GOX: {gox_valve_losses} Pa -> {gox_valve_losses / 1000} kPa -> {gox_valve_losses / 100000} bar\n")

print(
    f"Total Ethanol Loss: {total_ethanol_loss} Pa -> {total_ethanol_loss / 1000} kPa -> {total_ethanol_loss / 100000} bar")
print(
    f"Total GOX Loss: {total_gox_loss} Pa -> {total_gox_loss / 1000} kPa -> {total_gox_loss / 100000} bar\n")

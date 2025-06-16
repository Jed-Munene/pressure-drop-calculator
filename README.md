# **Pressure Drop Calculator for Ethanol and GOX Feed System**

This script is used to calculate the **pressure losses** in the **ethanol** and **gaseous oxygen (GOX)** feed lines of a propulsion system. It models pressure drops from pipe friction, fittings (elbows and T-sections), and valves (solenoid, check, and ball valves) using physical and empirical equations.

---

## What’s Included

This version includes significant enhancements and added functionality:

## New Features

- **Dual flow system modeling**: Separate calculations for ethanol (liquid) and GOX (gas)
- **Dynamic gas density**: GOX density calculated from temperature and pressure via ideal gas law
- **Flow conversions**: Converts LPM to m³/s, then to mass and volumetric flow
- **Valve losses**: Incorporates Cv-based losses for solenoid, check, and ball valves
- **Fitting losses**: Accounts for pressure losses from elbows and T-sections using loss coefficients (K-values)
- **Velocity calculations**: Computes flow velocity in the pipe for each fluid
- **Detailed breakdown**: Pipe, valve, and fitting losses printed separately in Pa, kPa, and bar


---

## How to Use

1. Clone the repository or copy the script.
2. Open the script file and review the input parameters at the top.
3. Run the script using Python:

```bash
python pressure_loss_calculator.py
```

---

## Input Parameters

These are defined at the top of the script:

| Variable                                               | Description                    |
| ------------------------------------------------------ | ------------------------------ |
| `ETHANOL_FLOW_RATE`                                    | Ethanol flow rate (LPM)        |
| `GOX_FLOW_RATE`                                        | Gaseous oxygen flow rate (LPM) |
| `DENSITY_ETHANOL`                                      | Ethanol density (kg/m³)        |
| `DIAMETER_PIPE`                                        | Pipe inner diameter (m)        |
| `PIPE_LENGTH_GOX`                                      | GOX line length (m)            |
| `PIPE_LENGTH_ETHANOL`                                  | Ethanol line length (m)        |
| `K_ELBOWS`, `K_T_SECTION`                              | Loss coefficients for fittings |
| `SOLENOID_VALVE_CV`, `CHECK_VALVE_CV`, `BALL_VALVE_CV` | Cv values for valve types      |
| `TEMPERATURE`                                          | Gas temperature (K)            |
| `FEED_PRESSURE`                                        | System feed pressure (Pa)      |
| `ELBOWS_*`, `T_SECTIONS_*`, `*_VALVES_*`               | Count of each fitting/valve    |

---

## Output

The script prints:

- Input summary (flow rates, geometry, system conditions)
- Volumetric flow rate in GPM and m³/s
- Velocity of ethanol and GOX in m/s
- Pipe pressure losses for each fluid
- Fitting pressure losses for each fluid
- Valve pressure losses for each fluid
- **Total pressure loss** per fluid line (in Pa, kPa, and bar)

---

## Notes

- The calculation assumes incompressible flow for ethanol and ideal gas behavior for GOX.
- The Darcy friction factor is assumed constant (0.024); refine as needed for Reynolds number-specific analysis.
- Unit conversions and coefficients are derived from standard engineering references.

---


import math
import LinearInterpolateFromFile as LI
import math as m

writeToTerminal = False


#########################################################################################

def NumericLinearInterpolate(x, x0, x1, y0, y1):
    y = (y0 * (x1 - x) + y1 * (x - x0)) / (x1 - x0)
    return y


def WTT(text):
    if writeToTerminal:
        print(text)


#########################################################################################


def HeatTransferCoefficientInTheAnnularSpaceBetweenPipes(t_fluid, t_pipeWall, D, d, L, w_fluid, flow, surface):
    """
    :param t_fluid: [°C]
    :param t_pipeWall: [°C]
    :param D: [m]
    :param d: [m]
    :param L: [m]
    :param w_fluid: [m/s]
    :param flow: laminar, turbulent, tryBoth
    :param surface: inside, outer, tryBoth
    """

    WTT("HeatTransferCoefficientInTheAnnularSpaceBetweenPipes")
    f, f1, f2, Nu = 0, 0, 0, 0
    alfa = [0, 0]
    X = D - d
    innerToOuterDiameterRatio = d / D
    t_ref = (t_fluid + t_pipeWall) / 2

    WaterPhyiscalProperties_file = r"C:\Users\juven\PycharmProjects\LinearInterpolator\DataSets\WaterPhysicalProperties_1bar.txt"
    WaterPhyiscalProperties = LI.AutomatedRequest("WaterPhysicalProperties_1bar.txt", 0, t_ref)
    # t(°C) | ρ(kg / m³) | cp(kJ / (kg·K)) | 10³ *β(1 / K) | 10³ *λ(W / (m·K)) | 10⁶ *μ(Pa·s) | 10⁶ *ν(m² / s) | Pr |
    # ID: 0 | ID: 1      | ID: 2           | ID: 3         | ID: 4             | ID: 5        | ID: 6          | ID: 7 |
    rho = WaterPhyiscalProperties[1]  # [kg/m^3]
    c_p = WaterPhyiscalProperties[2]  # [kJ/(kg*K)]
    Beta = WaterPhyiscalProperties[3] * 10e-3  # [1/K]
    lam = WaterPhyiscalProperties[4] * 10e-3  # [W/(m*K)]
    mu = WaterPhyiscalProperties[5] * 10e-6  # [Pa*s]
    nu = WaterPhyiscalProperties[6] * 10e-6  # [m^2/s]
    Pr = WaterPhyiscalProperties[7]  # Pr=mu*c_p/lam [-]
    Pr_w = LI.AutomatedRequest("WaterPhysicalProperties_1bar.txt", 0, t_pipeWall)[7]  # Pr=mu*c_p/lam [-]

    WTT("rho: " + str(rho))
    WTT("c_p: " + str(c_p))
    WTT("Beta: " + str(Beta))
    WTT("lam: " + str(lam))
    WTT("mu: " + str(mu))
    WTT("nu: " + str(nu))
    WTT("Pr: " + str(Pr))

    Phi_T_fluid = pow((Pr / Pr_w), 0.14)

    Re = w_fluid * X / nu

    WTT("Re: " + str(Re))

    if flow != "laminar" and flow != "turbulent" and flow != "tryBoth":
        WTT("Flow Definition Invalid")
        return
    if surface != "tryBoth" and surface != "inside" and surface != "outer":
        WTT("Surface Definition Invalid")
        return

    WTT("Testing for Flow Definition of: " + flow)
    WTT("Testing for Surface Definition of: " + surface)
    if flow == "laminar" or flow == "tryBoth":
        if Re > 2300: WTT("Flow Numerically Not Laminar!")
        if not (0 < innerToOuterDiameterRatio < 1000 and 0.1 <= Pr <= 1000):
            WTT("The Provided Combination of Values are Not Valid for Laminar Flow Conditions")
            return
        if surface == "tryBoth":
            f1 = (4 - 0.012 / (0.02 + innerToOuterDiameterRatio)) * pow(innerToOuterDiameterRatio, 0.04)
            f2 = 1 + 0.14 * pow(innerToOuterDiameterRatio, 0.1)
        elif surface == "outer":
            f1 = 1.2 * pow(innerToOuterDiameterRatio, 0.5)
            f2 = 1 + 0.14 / m.sqrt(innerToOuterDiameterRatio)
        elif surface == "inside":
            f1 = 1.2 * pow(innerToOuterDiameterRatio, -0.8)
            f2 = 1 + 0.14 / m.sqrt(innerToOuterDiameterRatio)
        NuLam = Phi_T_fluid * (3.66 + f1 +
                                  f2 * (0.19 * pow(Re * Pr * X / L, 0.8))
                                  / (1 + 0.117 * pow(Re * Pr * X / L, 0.467)))
        alfa[0] = NuLam * lam / X

    if flow == "turbulent" or flow == "tryBoth":
        if Re < 2300: WTT("Flow Numerically Not Turbulent!")
        if not (0.1 < X / L and 0.5 <= Pr <= 2000):
            WTT("The Provided Combination of Values are Not Valid for Turbulent Flow Conditions")
            return
        if surface == "tryBoth":
            f = 1 - 0.14 * pow(innerToOuterDiameterRatio, 0.6) + 0.86 * pow(innerToOuterDiameterRatio, 0.84) \
                / (1 + innerToOuterDiameterRatio)
        elif surface == "outer":
            f = 1 - 0.14 * pow(innerToOuterDiameterRatio, 0.6)
        elif surface == "inside":
            f = 0.86 * pow(innerToOuterDiameterRatio, -0.16)
        zeta = pow((1.82 * math.log10(Re) - 1.64), -2)
        NuTur = Phi_T_fluid * f * (1 + pow(X / L, 2 / 3)) * \
                (((zeta / 8) * (Re - 1000) * Pr) / (1 + 12.7 * math.sqrt(zeta / 8) * (pow(Pr, 2 / 3) - 1)))
        alfa[1] = NuTur * lam / X

    WTT("[laminaris, turbulent]")
    WTT("alfa: " + str(alfa))
    WTT("//////////////////////////////////////////")
    return alfa


def HydrogenRateOfHeatFlow(V_dot, t_initial, t_work):
    """
    :param V_dot: [Nl/min] #The volumetric flow of the hydrogen
    :param t_initial: [°C] #The initial temperature of the input hydrogen
    :param t_work: [°C] #The working temperature of the fuel cell
    """

    WTT("HydrogenRateOfHeatFlow")

    M_H2 = 2.016 * 10e-3  # [kg/mol]
    T_ref = (t_initial + t_work) / 2 + 273.15  # [K]
    HydrogenPhysicalProperties_file = r"C:\Users\juven\PycharmProjects\LinearInterpolator\DataSets\HydrogenPhysicalProperties_atm.txt"
    HydrogenPhysicalProperties_file1_5bar = r"C:\Users\juven\PycharmProjects\LinearInterpolator\DataSets\HydrogenPhysicalProperties_1_5bar.txt"
    c_p_mol = LI.AutomatedRequest("HydrogenPhysicalProperties_atm.txt", 0, T_ref)[1]  # [J/(mol*K)]
    c_p = c_p_mol / M_H2  # [J/(kg*K)]
    #Unfortunately the previous file (HydrogenPhysicalProperties_atm.txt) doesn't contain everything we need
    rho = LI.AutomatedRequest("HydrogenPhysicalProperties_1_5bar.txt", 0, T_ref)[1]  # [kg/m^3]

    m_dot = (V_dot * rho) * 10e-3 / 60  # [kg/s]

    Q_dot = m_dot * c_p * (t_work - t_initial)

    WTT("T_ref: " + str(T_ref))
    WTT("c_p: " + str(c_p))
    WTT("rho: " + str(rho))
    WTT("m_dot: " + str(m_dot))
    WTT("Q_dot: " + str(Q_dot))
    WTT("//////////////////////////////////////////")

    return Q_dot


def ConvectiveHeatResistanceOfCylindricalWall(d, delta, L, lam):
    """
    :param d: [m]
    :param delta: [m]
    :param L: [m]
    :param lam: [W/m*K]
    """
    WTT("ConvectiveHeatResistanceOfCylindricalWall")

    R_CylindricalWall = math.log((d + 2 * delta) / d) / (2 * math.pi * L * lam)

    WTT(R_CylindricalWall)
    WTT("//////////////////////////////////////////")
    return R_CylindricalWall


def AirFlowAroundIsolatedCylinder(t_fluid, t_cylindricalWall, L, d, delta, w_veg):
    """
    :param t_fluid: [°C]
    :param t_cylindricalWall: [°C]
    :param L: [m]
    :param d: [m]
    :param delta: [m]
    :param w_veg: [m/s]
    """
    WTT("AirFlowAroundIsolatedCylinder")

    T_cylindricalWall = t_cylindricalWall + 273.15
    T_fluid = t_fluid + 273.15
    t_ref = (t_cylindricalWall + t_fluid) / 2
    PhyiscalPropertiesOfAir_file = r"C:\Users\juven\PycharmProjects\LinearInterpolator\DataSets\AirPhysicalProperties_1bar.txt"
    PhyiscalPropertiesOfAir = LI.AutomatedRequest("AirPhysicalProperties_1bar.txt", 0, t_ref)
    # t(°C) | ρ(kg / m³) | cp(kJ / (kg·K)) | 10³ *β(1 / K) | 10³ *λ(W / (m·K)) | 10⁶ *μ(Pa·s) | 10⁶ *ν(m² / s) | Pr |
    # ID: 0 | ID: 1      | ID: 2           | ID: 3         | ID: 4             | ID: 5        | ID: 6          | ID: 7 |

    rho = PhyiscalPropertiesOfAir[1]  # [kg/m^3]
    c_p = PhyiscalPropertiesOfAir[2]  # [kJ/(kg*K)]
    Beta = PhyiscalPropertiesOfAir[3] * 10e-3  # [1/K]
    lam = PhyiscalPropertiesOfAir[4] * 10e-3  # [W/(m*K)]
    mu = PhyiscalPropertiesOfAir[5] * 10e-6  # [Pa*s]
    nu = PhyiscalPropertiesOfAir[6] * 10e-6  # [m^2/s]
    Pr = PhyiscalPropertiesOfAir[7]  # Pr=mu*c_p/lam [-]

    WTT("rho: " + str(rho))
    WTT("c_p: " + str(c_p))
    WTT("Beta: " + str(Beta))
    WTT("lam: " + str(lam))
    WTT("mu: " + str(mu))
    WTT("nu: " + str(nu))
    WTT("Pr: " + str(Pr))

    # Formalas from "hokozles feladatgyujetmeny" exercise 1.2.3 as provided by BME EGR 
    Phi_T_fluid = pow(T_fluid / T_cylindricalWall, 0.12)
    Phi_Psi = 1  # Assume flow is perpendicular to the cylinder wall
    Pr = nu / (lam / rho * c_p)
    Re = w_veg * (d + 2 * delta) / nu

    Nu = Phi_T_fluid * Phi_Psi * (0.3 + 0.62 * pow(Re, 0.5) * pow(Pr, 0.33) /
                                (1 + pow(0.4 / pow(Pr, 2 / 3), 0.25)) *
                                pow(1 + pow(Re / 2.82 * 10e5, 5 / 8), 0.8))
    alfa = Nu * lam / (d + 2 * delta)

    WTT("alfa: " + str(alfa))
    WTT("//////////////////////////////////////////")
    return alfa


#########################################################################################

def CylindricalBubbleHumidifierThermalContactDiagramm(t_heatSource, t_surrounding, v_insulation, lam_insulation, w_surrounding):
    """
    :param v_insulation: [m]
    :param lam_insulation: [[W/(m*k)]]
    :return:
    """
    t_fluid = 23
    D = 140 * 10e-3
    d = 56 * 10e-3
    v1 = 2 * 10e-3
    v2 = 4 * 10e-3
    L = 100 * 10e-3
    w_fluid = 0.04
    A1_b = pow(d, 2) / 4 * m.pi * L
    A1_k = pow(d + v1, 2) / 4 * m.pi * L
    A2_b = pow(D, 2) / 4 * m.pi * L
    A2_k = pow(D + v2, 2) / 4 * m.pi * L
    V_dot = 5  # [Nl/min]
    lam_steel = 46  # [W/(m*k)]
    t_insulation = (t_heatSource+t_surrounding)/2  # [°C] (Assumed for iter 1)
    d_insulation = D + v2
    t_surroundingPrevious = 0
    t_insulationPrevious = 0

    # We have two simple lines of heat transfer.
    # Both lines can be traced out from the heating fluid.
    # The first line supplies heat to the input hydrogen:
    # The necessary heat flow for heating up the room temperature input hydrogen
    Q_dotH2 = HydrogenRateOfHeatFlow(V_dot, t_fluid, t_heatSource)
    iter = 0

    while abs(t_surroundingPrevious - t_surrounding) > 0.5 and abs(t_insulationPrevious - t_insulation) > 0.5 or iter == 100:
        # The second line of heat transfer starts once again 
        # from the heating fluid and connects to the surrounding air 
        # through two convective and two conductive resistances

        #The first step is calulating the convective heat transfer coefficients of 
        # the heating fluid and the steel walls
        alfa_heatingFluid = HeatTransferCoefficientInTheAnnularSpaceBetweenPipes(t_heatSource, t_surrounding, D, d, L, w_fluid, "tryBoth", "tryBoth")
        alfa_heatingFluid = alfa_heatingFluid[1]
        R1_conv = 1 / (A1_k * alfa_heatingFluid)
        R2_conv = 1 / (A2_b * alfa_heatingFluid)

        # Two conductive heat resistances arise when the heatflow 
        # must go through the wall and insulation material
        R2_wall = ConvectiveHeatResistanceOfCylindricalWall(D, v2, L, lam_steel)
        R_insulation = ConvectiveHeatResistanceOfCylindricalWall(d_insulation, v_insulation, L, lam_insulation)

        # We can then determine the convective heat transfer coefficient of
        # the outer wall and the surrounding air
        alfa_air = AirFlowAroundIsolatedCylinder(t_fluid, t_insulation, L, d_insulation, v_insulation, w_surrounding)
        R_airConv = 1 / ((d_insulation + 2 * v_insulation) * math.pi * L * alfa_air)

        # We then sum the heat resistances of the entire system.
        R_szum = R2_conv + R2_wall + R_insulation + R_airConv
        # Determine the heatflow through the system
        Q_dotAir = (t_heatSource - t_fluid) / R_szum
        # Now we can calculate the next iteration for the insulation and air temperature
        t_surroundingPrevious = t_surrounding
        t_surrounding = t_heatSource - Q_dotAir * R2_conv
        t_insulationPrevious = t_insulation
        t_insulation = Q_dotAir * R_airConv + t_fluid
        iter += 1
        WTT("iter" + str(iter) + " - t_insulation: " + str(t_insulation) + ", t_surrounding: " + str(t_surrounding)+ ", Q_dotH2: " + str(Q_dotH2))
        return [t_surrounding, t_insulation]


#########################################################################################

#########################################################################################
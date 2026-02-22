import math

#This work is based on formulas provided by the book Fuel Cell Fundamentals by Ryan O'Hayre et al.
#The relevant chapter and page numbers are found next to each new introduced formula after the formula number.

"""
# Verified test values
t_FC = 80  # [°C]
T_FC = t_FC + 273.15  # [K]
p_FC_H2 = 1.5  # [bar]
p_FC_Air = 1.5  # [bar]
a_w_Air = 0.3  # [-] least relative humidity of Air at room temperature
a_w_H2 = (0.7, 1)
v_M=0.125 # [mm] membrane thickness (t_M in literature)
z=0.125
"""

# Water diffusivity plot samples (lambda:D_lambda)
WD_plot = {0:0, 2: 0.8, 3: 2.9, 4: 1.5}


def WaterVaporActivity(t, p, x_H2O):
    #(4.35) Chapter 4, pg.138
    p_SAT = pow(10, -2.1794 + 0.02953 * t
                - 9.1837 * pow(10, -5) * t * t
                + 1.4454 * pow(10, -7) * t * t * t)
    # (4.37) Chapter 4, pg.139
    x_H2O_SAT = p_SAT / p
    a_w = x_H2O / x_H2O_SAT
    return [p_SAT, x_H2O_SAT, a_w]


def LambdaBoundaryCondition(a_w):
    lam = 0
    # (4.45) Chapter 4, pg.142
    if a_w < 1:
        lam = 0.043 + 17.18 * a_w \
              - 39.85 * a_w * a_w \
              + 36 * a_w * a_w * a_w
    elif a_w < 3:
        lam = 14 + 4 * (a_w - 1)
    return lam


def Diffusivity(lam, T):
    # T in Kelvin
    # (4.51) Chapter 4, pg.143
    if lam > 4:
        D_lam = math.exp(2416 * ((1 / 303) - (1 / T))) \
                * (2.563 - 0.33 * lam
                   + 0.0264 * lam * lam
                   - 0.000671 * lam * lam * lam)
        return D_lam
    elif lam >= 0:
        x0=0
        x1=2
        if lam>2:
            x0 = 2
            x1 = 3
        if lam > 3:
            x0 = 3
            x1 = 4
        y0 = WD_plot.get(x0)
        y1 = WD_plot.get(x1)
        D_lam = LinearInterpolate(lam, x0, x1, y0, y1)
        #print("linear "+str(x0)+" "+str(x1)+", "+str(y0)+" "+str(y1)+": "+str(lam))
        return D_lam
    else:
        return None


def SolveWaterContentFunction(z, j, M_m, a_w_AN, a_w_CAT, T_AN, T_CAT):
    n_SAT_drag = 2.5  # [-]
    Rho_dry = 0.00197  # [kg/cm^3]
    F = 96.485  # [C/mol]

    lam_AN = LambdaBoundaryCondition(a_w_AN)
    lam_CAT = LambdaBoundaryCondition(a_w_CAT)
    lam_AVG=(lam_AN+lam_CAT)/2

    #D_lam_AN = float(Diffusivity(lam_AN, T_AN))
    D_lam_CAT = float(Diffusivity(lam_CAT, T_CAT))

    # (4.53) Chapter 4, pg.144
    exp_mult = j * M_m * n_SAT_drag / (22 * F * Rho_dry * D_lam_CAT) * 1000
    #print(exp_mult)
    try:
        C = (lam_CAT - lam_AN) / (math.exp(z * exp_mult) - 1)
        alfa = (lam_AN - C) / 4.4
    except: return 6
    thicknessList=[]
    lamList=[]
    for i in range(0,100):
        # (4.54) Chapter 4, pg.144
        lam = 4.4*alfa+C*math.exp(exp_mult*z/100*i)
        thicknessList.append(z/100*i)
        lamList.append(lam)

    """
    # Debugging variables
    print("exp_mlt: " + str(exp_mult))
    print("alfa: " + str(alfa))
    print("C: " + str(C))
    print("lam_AN: " + str(lam_AN))
    #print("lam_CAT: " + str(lam_CAT))
    #print("D_lam: " + str(D_lam_CAT))
    """

    return alfa


def LinearInterpolate(x, x0, x1, y0, y1):
    y = (y0 * (x1 - x) + y1 * (x - x0)) / (x1 - x0)
    return y

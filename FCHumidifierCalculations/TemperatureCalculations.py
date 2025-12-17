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


def HoatadasiTenyezoCsovekKozottiGyurusTerben(t_kozeg, t_csofal, D, d, L, w_kozeg, aramlas, felulet):
    """
    :param t_kozeg: [°C]
    :param t_csofal: [°C]
    :param D: [m]
    :param d: [m]
    :param L: [m]
    :param w_kozeg: [m/s]
    :param aramlas: lamináris, turbulens, mindkettő
    :param felulet: belső, külső, mindkettő
    """

    WTT("HoatadasiTenyezoCsovekKozottiGyurusTerben")
    f, f1, f2, Nu = 0, 0, 0, 0
    alfa = [0, 0]
    X = D - d
    atmeroArany = d / D
    t_ref = (t_kozeg + t_csofal) / 2

    VizTablazat = r"C:\Users\juven\PycharmProjects\LinearInterpolator\DataSets\WaterPhysicalProperties_1bar.txt"
    vizFizikaiTulajdonsagai = LI.AutomatedRequest(VizTablazat, 0, t_ref)
    # t(°C) | ρ(kg / m³) | cp(kJ / (kg·K)) | 10³ *β(1 / K) | 10³ *λ(W / (m·K)) | 10⁶ *μ(Pa·s) | 10⁶ *ν(m² / s) | Pr |
    # ID: 0 | ID: 1      | ID: 2           | ID: 3         | ID: 4             | ID: 5        | ID: 6          | ID: 7 |
    rho = vizFizikaiTulajdonsagai[1]  # [kg/m^3]
    c_p = vizFizikaiTulajdonsagai[2]  # [kJ/(kg*K)]
    Beta = vizFizikaiTulajdonsagai[3] * 10e-3  # [1/K]
    lam = vizFizikaiTulajdonsagai[4] * 10e-3  # [W/(m*K)]
    mu = vizFizikaiTulajdonsagai[5] * 10e-6  # [Pa*s]
    nu = vizFizikaiTulajdonsagai[6] * 10e-6  # [m^2/s]
    Pr = vizFizikaiTulajdonsagai[7]  # Pr=mu*c_p/lam [-]
    Pr_w = LI.AutomatedRequest(VizTablazat, 0, t_csofal)[7]  # Pr=mu*c_p/lam [-]

    WTT("rho: " + str(rho))
    WTT("c_p: " + str(c_p))
    WTT("Beta: " + str(Beta))
    WTT("lam: " + str(lam))
    WTT("mu: " + str(mu))
    WTT("nu: " + str(nu))
    WTT("Pr: " + str(Pr))

    Phi_T_folyadek = pow((Pr / Pr_w), 0.14)

    Re = w_kozeg * X / nu

    WTT("Re: " + str(Re))

    if aramlas != "lamináris" and aramlas != "turbulens" and aramlas != "mindkettő":
        WTT("Helytelen aramlás")
        return
    if felulet != "mindkettő" and felulet != "belső" and felulet != "külső":
        WTT("Helytelen felület")
        return

    WTT("Vizsgált áramlási forma: " + aramlas)
    WTT("Vizsgált hőatadó felület: " + felulet)
    if aramlas == "lamináris" or aramlas == "mindkettő":
        if Re > 2300: WTT("Numerikusan nem lamináris!!!")
        if not (0 < atmeroArany < 1000 and 0.1 <= Pr <= 1000):
            WTT("Az adatok nem érvényesek lamináris áramláshoz")
            return
        if felulet == "mindkettő":
            f1 = (4 - 0.012 / (0.02 + atmeroArany)) * pow(atmeroArany, 0.04)
            f2 = 1 + 0.14 * pow(atmeroArany, 0.1)
        elif felulet == "külső":
            f1 = 1.2 * pow(atmeroArany, 0.5)
            f2 = 1 + 0.14 / m.sqrt(atmeroArany)
        elif felulet == "belső":
            f1 = 1.2 * pow(atmeroArany, -0.8)
            f2 = 1 + 0.14 / m.sqrt(atmeroArany)
        NuLam = Phi_T_folyadek * (3.66 + f1 +
                                  f2 * (0.19 * pow(Re * Pr * X / L, 0.8))
                                  / (1 + 0.117 * pow(Re * Pr * X / L, 0.467)))
        alfa[0] = NuLam * lam / X

    if aramlas == "turbulens" or aramlas == "mindkettő":
        if Re < 2300: WTT("Numerikusan nem turbulens!!!")
        if not (0.1 < X / L and 0.5 <= Pr <= 2000):
            WTT("Az adatok nem érvényesek turbulens áramláshoz")
            return
        if felulet == "mindkettő":
            f = 1 - 0.14 * pow(atmeroArany, 0.6) + 0.86 * pow(atmeroArany, 0.84) \
                / (1 + atmeroArany)
        elif felulet == "külső":
            f = 1 - 0.14 * pow(atmeroArany, 0.6)
        elif felulet == "belső":
            f = 0.86 * pow(atmeroArany, -0.16)
        zeta = pow((1.82 * math.log10(Re) - 1.64), -2)
        NuTur = Phi_T_folyadek * f * (1 + pow(X / L, 2 / 3)) * \
                (((zeta / 8) * (Re - 1000) * Pr) / (1 + 12.7 * math.sqrt(zeta / 8) * (pow(Pr, 2 / 3) - 1)))
        alfa[1] = NuTur * lam / X

    WTT("[laminaris, turbulens]")
    WTT("alfa: " + str(alfa))
    WTT("//////////////////////////////////////////")
    return alfa


def HidrogenHoaramIgeny(V_pont, t_kezdo, t_elvart):
    """
    :param V_pont: [Nl/min]
    :param t_kezdo: [°C]
    :param t_elvart: [°C]
    """

    WTT("HidrogenHoaramIgeny")

    M_H2 = 2.016 * 10e-3  # [kg/mol]
    T_ref = (t_kezdo + t_elvart) / 2 + 273.15  # [K]
    hidrogenTablazat = r"C:\Users\juven\PycharmProjects\LinearInterpolator\DataSets\HydrogenPhysicalProperties_atm.txt"
    hidrogenTablazat1_5bar = r"C:\Users\juven\PycharmProjects\LinearInterpolator\DataSets\HydrogenPhysicalProperties_1_5bar.txt"
    c_p_mol = LI.AutomatedRequest(hidrogenTablazat, 0, T_ref)[1]  # [J/(mol*K)]
    c_p = c_p_mol / M_H2  # [J/(kg*K)]
    rho = LI.AutomatedRequest(hidrogenTablazat1_5bar, 0, T_ref)[1]  # [kg/m^3]

    m_pont = (V_pont * rho) * 10e-3 / 60  # [kg/s]

    Q_pont = m_pont * c_p * (t_elvart - t_kezdo)

    WTT("T_ref: " + str(T_ref))
    WTT("c_p: " + str(c_p))
    WTT("rho: " + str(rho))
    WTT("m_pont: " + str(m_pont))
    WTT("Q_pont: " + str(Q_pont))
    WTT("//////////////////////////////////////////")

    return Q_pont


def HengeresFalKonvektivHoellenallasa(d, delta, L, lam):
    """
    :param d: [m]
    :param delta: [m]
    :param L: [m]
    :param lam: [W/m*K]
    """
    WTT("HengeresFalKonvektivHoellenallasa")

    R_heng = math.log((d + 2 * delta) / d) / (2 * math.pi * L * lam)

    WTT(R_heng)
    WTT("//////////////////////////////////////////")
    return R_heng


def EgyedulalloHengerKorulAramloLevego(t_veg, t_fal, L, d, delta, w_veg):
    """
    :param t_veg: [°C]
    :param t_fal: [°C]
    :param L: [m]
    :param d: [m]
    :param delta: [m]
    :param w_veg: [m/s]
    """
    WTT("EgyedulalloHengerKoruliAramlasLevego")

    T_fal = t_fal + 273.15
    T_veg = t_veg + 273.15
    t_ref = (t_fal + t_veg) / 2
    levegoTablazat = r"C:\Users\juven\PycharmProjects\LinearInterpolator\DataSets\AirPhysicalProperties_1bar.txt"
    levegoFizikaiTulajdonsagai = LI.AutomatedRequest(levegoTablazat, 0, t_ref)
    # t(°C) | ρ(kg / m³) | cp(kJ / (kg·K)) | 10³ *β(1 / K) | 10³ *λ(W / (m·K)) | 10⁶ *μ(Pa·s) | 10⁶ *ν(m² / s) | Pr |
    # ID: 0 | ID: 1      | ID: 2           | ID: 3         | ID: 4             | ID: 5        | ID: 6          | ID: 7 |

    rho = levegoFizikaiTulajdonsagai[1]  # [kg/m^3]
    c_p = levegoFizikaiTulajdonsagai[2]  # [kJ/(kg*K)]
    Beta = levegoFizikaiTulajdonsagai[3] * 10e-3  # [1/K]
    lam = levegoFizikaiTulajdonsagai[4] * 10e-3  # [W/(m*K)]
    mu = levegoFizikaiTulajdonsagai[5] * 10e-6  # [Pa*s]
    nu = levegoFizikaiTulajdonsagai[6] * 10e-6  # [m^2/s]
    Pr = levegoFizikaiTulajdonsagai[7]  # Pr=mu*c_p/lam [-]

    WTT("rho: " + str(rho))
    WTT("c_p: " + str(c_p))
    WTT("Beta: " + str(Beta))
    WTT("lam: " + str(lam))
    WTT("mu: " + str(mu))
    WTT("nu: " + str(nu))
    WTT("Pr: " + str(Pr))

    # Képletek a "hokozles_fealdatgyujetmeny.pdf" 1.2.3 feladata alapján
    Phi_T_gaz = pow(T_veg / T_fal, 0.12)
    Phi_Psi = 1  # Az áramlás merőleges a hengerre
    Pr = nu / (lam / rho * c_p)
    Re = w_veg * (d + 2 * delta) / nu

    Nu = Phi_T_gaz * Phi_Psi * (0.3 + 0.62 * pow(Re, 0.5) * pow(Pr, 0.33) /
                                (1 + pow(0.4 / pow(Pr, 2 / 3), 0.25)) *
                                pow(1 + pow(Re / 2.82 * 10e5, 5 / 8), 0.8))
    alfa = Nu * lam / (d + 2 * delta)

    WTT("alfa: " + str(alfa))
    WTT("//////////////////////////////////////////")
    return alfa


#########################################################################################

def ParasitoHokapcsolas(v_szig, lam_szig):
    """
    :param v_szig: [m]
    :param lam_szig: [[W/(m*k)]]
    :return:
    """
    t_futo = 80
    t_k = 75  # [°C] (becsült)
    t_veg = 23
    T_futo = 273.15 + t_futo
    T_csofal = 273.15 + t_veg
    D = 140 * 10e-3
    d = 56 * 10e-3
    v1 = 2 * 10e-3
    v2 = 4 * 10e-3
    L = 100 * 10e-3
    w_kozeg = 0.04
    A1_b = pow(d, 2) / 4 * m.pi * L
    A1_k = pow(d + v1, 2) / 4 * m.pi * L
    A2_b = pow(D, 2) / 4 * m.pi * L
    A2_k = pow(D + v2, 2) / 4 * m.pi * L
    V_pont = 5  # [Nl/min]
    lam_acel = 46  # [W/(m*k)]

    d_szig = D + v2
    # v_szig = 40 * 10e-3 #[W/(m*k)]
    # lam_szig = 0.35  # [W/(m*k)] Egyenlőre ez a polietilén
    w_lev = 0.004  # [m/s]
    t_szig = 40  # [°C] (becsült)

    t_k_elozo = 0
    t_szig_elozo = 0

    # szoba hőmérsékletű hidrogén felveszi a víz hőmérsékletét
    Q_pontH2 = HidrogenHoaramIgeny(V_pont, t_veg, t_futo)
    iter = 0

    while abs(t_k_elozo - t_k) > 0.5 and abs(t_szig_elozo - t_szig) > 0.5 or iter == 100:
        # Az áramló közeg átadja a hőmérsékletét a két falnak konvektív hőátadással
        alfa_fut = HoatadasiTenyezoCsovekKozottiGyurusTerben(t_futo, t_k, D, d, L, w_kozeg, "mindkettő", "mindkettő")
        alfa_fut = alfa_fut[1]
        R1_konv = 1 / (A1_k * alfa_fut)
        R2_konv = 1 / (A2_b * alfa_fut)

        R2_fal = HengeresFalKonvektivHoellenallasa(D, v2, L, lam_acel)
        R_szig = HengeresFalKonvektivHoellenallasa(d_szig, v_szig, L, lam_szig)

        # A készülék külső falán a szigetelés és a levegő között is konvektív hőátadás történik
        alfa_lev = EgyedulalloHengerKorulAramloLevego(t_veg, t_szig, L, d_szig, v_szig, w_lev)
        Rlev_konv = 1 / ((d_szig + 2 * v_szig) * math.pi * L * alfa_lev)

        # A teljes hőkapcsolás összeáll
        R_szum = R2_konv + R2_fal + R_szig + Rlev_konv
        # Kiértékeljük a közeg és a fűtőtér közötti össz hőáramot
        Q_pont_lev = (t_futo - t_veg) / R_szum
        # Új becsléseket veszünk fel a konvektív hőátadás falaira
        t_k_elozo = t_k
        t_k = t_futo - Q_pont_lev * R2_konv
        t_szig_elozo = t_szig
        t_szig = Q_pont_lev * Rlev_konv + t_veg
        iter += 1
        WTT("iter" + str(iter) + " - t_szig: " + str(t_szig) + ", t_k: " + str(t_k))
    return [t_k, t_szig]


#########################################################################################

ParasitoHokapcsolas(4 * 10e-3, 35)

#########################################################################################
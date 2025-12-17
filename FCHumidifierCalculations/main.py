import OptimalWaterVolumeCalculations as OC
import TemperatureCalculations as TC
import DisplayPlots as DP

T = 353  # [K]
z = 0.0125  # [cm]
j = 0.7  # [A/cm^2]
M_m = 1  # [kg/mol]


def IterateWaterActivityForAlfa():
    aw_H2List = list(range(0, 101))
    aw_airList = list(range(0, 101))
    alfaList = []
    for i in range(0,101):alfaList.append([0]*101)
    for aw_H2 in aw_H2List:
        for aw_air in aw_airList:
            alfa = OC.SolveWaterContentFunction(z, j, M_m, aw_H2 / 100, aw_air / 100, T, T)
            alfaList[aw_H2][aw_air]=alfa
            if aw_H2==30 and aw_air==30: print("aw_H2=30% esetén alfa="+str(alfa))
            if aw_H2 == 65 and aw_air == 30: print("aw_H2=65% esetén alfa="+str(alfa))
    #for i in range(0, 101): print(alfaList[i])
    #print(OC.SolveWaterContentFunction(z, j, M_m, 0.8, 1, T, T))
    DP.Display3DPlot(aw_airList, "aw_air [%]", aw_H2List, "aw_H2 [%]", alfaList, "alfa [1]")


def CriticalThicknessOfInsulation():
    v_szig_list = []
    t_szig_list = []
    lam_szig = 0.35 * 10e-3  # [m] üveggyapot
    critFound = False
    for v_szig in range(1, 1000):  # 0.01-10 [mm]
        v_szig_list.append(v_szig / 100)  # int->[mm]
        v_szig = v_szig / 100 / 1000  # int->[mm]->[m]
        wallTemps = TC.ParasitoHokapcsolas(v_szig, lam_szig)
        t_szig = wallTemps[1]
        t_szig_list.append(t_szig)
        # print(t_szig)
        if t_szig < 24 and critFound == False:
            print(v_szig * 1000)    # [m] -> [mm]
            critFound = True

    DP.Display2DPlot(v_szig_list, "v_szig [mm]", t_szig_list, "t_szig [°C]")


def main():
    IterateWaterActivityForAlfa()
    CriticalThicknessOfInsulation()


if __name__ == '__main__':
    main()

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
    print("Validation samples: ")
    for i in range(0, 101): alfaList.append([0] * 101)
    for aw_H2 in aw_H2List:
        for aw_air in aw_airList:
            alfa = OC.SolveWaterContentFunction(z, j, M_m, aw_H2 / 100, aw_air / 100, T, T)
            alfaList[aw_H2][aw_air] = alfa
            if aw_H2 == 30 and aw_air == 30: print("at aw_H2=30%: alfa=" + str(alfa))
            if aw_H2 == 65 and aw_air == 30: print("at aw_H2=65%: alfa=" + str(alfa))
    # for i in range(0, 101): print(alfaList[i])
    # print(OC.SolveWaterContentFunction(z, j, M_m, 0.8, 1, T, T))
    DP.Display3DPlot(aw_airList, "aw_air [%]", aw_H2List, "aw_H2 [%]", alfaList, "alfa [1]")


def CriticalThicknessOfInsulation():
    v_insulation_list = []
    t_insulation_list = []
    lam_insulation = 0.35 * 10e-3  # [W/(m*k)] fiberglass
    t_heatSource = 80  # [°C]
    #t_surrounding = 45  # [°C] Assume avg. summer weather and temperature inside of black chassis
    t_surrounding = 23  # [°C] Assume room temperature
    #t_surrounding = 12  # [°C] Assume worst case for summer (Very high temp difference)
    t_bleed = 1 # [°C] Permissible temperature difference between the surrounding air and the insulation surface
    w_surrounding = 0.004  # [m/s] Assume inside chassis
    critFound = False
    for v_insulation in range(1, 1000):  # 0.01-10 [mm]
        v_insulation_list.append(v_insulation / 100)  # int->[mm]
        v_insulation = v_insulation / 100 / 1000  # int->[mm]->[m]
        wallTemps = TC.CylindricalBubbleHumidifierThermalContactDiagramm(t_heatSource, t_surrounding, v_insulation,
                                                                         lam_insulation, w_surrounding)
        t_insulation = wallTemps[1]
        t_insulation_list.append(t_insulation)
        # print(t_insulation)
        if t_insulation < t_surrounding + t_bleed and critFound == False:
            print("optimal insulation thickness: " + str(v_insulation * 1000) + " [mm]")  # [m] -> [mm]
            critFound = True

    DP.Display2DPlot(v_insulation_list, "v_insulation [mm]", t_insulation_list, "t_insulation [°C]")


def main():
    IterateWaterActivityForAlfa()
    CriticalThicknessOfInsulation()


if __name__ == '__main__':
    main()

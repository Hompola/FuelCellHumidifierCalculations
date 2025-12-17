import os
import sys
import displayData as d
from displayData import DisplayData
from displayData import WTT

g = d.g
w = d.w


def ListFiles():
    global parentDirectory
    global dataDirectory
    global dataSets
    # First lets find the folder that main.py is in:
    parentDirectory = os.getcwd()
    print(parentDirectory)
    # Now lets check for a DataSets folder, and create one if there is none:
    dataDirectory = os.path.join(parentDirectory + "\DataSets")
    try:
        os.mkdir(dataDirectory)
        WTT("A" + g + "'DataSets'" + w + " directory has been created in the folder of this programme!")
    except FileExistsError:
        pass

    # Then lets see if there are any data sets in the folder. if not, then let's set up some help
    dataSets = os.listdir(dataDirectory)
    if len(dataSets) == 0:
        SetUp()
        ListFiles()
    else:
        # List the data sets in the folder with and ID
        for i in range(0, len(dataSets)):
            WTT("ID:" + str(i) + " | " + dataSets[i] + "\n")


def ProcessFile(filePath):
    global dataArray
    global ColumnWidth
    manual = False
    if filePath == "choose":
        manual = True
        file = ProcessFileInputHandler()
    else:
        try:
            file = open(filePath, encoding="utf-8", errors="ignore")
        except:
            WTT("not a valid path")
            return

    data = file.read().strip().split("\n")
    file.close()
    # Remove all empty strings
    data = [x for x in data if x!='']

    dataArray = []
    # Fill up dataArray with each row
    for line in data:
        dataArray.append(line.split("\t"))

    # Since there are multiple options for a table to be made up, the code should handle a couple of structures.
    # One of these options is that the table contains headers at the top row and multiple columns
    # Data could also be ordered in a single column below each other, in a list format
    # Both will be converted into a matrix format, so only the list structure has to be altered

    # For nicely displaying the data with the monospace font of the terminal means we need the widest column
    # The starting value is set to 6, as a column ID will also be displayed following the "ID: ##" structure
    columnWidth = 6
    # We also determine  the length of the columns
    ColumnHeigth = 0
    if len(dataArray[0]) == 1:
        for i in range(1, len(data)):
            if len(data[i]) > columnWidth:
                columnWidth = len(data[i])
            try:
                dataArray[i] = float(data[i])
            except:
                ColumnHeigth = len(data) - i
                pass

        # Create an array of empty lists instead of the previous dataArray
        dataArray = [[] for i in range(ColumnHeigth)]
        # Fill the array with each row
        for i in range(0, len(data)):
            dataArray[i % ColumnHeigth].append(data[i])
    elif manual:
        DisplayData(dataArray)
    return dataArray


def ProcessFileInputHandler():
    while True:
        WTT("\nInput the ID of the " + g + "data set" + w + " you'd like to work in: ")
        fileID = input()
        try:
            fileID = int(fileID)
            file = open(os.path.join(dataDirectory, dataSets[fileID]), encoding="utf-8", errors="ignore")
        except ValueError:
            WTT("Looks like the ID isn't an integer! Try again.")
            continue  # Restart the loop
        except IndexError:
            WTT("The specified ID is out of range! Try again.")
            continue  # Restart the loop
        except Exception as e:
            WTT(f"An unexpected error occurred: {e}")
            WTT("Would you like to retry? [y/n]: ")
            if input().lower() != "y":
                quit()
            continue
        # Runs only if ID is valid
        return file


def Interpolate(array, columnID=None, rowValue=None):
    manual = False
    if columnID == "choose" or rowValue == "choose":
        manual = True
        arrayInfo = InterpolateInputHandler()
        ID = arrayInfo[0]
        x = arrayInfo[1]
    else:
        ID = columnID
        x = rowValue

    dataArray = array
    x0 = 0
    x0index = 0
    y0 = 0
    x1 = 0
    x1index = 0
    y1 = 0
    result = []

    # We'll convert all the values into floating point numbers to work with from now on
    values = [list(map(float, i)) for i in dataArray[1:len(dataArray)]]  # I'm quite proud of this one :3
    for i in range(len(values)):
        if x < values[i][ID]:
            x0index = i - 1
            x1index = i
            break
    x0 = values[x0index][ID]
    x1 = values[x1index][ID]
    for i in range(len(dataArray[0])):
        y0 = values[x0index][i]
        y1 = values[x1index][i]
        if x1 - x0 == 0:
            y = y0
        else:
            y = round((y0 * (x1 - x) + y1 * (x - x0)) / (x1 - x0), 5)
        result.append(y)
    if manual:
        DisplayData([dataArray[0], result])
    return result


def InterpolateInputHandler():
    while True:
        WTT("\nSpecify the ID of the " + g + "column" + w + " your search will be in: ")
        ID = input()
        try:
            ID = int(ID)
            WTT("Chosen column: " + dataArray[0][ID])
        except ValueError:
            WTT("Looks like the ID isn't an integer! Try again.")
            continue  # Restart the loop
        except IndexError:
            WTT("The specified ID is out of range! Try again.")
            continue  # Restart the loop
        except Exception as e:
            WTT(f"An unexpected error occurred: {e}")
            WTT("Would you like to retry? [y/n]: ")
            if input().lower() != "y":
                quit()
            continue
        # Runs only if ID is valid
        break  # Exit the loop after successful input

    while True:
        WTT("\nSpecify the " + g + "value" + w + ": ")
        value = input()
        try:
            value = float(value)
        except ValueError:
            WTT("Looks like the value isn't a number! Try again.")
            continue  # Restart the loop
        except Exception as e:
            WTT(f"An unexpected error occurred: {e}")
            WTT("Would you like to retry? [y/n]: ")
            if input().lower() != "y":
                quit()
            continue
        # Runs only if the value is valid
        return ID, value  # Exit the loop after successful input


def ManualRequest():
    ListFiles()
    dataArray = ProcessFile("choose")
    Interpolate(dataArray, "choose", "choose")


def AutomatedRequest(filePath, columnID, rowValue):
    dataArray = ProcessFile(filePath)
    result = Interpolate(dataArray, columnID, rowValue)
    return result


def main():
    ManualRequest()
    # filePath = r"C:\Users\juven\PycharmProjects\LinearInterpolator\DataSets\HydrogenPhysicalProperties1_5bar.txt"
    # print(AutomatedRequest(filePath, 0,100))
    return


def SetUp():
    WTT("\nLooks like there aren't any data sets to read in yet!")
    WTT("\nGo to this folder in order to start expanding it: " + dataDirectory)
    WTT("\nThere should be some example files there now, called DataSet1 and Dataset2!"
        "\nIf you dont find one, this programme didn't have the permission to create a folder! "
        "\nFeel free to do so yourself, make sure its called DataSets!"
        "\nYou will also find an INFO file in the programme folder\n\n\n")
    infofile = open(os.path.join(parentDirectory, "INFO.txt"), "w", encoding="utf-8", errors="ignore")
    infofile.write("Try opening the DataSets directory and filling it with .txt files that have tables of values "
                   "to interpolate between! \nUse text headers for each column! \nYou may also use a list format, "
                   "with the headers breaking up each column.\n\nYou have a couple of examples in your DataSets "
                   "folder already, if the programme made your directory. "
                   "\nOtherwise, with an empty DataSets folder, run the programme again to get some examples!"
                   "\nBoth DataSet1 and DataSet2 contain the Physical Properties of Dry Air at 1 Bar Pressure")
    infofile.close()

    DataSet1 = open(os.path.join(dataDirectory, "DataSet1.txt"), "w", encoding="utf-8", errors="ignore")
    DataSet1.write("t (°C)\tρ (kg/m³)\tcp (kJ/(kg·K))\t10³ * β (1/K)\t10³ * λ (W/(m·K))\t10⁶ * μ (Pa·s)\t10⁶ * ν ("
                   "m²/s)\tPr\n" +
                   "-180\t3.8515\t1.071\t11.071\t9.00\t6.44\t1.67\t0.77\n" +
                   "-160\t3.1258\t1.036\t9.320\t10.70\t7.58\t2.51\t0.75\n" +
                   "-140\t2.6391\t1.021\t7.758\t12.90\t9.20\t3.48\t0.74\n" +
                   "-120\t2.2867\t1.014\t6.659\t14.60\t10.49\t4.587\t0.73\n" +
                   "-100\t2.0186\t1.011\t5.846\t16.40\t11.76\t5.806\t0.72\n" +
                   "-80\t1.8073\t1.009\t5.219\t18.16\t12.89\t7.132\t0.72\n" +
                   "-60\t1.6364\t1.007\t4.719\t19.83\t14.02\t8.567\t0.71\n" +
                   "-40\t1.4952\t1.006\t4.304\t21.45\t15.09\t10.09\t0.71\n" +
                   "-20\t1.3765\t1.006\t3.962\t23.01\t16.15\t11.73\t0.71\n" +
                   "0\t1.2754\t1.006\t3.671\t24.54\t17.10\t13.41\t0.70\n" +
                   "20\t1.1881\t1.007\t3.419\t26.03\t17.98\t15.13\t0.70\n" +
                   "40\t1.1120\t1.008\t3.200\t27.49\t18.81\t16.92\t0.69\n" +
                   "60\t1.0452\t1.009\t3.007\t28.94\t19.73\t18.88\t0.69\n" +
                   "80\t0.9859\t1.010\t2.836\t30.38\t20.73\t21.30\t0.69\n" +
                   "100\t0.9329\t1.012\t2.684\t31.81\t21.60\t23.15\t0.69\n" +
                   "120\t0.8854\t1.014\t2.547\t33.22\t22.43\t25.33\t0.68\n" +
                   "140\t0.8425\t1.017\t2.423\t34.66\t23.19\t27.53\t0.68\n" +
                   "160\t0.8036\t1.020\t2.311\t36.07\t24.01\t29.88\t0.68\n" +
                   "180\t0.7681\t1.023\t2.209\t37.49\t24.91\t32.43\t0.68\n" +
                   "200\t0.7356\t1.026\t2.115\t38.91\t25.70\t34.94\t0.68\n" +
                   "250\t0.6653\t1.035\t1.912\t42.43\t27.40\t41.18\t0.67\n" +
                   "300\t0.6072\t1.046\t1.745\t45.91\t29.20\t48.09\t0.67\n" +
                   "350\t0.5585\t1.057\t1.605\t49.31\t30.90\t55.33\t0.66\n" +
                   "400\t0.5170\t1.069\t1.485\t52.57\t32.55\t62.95\t0.66\n" +
                   "450\t0.4813\t1.081\t1.383\t55.64\t34.00\t70.64\t0.66\n" +
                   "500\t0.4502\t1.093\t1.293\t58.48\t35.50\t78.86\t0.66\n" +
                   "600\t0.3968\t1.116\t1.145\t63.50\t38.30\t96.08\t0.67\n" +
                   "700\t0.3577\t1.137\t1.027\t67.80\t40.87\t114.3\t0.69\n" +
                   "800\t0.3243\t1.155\t0.932\t71.30\t43.32\t133.6\t0.70\n" +
                   "900\t0.2967\t1.171\t0.852\t74.30\t45.65\t153.9\t0.72\n" +
                   "1000\t0.2734\t1.185\t0.786\t76.80\t47.88\t175.1\t0.74\n")
    DataSet1.close()

    DataSet2 = open(os.path.join(dataDirectory, "DataSet2.txt"), "w", encoding="utf-8", errors="ignore")
    DataSet2.write("""t (°C)\n-180\n-160\n-140\n-120\n-100\n-80\n-60\n-40\n-20\n0\n20\n40\n60\n80\n100\n120\n140\n160
    \n180\n200\n250\n300\n350\n400\n450\n500\n600\n700\n800\n900\n1000\n\nρ (
    kg/m³)\n3.8515\n3.1258\n2.6391\n2.2867\n2.0186\n1.8073\n1.6364\n1.4952\n1.3765\n1.2754\n1.1881\n1.1120\n1.0452\n0
    .9859\n0.9329\n0.8854\n0.8425\n0.8036\n0.7681\n0.7356\n0.6653\n0.6072\n0.5585\n0.5170\n0.4813\n0.4502\n0.3968\n0
    .3577\n0.3243\n0.2967\n0.2734\n\ncp (kJ/(kg*K))\n1.071\n1.036\n1.021\n1.014\n1.011\n1.009\n1.007\n1.006\n1.006\n1
    .006\n1.007\n1.008\n1.009\n1.010\n1.012\n1.014\n1.017\n1.020\n1.023\n1.026\n1.035\n1.046\n1.057\n1.069\n1.081\n1
    .093\n1.116\n1.137\n1.155\n1.171\n1.185\n\n10³ * β (
    1/K)\n11.071\n9.320\n7.758\n6.659\n5.846\n5.219\n4.719\n4.304\n3.962\n3.671\n3.419\n3.200\n3.007\n2.836\n2.684\n2
    .547\n2.423\n2.311\n2.209\n2.115\n1.912\n1.745\n1.605\n1.485\n1.383\n1.293\n1.145\n1.027\n0.932\n0.852\n0.786\n
    \n10³ * λ (W/(m*K))\n9.00\n10.70\n12.70\n14.60\n16.40\n18.16\n19.83\n21.45\n23.01\n24.54\n26.03\n27.49\n28.94\n30
    .38\n31.81\n33.22\n34.66\n36.07\n37.49\n38.91\n42.43\n45.91\n49.31\n52.57\n55.64\n58.48\n63.50\n67.80\n71.30\n74
    .30\n76.80\n\n10⁶ * μ (Pa*s)\n6.44\n7.58\n9.20\n10.49\n11.76\n12.89\n14.02\n15.09\n16.15\n17.10\n17.98\n18.81\n19
    .73\n20.73\n21.60\n22.43\n23.19\n24.01\n24.91\n25.70\n27.40\n29.20\n30.90\n32.55\n34.00\n35.50\n38.30\n40.87\n43
    .32\n45.65\n47.88\n\n10⁶ * ν (m²/s)\n1.67\n2.51\n3.48\n4.587\n5.806\n7.132\n8.567\n10.09\n11.73\n13.41\n15.13\n16
    .92\n18.88\n21.30\n23.15\n25.33\n27.53\n29.88\n32.43\n34.94\n41.18\n48.09\n55.33\n62.95\n70.64\n78.86\n96.08\n114
    .3\n133.6\n153.9\n175.1\n\nPr\n0.77\n0.75\n0.74\n0.73\n0.72\n0.72\n0.71\n0.71\n0.71\n0.70\n0.70\n0.69\n0.69\n0.69
    \n0.69\n0.68\n0.68\n0.68\n0.68\n0.68\n0.67\n0.67\n0.66\n0.66\n0.66\n0.66\n0.67\n0.69\n0.70\n0.72\n0.74""")
    DataSet2.close()


if __name__ == '__main__':
    main()

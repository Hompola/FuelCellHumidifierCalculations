import sys

w = "\033[0m"  # white color
g = "\033[0;32m"  # green color

def WTT(message):  # Write to terminal
    sys.stdout.write(message)


def DisplayData(dataArray):
    # Parse the array for columnWidth
    columnWidth = 0
    for i in range(len(dataArray)):
        for j in range(len(dataArray[i])):
            if len(str(dataArray[i][j])) > columnWidth:
                columnWidth = len(str(dataArray[i][j]))

    # Display the data
    WTT("_" * columnWidth * (len(dataArray[0]) + 1) + "\n")
    for i in range(len(dataArray[0])):
        ID = str(i)
        WTT("ID: " + ID + " " * (columnWidth - 4 - len(ID)) + " | ")
    WTT("\n")
    for line in dataArray:
        for value in line:
            value=str(value)
            WTT(value + " " * (columnWidth - len(value)) + " | ")
        WTT("\n")
    for value in dataArray[0]:
        WTT(value + " " * (columnWidth - len(value)) + " | ")
    WTT("\n")
    for i in range(len(dataArray[0])):
        ID = str(i)
        WTT("ID: " + ID + " " * (columnWidth - 4 - len(ID)) + " | ")
    WTT("\n")
    WTT("_" * columnWidth * (len(dataArray[0]) + 1) + "\n")
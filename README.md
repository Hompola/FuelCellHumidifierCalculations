Fuel Cell Humidifier Calculations
-
These files were created as a part of my bachelor's degree final project.
The code is written entirely in python, and the goal was to make it functional enough to use instead of focusing on aesthetics or semantic quality.
In that, I would recommend using the code for its practical application. These are in detail:
<ol>
  <li>OptimalWaterVolumeCalculations.py module
    <ol>
      <li>TLDR: Determines the number of water molecules consumed by the fuel cell at the given working parameters</li>
      <li>The theory is based on the book Fuel Fundamentals by Ryan O'Hayre et al. and therefore contains implementations of the book's formulas. Used formulas are marked by formula number, chapter and page numbers.</li>
      <li> Contains the necessary functions for calculating the values for a graph comparing the water content of the input air and input hydrogen to the number of water molecules dragged through the membrane. Furthermore, the effect of the inputs' water contents on the membrane's ohmic resistance can be calculated with the results. </li>
    </ol>
  </li>
  
  <li>TemperatureCalculations.py module
    <ol>
      <li>TLDR: Determines the heat loss experienced by the heating fluid in the humidifier's mantle</li>
      <li>The theory is based on classic thermodynamic formulas and technical examples. Relies heavily on text provided by BME EGR as a part of the BME Mechanical Engineering BSc</li>
      <li>Contains formulas for calculating the complete thermodynamic state of a bubble column humidifier. The heat flow through the system, the heat resistance of the components and the temperature of the outer insulation layer. Additionally, it contains formulas for calculating the critical thickness of insulation so that the outer insulation temperature is within a certain threshold above the surrounding air's temperature</li>
    </ol>
  </li>
  
  <li>The main.py module
    <ol>
      <li>TLDR: Uses the above to get and display the results</li>
      <li>Solves for alfa (number of water molecules pulled through the PEM with each hydrogen atom) based on input water contents of hydrogen and air. Displays results in a 3D graph.</li>
      <li>Solves for critical thickness of insulation within a range of 0.1-10 mm at which the insulation temperature is acceptable as described above. Displays the temperatures for every thickness in a 2D plot. </li>
    </ol>
  </li>

<li>Other auxiliary modules
    <ol>
      <li>DisplayData.py and DisplayPlots: Contain independent functions for displaying arrays on terminal and on matplotlib plots.</li>
      <li>LinearInterpolateFromFile.py: uses files in the DataSets folder to interpolate values in a table based on a chosen column and value. Used for example to get the physical properties of water at a given temperature. Able to issue manual and automatic requests.</li>
</ol>

All code is commented and some basic error handling is implemented. Making changes to variables in the files should go without issues. Changing the code requires care, especially when creating datasets! I hope that someone might find this useful and enjoy! 

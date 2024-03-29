# Program Structure

The program has six general sections. One is the atmospheric module, where density, pressure, temperature, and air dynamic viscosity are calculated according to the ISA. The second module takes user inputs in a GUI, such as the Reynolds number, the free stream velocity, the data source, and wing geometric properties. The third module reads airfoil polar data from the corresponding files containing airfoil data, based on the user data source selection. The fourth module is the implementation of the standard LTT, based on the work of Anderson (2017) and Bertin (2021). Outputs of that module contain the lift, drag, and bending moment experienced by the wing. The fifth module is an extended LTT that accounts for non-linearities in the stall region at high angles of attack. The sixth and last module takes care of writing the outputs computed with the classical LTT into an output text file. 

# User Manual

Initialize the model by running the code in any Python environment and a GUI for the necessary inputs will appear (a more user-friendly GUI including a calculation terminal is currently being worked on). Input the wing geometry and initial aerodynamic parameters as well as your desired source for the airfoil data. In the latter you can distinguish between XFoil or Windtunnel data as source, making sure that in the path you give as working directory you have two folders called "XFoil" and "Windtunnel" which contain text files with aerodynamic coefficients in the output format of XFoil, as shown below. 

![grafik](https://github.com/LOMACA/LTT_High_AR/assets/150819500/9c392d1a-979b-4aac-a007-bf10a66404b5)




























These files should be called as follows: the initials "Re" for the selected Reynolds number followed by the integer for the selected Reynolds number. An example would be "Re1200000.txt", indicating a selected Reynolds number of 1'200'000. Some example files containing aerodynamic coefficients for a NACA2412 airfoil created with XFoil are uploaded here and can be used as trial inputs. 

Once the necessary files are in the folder and the GUI has been fed with the required data, click "Submit" on the GUI and close the GUI, after which the program will run and store the wing aerodynamic data in an output file called "output.txt" in the set working directory. Currently, only one angle of attack of the wing can be simulated at a time. Therefore, should the user want a sequence of angles of attack simulated, she or he should re-run the model with a different angle of attack and the new data will be stored in the same output text file under the data of the first run. An extension of the model, where a user can select a range of angles of attack is currently worked on. 

The model also allows interpolation of airfoil data and the user does not necessarily have to select the exact Reynolds number of the stored files. Should a Reynolds number in between two simulated Reynolds numbers be chosen, in the wind tunnel case the model automatically interpolates the data based on the closest two Reynolds number to the one chosen by the user. In the XFoil data selection case, the model calls XFoil and runs an angle of attack sequence based on the Reynolds number selected by the user. 

## Error avoidance

* make sure you run the program with an angle of attack that is featured in the 2-D airfoil data or create a new text file with the data at the desired angle of attack
* make sure to select a Reynolds number either contained in the files or else in between two existing files
* check for the correct file name. Reminder: "Re1200000.txt" is an example for a valid file name

## Sources

- Drela, M. (2000). XFOIL (Version 6.91) [Software]. MIT AeroAstro. https://web.mit.edu/drela/Public/web/xfoil/
- Anderson, J. D. (2017). Fundamentals of aerodynamics (Sixth edition). McGraw-Hill Education. https://www.mheducation.com/highered/product/fundamentals-aerodynamics-anderson/M9781264151929.html
- Bertin, J. J., & Cummings, R. M. (2021, August 12). Aerodynamics for Engineers. Higher Education from Cambridge University Press; Cambridge University Press. https://doi.org/10.1017/9781009105842





# User Manual

This program computes aerodynamic data using the Lifting Linear Theory as well as a Non-Linear Numerical Lifting Line Model. 

Initialize the model by running the code in any Python environment and a GUI for the necessary inputs will appear. Input the wing geometry and initial aerodynamic parameters as well as your desired source for the airfoil data. In the latter you can distinguish between XFoil or Windtunnel data as source, making sure that in the path you give as working directory you have two folders called "XFoil" and "Windtunnel" which contain text files with aerodynamic coefficients in the output format of XFoil, as shown below. 

![grafik](https://github.com/LOMACA/LTT_High_AR/assets/150819500/9c392d1a-979b-4aac-a007-bf10a66404b5)




























These files should be called as follows: the initials "Re" for the selected Reynolds number followed by the actual integer for the selected Reynolds number. An example would be "Re1200000.txt". Some example files containing aerodynamic coefficients for a NACA2412 airfoil created with XFoil are uploaded here and can be used as trial inputs. 

Once the necessary files are in the folder and the GUI has been fed with the required data, click "Submit" on the GUI and close the GUI, after which the program will run and store the wing aerodynamic data in an output file called "output.txt" in the set working directory. Currently, only one angle of attack of the wing can be simulated at a time. Therefore, should the user want a sequence of angles of attack simulated, she or he should re-run the model with a different angle of attack and the new data will be stored in the same output text file under the data of the first run. An extension of the model, where a user can select a range of angles of attack is currently worked on. 

The model also allows interpolation of airfoil data and the user does not necessarily have to select the exact Reynolds number of the stored files. Should a Reynolds number in between two simulated Reynolds numbers be chosen, in the wind tunnel case the model automatically interpolates the data based on the closest two Reynolds number to the one chosen by the user. In the XFoil data selection case, the model calls XFoil and runs an angle of attack sequence based on the Reynolds number selected by the user. 

## Error avoidance

* make sure you run the program with an angle of attack that is featured in the 2-D airfoil data or create a new text file with the data at the desired angle of attack
* make sure to select a Reynolds number either contained in the files or else in between two existing files
* check for the correct file name. Reminder: "Re1200000.txt" is an example for a valid file name

## Sources

- Drela, M. (2000). XFOIL (Version 6.91) [Software]. MIT AeroAstro. https://web.mit.edu/drela/Public/web/xfoil/
- Anderson, J. D. (2017). Fundamentals of aerodynamics (Sixth edition). McGraw-Hill Education.
- Bertin, J. J., & Cummings, R. M. (2021, August 12). Aerodynamics for Engineers. Higher Education from Cambridge University Press; Cambridge University Press. https://doi.org/10.1017/9781009105842





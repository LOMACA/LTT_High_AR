This program computes aerodynamic data using the Lifting Linear Theory as well as a Non-Linear Numerical Lifting Line Model. 


User manual:

Initialize the model by running the code in any Python environment and a GUI for the necessary inputs will appear. Input the wing geometry and initial aerodynamic parameters as well as your desired source for the airfoil data. In the latter you can distinguish between XFoil or Windtunnel data as source, making sure that in the path you give as working directory you have two folders called "XFoil" and "Windtunnel" which contain text files with aerodynamic coefficients in the output format of XFoil, as shown below. 

![grafik](https://github.com/LOMACA/LTT_High_AR/assets/150819500/9c392d1a-979b-4aac-a007-bf10a66404b5)




























These files should be called as follows: the initials "Re" for the selected Reynolds number followed by the actual integer for the selected Reynolds number. An example would be "Re1200000.txt". Some example files containing aerodynamic coefficients for a NACA2412 airfoil created with XFoil are uploaded here and can be used as trial inputs. 

Once the necessary files are in the folder and the GUI has been fed with the required data, click "Submit" on the GUI and close the GUI, after which the program will run and store the wing aerodynamic data in an output file called "output.txt" in the set working directory. Currently, only one angle of attack of the wing can be simulated at a time. Therefore, should the user want a sequence of angles of attack simulated, she or he should re-run the model with a different angle of attack and the new data will be stored in the same output text file under the data of the first run. 


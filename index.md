![surilogo.PNG](attachment:surilogo.PNG)

# Interactive Predictive Storage Capacity Calculator

## To those viewing this webpage, this is a preview of what my jupyter notebook looks like. The interactive part of the calculator is unavaliable in this format as it requires python functionality. Thank you!

### Carlos Figueroa-Diaz, SURI Researcher, University of Texas at Austin

### The Objective

The purpose of this notebook is to provide a calculator that inputs various production data and outputs predicted CO<sub>2 </sub>  storage capacity. This calculator takes in the following inputs:

**Pressure** - The reservoir pressure in psia

**CUMO** - Cumulative Oil Production in mmbbl ($10^6$ barrels)

**CUMG** - Cumulative Gas Production in Bcf (($10^9$ cubic ft)

**Water Volume** - Initial water volume in mmbbl ($10^6$ barrels)

**Injection Volume** - Volume of injected water in mmbbl ($10^6$ barrels)

**Depth** - The resevoir depth in ft

**Geothermal Gradient** - The temperature gradient of the reservoir in °C/km

**Mean Surface Temperature** - The average surface temperature in °C

**API** - The API gravity of the oil

**Efficiency** - The percent storage efficiency 

**GOR** - The gas oil ratio

**Density of CO2** - Density of CO2 in kg/m<sup>3</sup>

###  Obtaining CO2 Density

The code below approximates the density of CO2 using a given temperature and pressure. This is done by using a multivariable regression model of an equation of state for carbon dioxide. The data for this EOS model is provided by the National Institute of Standards and Technology.

### Key Assumptions

This calculator assumes the following:

1. The gas in the reservoir is pure methane.

2. At discovery, the oil is undersaturated and above the bubble point so there is no gas cap. Therefore the GOR reflects the solution gas ratio.

3. The gas deviation factor is 1.

With these assumptions in mind, note that this calculator is best suited for pure gas or pure oil reservoirs. A calculator suited for mixed reservoirs requires more complicated equations that are better suited for interactions between gas and oil.

### The Code
Now that the inputs for this calculator have been defined, let's get into the code. There are sliders for each input. Each slider can be moved to a range of values by the user, to see how changing a certain feature affects predicted storage capacity in real time. These inputs are then used in a function called carbonCalc. 

carbonCalc takes in the inputs from the sliders and through a series of unit conversions and formulas, outputs the predicted storage of CO<sub>2</sub> in kg/m<sup>3</sup>. It calculates the bubble point pressure, formation volume factor for oil, the gas expansion factor, and uses this information to calculate the voided reservoir volume in millions of m<sup>3</sup>. 

The density of CO<sub>2 </sub>, which is approximated using a multivariable regression of a CO2 equation of state, is then multiplied by this reservoir volume to determine the approximate amount of CO<sub>2 </sub> in megatons you could store in the reservoir. Once the inputs are selected, the predicted storage capacity is displayed under the sliders. At any time the user is free to change inputs and observe how a change affects predicted storage capacity. Make sure to run all of the cells below. 


```python
# -*- coding: utf-8 -*-
"""
@author: Carlos Figueroa-Diaz @ The University of Texas at Austin
"""
import numpy as np
from IPython.display import display
from ipywidgets import interactive                      # widgets and interactivity
from ipywidgets import widgets                            
from ipywidgets import Layout
from ipywidgets import Label
from ipywidgets import VBox, HBox
from ipywidgets import BoundedFloatText

name = 'Predicted Storage Capacity Calculator, Carlos Figueroa-Diaz,  The University of Texas at Austin'

#Making an author label
l = widgets.Text(value=name.center(200),layout=Layout(width='950px', height='30px'))

#Creating sliders for all of the CO2 params
P = widgets.FloatSlider(min = 500, max = 16000, value = 2000 , step = 1, description = 'Pressure (psia)',
                        orientation='vertical', continuous_update=True, layout=Layout(width='120px', height='300px'))

CUMO = ndata = widgets.FloatSlider(min = 0, max = 1e2, value = 1, step = 1, description = 'CUMO ($10^{6}$ bbl)',
                        orientation='vertical', continuous_update=True, layout=Layout(width='120px', height='300px'))

CUMG = ndata = widgets.FloatSlider(min = 0, max = 5e3, value = 1, step = 1, description = 'CUMG ($10^{9}$ ft$^3$)',
                        orientation='vertical', continuous_update=True, layout=Layout(width='120px', height='300px'))

water = ndata = widgets.FloatSlider(min = 0, max = 1e4, value = 0, step = 1, description = 'Water Volume ($10^{6}$ bbl)',
                        orientation='vertical', continuous_update=True, layout=Layout(width='140px', height='300px'))

inject = ndata = widgets.FloatSlider(min = 0, max = 1e4, value = 0, step = 1, description = 'Injection Volume ($10^{6}$ bbl)',
                        orientation='vertical', continuous_update=True, layout=Layout(width='160px', height='300px'))

depth_ft = ndata = widgets.FloatSlider(min = 1500, max = 17000, value = 2000, step = 1, description = 'Depth (ft)',
                        orientation='vertical', continuous_update=True, layout=Layout(width='120px', height='300px'))

geoGrad = ndata = widgets.FloatSlider(min = 14, max = 52, value = 25, step = 1, description = 'Geothermal Gradient (C/km)',
                        orientation='vertical',continuous_update=True, layout=Layout(width='180px', height='300px'))

meanST = ndata = widgets.FloatSlider(min = -20, max = 50, value = 20, step = 1, description = 'Mean Surface Temp (C)',
                        orientation='vertical',continuous_update=True, layout=Layout(width='140px', height='300px'))

API = ndata = widgets.FloatSlider(min = 0, max = 80, value = 45, step = 1, description = 'API', orientation='vertical',
                        continuous_update=True, layout=Layout(width='120px', height='300px'))

efficiency = ndata = widgets.FloatSlider(min = .01, max = 1, value = .75, step = .01, description = 'Efficiency',
                        orientation='vertical',continuous_update=True, layout=Layout(width='120px', height='300px'))

#initializes the user interface
uipars = widgets.HBox([P,CUMO, CUMG, water, inject, depth_ft, geoGrad, meanST, API, efficiency])
uik = widgets.VBox([l,uipars])

#A quintic multivariable regression model for an EOS of CO2 density with respect to temperature and pressure
def EOS_quintic(X,Y):
    #convert pressure to MPa
    Y /= 145.038
    var = np.array([1, X, Y, X**2, X*Y, Y**2, X**3, X**2 * Y, X * Y**2, Y**3, X**4, X**3 * Y, X**2 * Y**2, X * Y**3,
                    Y**4, X**5, X**4 * Y, X**3 * Y**2, X**2 * Y**3, X * Y**4, Y**5])
    coeff = np.array([8.88392845e+02, -1.95891892e+01,  3.32715613e+01,  4.72967153e-02, 7.35818191e-01, -1.43785734e+00,
                      5.33269784e-04, -5.13514242e-03, -5.55746959e-03,  2.32765649e-02, -2.69521096e-06,  4.36857231e-06,
                      6.47348502e-05, -3.18001817e-05, -1.51807659e-04,  3.46984410e-09, 8.51630513e-09, -6.49695396e-08,
                      -1.98531352e-07,  3.12705430e-07, 3.26340650e-07])
    return np.dot(var, coeff)

#prints the storage capacity of CO2 based on field data
def carbonCalc(P, CUMG, CUMO, water, inject, depth_ft, geoGrad, meanST, API, efficiency):
    #converts production volumes to m3
    oProd = CUMO * 0.16
    gProd = CUMG *0.0283*1000
    wProd = water * 0.16
    iProd = inject * 0.16
    #preliminary unit conversions and calculations 
    if CUMO > 0:
        GOR = CUMG / CUMO * 1000
    else:
        GOR = 0
    gasSC = 0.0991
    virginP = P / 145.0
    depth_m = depth_ft * 0.3048
    resT = meanST +  geoGrad * depth_m / 1000.0
    resTR = resT * 9.0/5.0 + 492
    #hydrostatic pressure has yet to be implemented in this code
    hydroP = depth_m * 10.52 / 1000.0
    oSG = 141.5/ (API + 131)  
    cO = (55.233 * 1e-6) - ((60.588 * 1e-6) * oSG)
    alpha = 0.00091 * (resTR -460) - 0.0125 * API
    #calculates the bubblepoint pressure
    pBub = 18.2 *(((GOR/gasSC)**.83)*(10 ** alpha) -1.4)
    #calculates the formation volume factor for oil and the gas expansion factor
    gasExpFac = 0.101 * (resT + 273) * 1/ (virginP * 273)
    bubFac = 1.0113 + 7.2046e-5 * ((((GOR**0.3738)*((gasSC ** 0.2914) / (oSG ** 0.6265))) + 0.24626 * ((resTR-460)**0.5371)) ** 3.0936)
    volFac = bubFac ** (-cO *(P - pBub))
    #calculates total volume of production
    resVol = oProd *volFac + gProd * gasExpFac + wProd - iProd
    #converts this to megatons of CO2
    rhoCO2 = EOS_quintic(resT, P)
    capacity = resVol * efficiency * rhoCO2 / 1000
    print("Density of CO2 = " + str(round(rhoCO2, 3)) + " kg/m3")
    print("Voided Reservoir Volume = " + str(round(resVol, 3)) +  " million m3")
    print("Storage Capacity = " + str(round(capacity, 3)) + " Megatons of CO2")
    
       
#Plots all of the sliders for the parameters, allowing a user to change parameters in real time
interactive_plot = widgets.interactive_output(carbonCalc, {'P' : P ,'CUMG' : CUMG , 'CUMO' : CUMO, 'water' : water, 'inject' : inject,'depth_ft' : depth_ft
                                                           , 'geoGrad' : geoGrad , 'meanST' : meanST,'API' : API,'efficiency': efficiency})


```

### Give it a try!

Run the cell below. Using production data from, change the sliders accordingly. For exact inputs, I recommend double clicking the number slots to type your values instead of using the slider, but it is up to you. As you change inputs, you'll see the predicted storage capacity update in real time.

![surilogo.PNG](attachment:surilogo.PNG)

### References

Thermophysical properties of CO<sub>2 </sub>: https://webbook.nist.gov/cgi/fluid.cgi?ID=C124389&Action=Page

Equations for formation volume factor, bubble point pressure, and oil specific gravity: https://www-sciencedirect-com.ezproxy.lib.utexas.edu/science/article/pii/B9780128002193000048

Oil Compressibility Equation: http://dx.doi.org/10.2118/92-03-02



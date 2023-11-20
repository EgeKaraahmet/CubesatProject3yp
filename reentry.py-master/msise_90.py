
import numpy as np
from io import StringIO
from bisect import bisect_left

class database:
    # data
    __atm_data = """
    MSISE-90 - Mean Solar Activity
    Altitude(km)|Temp(K)|Density(kg/m3)|Pressure(Pa)|Mol. Wt.(kg/kmol)
    0|300.2511|1.17E+00|1.01E+05|28.9502
    20|206.2085|9.49E-02|5.62E+03|28.9502
    40|257.6979|4.07E-03|3.02E+02|28.9502
    60|244.1212|3.31E-04|2.32E+01|28.9502
    80|196.3636|1.68E-05|9.45E-01|29.0175
    100|184.0160|5.08E-07|2.81E-02|27.7137
    120|374.9715|1.80E-08|2.17E-03|25.8745
    140|635.5703|3.26E-09|7.03E-04|24.5349
    160|787.5532|1.18E-09|3.31E-04|23.4225
    180|877.6729|5.51E-10|1.80E-04|22.4106
    200|931.2806|2.91E-10|1.05E-04|21.4734
    220|963.2701|1.66E-10|6.44E-05|20.6108
    240|982.4191|9.91E-11|4.09E-05|19.8292
    260|993.9173|6.16E-11|2.66E-05|19.1337
    280|1000.8427|3.94E-11|1.77E-05|18.5256
    300|1005.0267|2.58E-11|1.20E-05|18.0015
    320|1007.5620|1.72E-11|8.20E-06|17.5537
    340|1009.1030|1.16E-11|5.69E-06|17.1721
    360|1010.0423|7.99E-12|3.98E-06|16.8449
    380|1010.6166|5.55E-12|2.81E-06|16.5597
    400|1010.9688|3.89E-12|2.01E-06|16.3044
    420|1011.1853|2.75E-12|1.44E-06|16.0669
    440|1011.3190|1.96E-12|1.04E-06|15.8360
    460|1011.4014|1.40E-12|7.55E-07|15.6008
    480|1011.4526|1.01E-12|5.53E-07|15.3508
    500|1011.4845|7.30E-13|4.07E-07|15.0760
    520|1011.5043|5.31E-13|3.03E-07|14.7669
    540|1011.5168|3.88E-13|2.27E-07|14.4148
    560|1011.5245|2.85E-13|1.71E-07|14.0125
    580|1011.5294|2.11E-13|1.31E-07|13.5547
    600|1011.5325|1.56E-13|1.01E-07|13.0389
    620|1011.5345|1.17E-13|7.89E-08|12.4665
    640|1011.5357|8.79E-14|6.24E-08|11.8428
    660|1011.5365|6.65E-14|5.01E-08|11.1779
    680|1011.5370|5.08E-14|4.07E-08|10.4854
    700|1011.5374|3.91E-14|3.36E-08|9.7818
    720|1011.5375|3.04E-14|2.82E-08|9.0847
    740|1011.5377|2.39E-14|2.39E-08|8.4111
    760|1011.5377|1.90E-14|2.06E-08|7.7753
    780|1011.5378|1.53E-14|1.79E-08|7.1884
    800|1011.5378|1.25E-14|1.58E-08|6.6572
    820|1011.5378|1.03E-14|1.40E-08|6.1849
    840|1011.5379|8.64E-15|1.26E-08|5.7711
    860|1011.5379|7.32E-15|1.14E-08|5.4132
    880|1011.5379|6.28E-15|1.04E-08|5.1066
    900|1011.5379|5.46E-15|9.47E-09|4.8460
    """
 
    # read parse and store data in ndr object
    def __init__(self):
        self.__ndr = self.__read_data()
        
    # parse atmospheric data (from data string)
    def __read_data(self):
        atw = np.genfromtxt(StringIO(self.__atm_data), delimiter='|', skip_header=2)
        return atw

    # simple linear interpolation
    def __interpolate(self, x_list, y_list, x):
        if any(y - x <= 0 for x, y in zip(x_list, x_list[1:])):
            raise ValueError("x_list must be in strictly ascending order!")
        intervals = zip(x_list, x_list[1:], y_list, y_list[1:])
        slopes = [(y2 - y1) / (x2 - x1) for x1, x2, y1, y2 in intervals]

        if x <= x_list[0]:
            return y_list[0]
        elif x >= x_list[-1]:
            return y_list[-1]
        else:
            i = bisect_left(x_list, x) - 1
            return y_list[i] + slopes[i] * (x - x_list[i])

    # returns interpolated atmos data for each alt (altitude value)
    def get_atmospheric_data(self, alt):
        # split the arrays
        altitude = self.__ndr[:, 0]
        temperature = self.__ndr[:, 1]
        density = self.__ndr[:, 2]
        pressure = self.__ndr[:, 3]
        molwt = self.__ndr[:, 4]
        # interpolate values
        t = self.__interpolate(altitude, temperature, alt)
        d = self.__interpolate(altitude, density, alt)
        p = self.__interpolate(altitude, pressure, alt)
        m = self.__interpolate(altitude, molwt, alt)

        # Compute dynamic viscosity using Sutherland's formula
        # https://doc.comsol.com/5.5/doc/com.comsol.help.cfd/cfd_ug_fluidflow_high_mach.08.27.html
        # miu_oxygen = 1.919*1e-5
        # miu = miu_oxygen *((273+139)/(t+139))*(t/273)**(1.5)
        # Compute mean free path
        # l = np.sqrt(np.pi/8) * (0.498/miu) * (d*p)**(-1/2)
        # particle diameter --> Oxygen atoms and nitrogen atoms--> consider the mean of the two
        # dp = 0.5 *(0.74+0.71)*1e-10
        # Compute Knudsen number (Kn)
        # https://doc.comsol.com/5.6/doc/com.comsol.help.particle/particle_ug_fluid_flow.08.37.html
        # Kn = l/dp

        # Boltzmann constant --> assume at high altitude, boltzmann gas
        kb = 1.380649*1e-23
        # particle diameter --> Oxygen atoms and nitrogen atoms--> consider the mean of the two
        dp = (2.94+3)/2*1e-10

        # characteristic length of our cubesat
        L = 0.3 # 0.3

        # Kundsen number
        Kn = kb*t / (np.sqrt(2)*L*np.pi*p*dp**2)

        # Correction factor:
        S = 1 + Kn * (2.514 + 0.8*np.exp(-0.55/Kn))

        # return t, d, p, m,Kn, S
        return d



# main function
if __name__ == "__main__":
    # create instance
    atm = database()
    
    # get data for sample altitude (given in km)
    alt = 200  # km

    # t, d, p, m,Kn, S= atm.get_atmospheric_data(alt)
    d = atm.get_atmospheric_data(alt)
    
    # print atm data for that altitude
    from tabulate import tabulate
    data = []
    # headers=['Altitude(km)', 'Temp(K)', 'Density(kg/m3)', 'Pressure(Pa)', 'Mol. Wt.(kg/kmol)', 'Kundsen no.', 'Slip factor']
    # headers = [ 'Kundsen no.','Slip factor']
    headers = ['density']
    # data.append( [alt, t, d, p, m, Kn, S] )
    data.append([d])
    print ( tabulate(data, headers) )
   


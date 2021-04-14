Hi Jinkai, 

This folder probably contains almost of the information you need to verify the scatter kernel. 


>>The velocities of all the molecules are initialised by the Box-Muller transform so they will follow a Maxwellian distribution. 

>>The temperature of the wall are set as 423K. 


>>Take reflected_Tx_0.1.txt as example, this file contains the reflected information for the molecules with initial velocity around 0.1*Cm in tangential x component, where Cm is the most probable speed Cm = (2 kB Tw / m). 

>> reflected_normal_19.txt contains the reflected information for the molecules with initial velocity in normal component. 


>> You might find the accommodation coefficients are essential in the scattering kernel. Here are some example for it. 
>> Tangential_AC.txt: the first column is the non-dimensionlised velocity, the second colume is the tangential accommodation coefficients. 



ySI=2*pi*ThisField.CentralFrequency/1.60218e-19*sqrt(m*1e-3*MediumParameters.IonizationEnergy*3e8*epsilon0*Element.RefractiveIndex/(2*max(max(I))))

yAI=(2*pi*ThisField.CentralFrequency/4.13e16)*sqrt(0.5*MediumParameters.IonizationEnergy*6.242e+18/27.21)/(max(max(abs(ThisField.Data)))/5.1422e11)


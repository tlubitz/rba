#import pandas
import matplotlib.pyplot as plt
import numpy
from rbatools.rba_Session import RBA_Session

dirInput2 = 'IND/iMM904_RBA/iMM904_RBA_DefaultKappEstim_Lahtvee'
targetMu = 0.1
Simulation = RBA_Session(dirInput2)
orig_enz = Simulation.model.parameters.functions._elements_by_id[
    'default_efficiency'].parameters._elements_by_id['CONSTANT'].value
kapp_LB = 0
kapp_UB = orig_enz*20
new_kapp = orig_enz
Mu_pred = Simulation.findMaxGrowthRate()
Mus = []
Kapps = []

while abs(targetMu - Mu_pred) > 0.01:
    print(new_kapp, abs(targetMu - Mu_pred))
    Simulation.model.parameters.functions._elements_by_id[
        'default_efficiency'].parameters._elements_by_id['CONSTANT'].value = new_kapp
    Simulation.model.parameters.functions._elements_by_id[
        'default_transporter_efficiency'].parameters._elements_by_id['CONSTANT'].value = new_kapp
    Simulation.rebuild_from_model()
    Mu_pred = Simulation.findMaxGrowthRate()
    Mus.append(abs(targetMu - Mu_pred))
    Kapps.append(new_kapp)
    if Mu_pred > targetMu:
        new_kapp_prelim = kapp_LB+(0.5*abs(kapp_LB-new_kapp))
    elif Mu_pred < targetMu:
        new_kapp_prelim = kapp_UB-(0.5*abs(new_kapp-kapp_UB))
    new_kapp = new_kapp_prelim

plt.scatter(Kapps, Mus)
plt.show()

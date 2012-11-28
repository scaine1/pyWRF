#module: pertibation_variables

pert_variable_dict={}
pert_variable_dict.update({'GEOP':('PH','PHB')})
pert_variable_dict.update({'Z':('PH','PHB')})

#pert_variable_dict.update({'T':('T',300.)})
pert_variable_dict.update({'THETA':('T',300.)})

#pert_variable_dict.update({'P':('P','PB')})
pert_variable_dict.update({'PRES':('P','PB')})



calc_variable_dict={}
#calc_variable_dict.update({'TK':('PRES','T')})
calc_variable_dict.update({'TEMP':('PRES','THETA')})

calc_variable_dict.update({'RH':('QVAPOR','PRES','TEMP')})

calc_variable_dict.update({'TD':('QVAPOR','PRES')})

calc_variable_dict.update({'PRES':''})		#used to divide by 100. to get into mb

calc_variable_dict.update({'Z':''})		#

calc_variable_dict.update({'SPH':('QVAPOR',)})

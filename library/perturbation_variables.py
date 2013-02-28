#module: pertibation_variables


#####################FOR THE ARW CORE##############
#####################FOR THE ARW CORE##############
#####################FOR THE ARW CORE##############
pert_variable_dict_ARW={}
pert_variable_dict_ARW.update({'GEOP':('PH','PHB')})
pert_variable_dict_ARW.update({'Z':('PH','PHB')})
pert_variable_dict_ARW.update({'THETA':('T',300.)})
pert_variable_dict_ARW.update({'PRES':('P','PB')})

calc_variable_dict_ARW={}
calc_variable_dict_ARW.update({'TEMP':('PRES','THETA')})
calc_variable_dict_ARW.update({'RH':('QVAPOR','PRES','TEMP')})
calc_variable_dict_ARW.update({'TD':('QVAPOR','PRES')})
calc_variable_dict_ARW.update({'PRES':''})		#used to divide by 100. to get into mb
calc_variable_dict_ARW.update({'Z':''})		#
calc_variable_dict_ARW.update({'SPH':('QVAPOR',)})


#####################FOR THE NMM CORE##############
#####################FOR THE NMM CORE##############
#####################FOR THE NMM CORE##############

pert_variable_dict_NMM={}
calc_variable_dict_NMM={}
calc_variable_dict_NMM.update({'PRES':('PINT',)})		#used to divide by 100. to get into mb
calc_variable_dict_NMM.update({'TEMP':('T',)})		#used to divide by 100. to get into mb
calc_variable_dict_NMM.update({'SPH':('QVAPOR',)})
calc_variable_dict_NMM.update({'RH':('QVAPOR','PRES','TEMP')})
calc_variable_dict_NMM.update({'TD':('QVAPOR','PRES')})
calc_variable_dict_NMM.update({'Z':('',)})
calc_variable_dict_NMM.update({'VTMK':('TEMP','QVAPOR')})
calc_variable_dict_NMM.update({'XLAT':('GLAT',)})
calc_variable_dict_NMM.update({'XLONG':('GLON',)})
#calc_variable_dict_NMM.update=({'XLAT':'GLAT'})
#calc_variable_dict_NMM.update=({'XLONG':'GLON'})
#calc_variable_dict_NMM.update({'PRES':''})		#used to divide by 100. to get into mb

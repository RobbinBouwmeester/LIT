def create_acyls(carbon,db,link="A",
				 c_mass=12.000000,
				 h_mass=1.007825037,
				 o_mass=15.99491464):
	chem_form = "C"+str(carbon)+"H"+str(((carbon*2)-(2*db)))+"O"+str(2)

	c_tot_mass = carbon*c_mass
	h_tot_mass = ((carbon*2)-(2*db))*h_mass
	o_tot_mass = o_mass*2

	if link == "O":
		chem_form = "C"+str(carbon)+"H"+str(((carbon*2)-(2*db))+2)+"O"+str(1)
		h_tot_mass = h_tot_mass + 2*h_mass
		o_tot_mass -= o_mass
	if link == "P":
		chem_form = "C"+str(carbon)+"H"+str(((carbon*2)-(2*db)))+"O"+str(1)
		o_tot_mass -= o_mass
	tot_mass = sum([c_tot_mass,h_tot_mass,o_tot_mass])
	deproto_mass = tot_mass-h_mass
	water_mass = tot_mass-(h_mass*2)-o_mass
	deproto_water_mass = tot_mass-(h_mass*3)-o_mass

	return(tot_mass,deproto_mass,water_mass,deproto_water_mass,chem_form)

def get_acyl_dict(min_acyl=12,max_axyl=25,min_db=0,max_db=6,link_types=["A","P","O"]):
	fa_mw_dict = {}
	fa_form_dict = {}

	for i in range(12,25):
		for j in range(0,6):
			for l in link_types:
				tot_mass,deproto_mass,water_mass,deproto_water_mass,chem_form = create_acyls(i,j,link=l)
				if l == "A": 
					fa_mw_dict[str(i)+":"+str(j)+"|[M]"] = tot_mass
					fa_mw_dict[str(i)+":"+str(j)+"|[M-H2O]"] = water_mass
					fa_mw_dict[str(i)+":"+str(j)+"|[M-H]-"] = deproto_mass
					fa_mw_dict[str(i)+":"+str(j)+"|[M-H2O-H]-"] = deproto_water_mass
					fa_form_dict[str(i)+":"+str(j)] = chem_form
				else: 
					fa_mw_dict[l+"-"+str(i)+":"+str(j)+"|[M-H]-"] = deproto_mass
					fa_mw_dict[l+"-"+str(i)+":"+str(j)+"|[M-H2O-H]-"] = deproto_water_mass
					fa_form_dict[l+"-"+str(i)+":"+str(j)] = chem_form
	return(fa_mw_dict,fa_form_dict)

if __name__ == "__main__":
	print(get_acyl_dict())
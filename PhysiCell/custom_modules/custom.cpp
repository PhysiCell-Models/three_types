/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2018, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/

#include "./custom.h"

void create_cell_types( void )
{
	// set the random seed 
	SeedRandom( parameters.ints("random_seed") );  
	
	/* 
	   Put any modifications to default cell definition here if you 
	   want to have "inherited" by other cell types. 
	   
	   This is a good place to set default functions. 
	*/ 
	
	cell_defaults.functions.volume_update_function = standard_volume_update_function;
	cell_defaults.functions.update_velocity = standard_update_cell_velocity;

	cell_defaults.functions.update_migration_bias = NULL; 
	cell_defaults.functions.update_phenotype = NULL; // update_cell_and_death_parameters_O2_based; 
	cell_defaults.functions.custom_cell_rule = NULL; 
	
	cell_defaults.functions.add_cell_basement_membrane_interactions = NULL; 
	cell_defaults.functions.calculate_distance_to_membrane = NULL; 
	
	/*
	   This parses the cell definitions in the XML config file. 
	*/
	
	initialize_cell_definitions_from_pugixml(); 
	
	/* 
	   Put any modifications to individual cell definitions here. 
	   
	   This is a good place to set custom functions. 
	*/ 

	get_cell_definition("A").functions.update_phenotype = A_phenotype; 
	get_cell_definition("B").functions.update_phenotype = B_phenotype; 
	get_cell_definition("C").functions.update_phenotype = C_phenotype; 
		
	/*
	   This builds the map of cell definitions and summarizes the setup. 
	*/
		
	build_cell_definitions_maps(); 
	display_cell_definitions( std::cout ); 
	
	return; 
}

void setup_microenvironment( void )
{
	// set domain parameters 
	
	// put any custom code to set non-homogeneous initial conditions or 
	// extra Dirichlet nodes here. 
	
	// initialize BioFVM 
	
	initialize_microenvironment(); 	
	
	return; 
}

void setup_tissue( void )
{
	double Xmin = microenvironment.mesh.bounding_box[0]; 
	double Ymin = microenvironment.mesh.bounding_box[1]; 
	double Zmin = microenvironment.mesh.bounding_box[2]; 

	double Xmax = microenvironment.mesh.bounding_box[3]; 
	double Ymax = microenvironment.mesh.bounding_box[4]; 
	double Zmax = microenvironment.mesh.bounding_box[5]; 
	
	if( default_microenvironment_options.simulate_2D == true )
	{
		Zmin = 0.0; 
		Zmax = 0.0; 
	}
	
	double Xrange = Xmax - Xmin; 
	Xmin += 0.1*Xrange; 
	Xrange *= 0.8;
	double Yrange = Ymax - Ymin; 
	Ymin += 0.1*Yrange; 
	Yrange *= 0.8;
	double Zrange = Zmax - Zmin; 
	Zmin += 0.1*Zrange; 
	Zrange *= 0.8; 
	// create some of each type of cell 
	
	Cell* pC;
	
	// place A
	
	for( int n = 0 ; n < parameters.ints("number_of_A") ; n++ )
	{
		std::vector<double> position = {0,0,0}; 
		position[0] = Xmin + UniformRandom()*Xrange; 
		position[1] = Ymin + UniformRandom()*Yrange; 
		position[2] = Zmin + UniformRandom()*Zrange; 
		
		pC = create_cell( get_cell_definition("A") ); 
		pC->assign_position( position );
		for( int k=0 ; k < pC->phenotype.death.rates.size() ; k++ )
		{ pC->phenotype.death.rates[k] = 0.0; }
	}
	
	// place B
	
	for( int n = 0 ; n < parameters.ints("number_of_B") ; n++ )
	{
		std::vector<double> position = {0,0,0}; 
		position[0] = Xmin + UniformRandom()*Xrange; 
		position[1] = Ymin + UniformRandom()*Yrange; 
		position[2] = Zmin + UniformRandom()*Zrange; 
		
		pC = create_cell( get_cell_definition("B") ); 
		pC->assign_position( position );
		for( int k=0 ; k < pC->phenotype.death.rates.size() ; k++ )
		{ pC->phenotype.death.rates[k] = 0.0; }
	}

	// place C
	
	for( int n = 0 ; n < parameters.ints("number_of_C") ; n++ )
	{
		std::vector<double> position = {0,0,0}; 
		position[0] = Xmin + UniformRandom()*Xrange; 
		position[1] = Ymin + UniformRandom()*Yrange; 
		position[2] = Zmin + UniformRandom()*Zrange; 
		
		pC = create_cell( get_cell_definition("C") ); 
		pC->assign_position( position );
		for( int k=0 ; k < pC->phenotype.death.rates.size() ; k++ )
		{ pC->phenotype.death.rates[k] = 0.0; }
	}

	return; 
}

std::vector<std::string> my_coloring_function( Cell* pCell )
{
	static int A_type = get_cell_definition( "A" ).type; 
	static int B_type = get_cell_definition( "B" ).type; 
	static int C_type = get_cell_definition( "C" ).type; 
	
	// start with flow cytometry coloring 
	
	std::vector<std::string> output = {"black" , "black" , "black" , "black"} ;

	// color live C 
		
	if( pCell->phenotype.death.dead == false && pCell->type == A_type )
	{
		 output[0] = parameters.strings("A_color");  
		 output[2] = parameters.strings("A_color");  
	}
	
	// color live B

	if( pCell->phenotype.death.dead == false && pCell->type == B_type )
	{
		 output[0] = parameters.strings("B_color");  
		 output[2] = parameters.strings("B_color");  
	}
	
	// color live C

	if( pCell->phenotype.death.dead == false && pCell->type == C_type )
	{
		 output[0] = parameters.strings("C_color");  
		 output[2] = parameters.strings("C_color");  
	}

	return output; 
}


void predator_hunting_function( Cell* pCell, Phenotype& phenotype, double dt )
{
	Cell* pTestCell = NULL; 
	
	double sated_volume = pCell->parameters.pReference_live_phenotype->volume.total * 
		parameters.doubles("relative_sated_volume" ); 
	
	for( int n=0; n < pCell->cells_in_my_container().size() ; n++ )
	{
		pTestCell = pCell->cells_in_my_container()[n]; 
		// if it's not me, not dead, and not my type, eat it 
		
		if( pTestCell != pCell && pTestCell->type != pCell->type && pTestCell->phenotype.death.dead == false )
		{
			// only eat if I'm not full 
			if( phenotype.volume.total < sated_volume )
			{
				pCell->ingest_cell(pTestCell); 
				return; 
			}
	
		}
	}
	
	return; 
}

void predator_cycling_function( Cell* pCell, Phenotype& phenotype, double dt )
{
	double sated_volume = pCell->parameters.pReference_live_phenotype->volume.total * 
		parameters.doubles("relative_sated_volume" ); 
	
	if( phenotype.volume.total > sated_volume )
	{ phenotype.cycle.data.transition_rate(0,1) = get_cell_definition("prey").phenotype.cycle.data.transition_rate(0,1) * 0.01; }
	else
	{ phenotype.cycle.data.transition_rate(0,1) = 0; }
	return; 
}

void prey_cycling_function( Cell* pCell , Phenotype& phenotype, double dt )
{
	static int signal_index = microenvironment.find_density_index( "prey signal" ); 
	
	double threshold = parameters.doubles("prey_quorom_threshold" ) + 1e-16 ; 
	double factor = (threshold - pCell->nearest_density_vector()[signal_index] )/threshold; 
	if( factor < 0 )
	{ factor = 0.0; } 
	
	phenotype.cycle.data.transition_rate(0,1) = get_cell_definition("prey").phenotype.cycle.data.transition_rate(0,1); 
	phenotype.cycle.data.transition_rate(0,1) *= factor; 
	
	return; 
}

/*
up_down_signal add_effect( up_down_signal start_signal , double factor, char factor_type )
{
	up_down_signal out = start_signal; 
	// neutral signal 
	if( factor_type == 'N' || factor_type == 'n' )
	{ return out; }

	// promoter signal 
	if( factor_type == 'P' || factor_type == 'p' )
	{
		// up = sum of all (scaled) promoter signals 
		out.up += factor; 
		// cap up signal at 1.0 
		if( out.up > 1.0 )
		{ out.up = 1.0; }
		return out; 
	}

	// inhibitor signal 
	// down is max of all (scaled) inhibitor signals 
	if( factor > out.down )
	{ out.down = factor; }
	if( out.down > 1 )
	{ out.down = 1.0; }
	return out; 
}
*/

up_down_signal::up_down_signal()
{
    up = 0.0;
	down = 0.0; 
    no_promoters = true; 
    no_inhibitors = true; 
	return; 
} 

void up_down_signal::add_effect( double factor, char factor_type )
{
	// neutral signal 
	if( factor_type == 'N' || factor_type == 'n' )
	{ return; }

	// promoter signal 
	if( factor_type == 'P' || factor_type == 'p' )
	{
		// up = sum of all (scaled) promoter signals 
		up += factor; 
		// cap up signal at 1.0 
		if( up > 1.0 )
		{ up = 1.0; }
		no_promoters = false; 
		return; 
	}

	// inhibitor signal 
	if( factor_type == 'I' || factor_type == 'i' )
	{
		if( factor > down )
		{ down = factor; }
		if( down > 1 )
		{ down = 1.0; }
		no_inhibitors = false; 
		return; 
	}

	return; 
}

void up_down_signal::add_effect( double factor, std::string factor_type )
{
	this->add_effect( factor, factor_type[0] ); 
	return; 
}

double up_down_signal::compute_effect( void )
{
	double UP = up;
	if( no_promoters )
	{ UP = 1.0; }
	return UP * (1.0 - down ); 
}

void up_down_signal::reset( void )
{ 
	up = 0.0;
	down = 0.0; 
	no_promoters = true; 
	no_inhibitors = true; 
	return; 
}

void up_down_signal::display( void )
{
	std::cout << "up    : " << up <<   " (no promoters : " << (int) no_promoters << ")" << std::endl;
	std::cout << "down  : " << down << " (no inhibiters: " << (int) no_inhibitors << ")" << std::endl;

	std::cout << "effect: " << compute_effect() << std::endl; 

	std::cout << std::endl; 
	return; 
}



void A_phenotype( Cell* pCell, Phenotype& phenotype, double dt )
{
	// housekeeping 
	static Cell_Definition* pCD  = find_cell_definition("A");
	static int nApoptosis = pCD->phenotype.death.find_death_model_index( "Apoptosis"); 
	static int nNecrosis  = pCD->phenotype.death.find_death_model_index( "Necrosis"); 

/*
	if( phenotype.death.dead == true )
	{
		#pragma omp critical 
		{
			std::cout << (long int) pCell << std::endl ; 
			std::cout << phenotype.cycle.model().name << std::endl; 
			std::cout << "apop : " << phenotype.death.rates[nApoptosis] << std::endl; 
			std::cout << "necro: " << phenotype.death.rates[nNecrosis] << std::endl << std::endl; 
			system("sleep 1");
		} 
	}
*/

	if( phenotype.death.dead == true )
	{

		phenotype.secretion.set_all_secretion_to_zero(); 
		phenotype.secretion.set_all_uptake_to_zero(); 
		phenotype.motility.is_motile = false; 

		pCell->functions.update_phenotype = NULL; 
		return; 
	}
	
	// sample A, B, C, resource, and pressure 
	static int nA = microenvironment.find_density_index( "signal A" ); 
	static int nB = microenvironment.find_density_index( "signal B" ); 
	static int nC = microenvironment.find_density_index( "signal C" ); 
	static int nR = microenvironment.find_density_index( "resource" ); 

	double A = pCell->nearest_density_vector()[nA];
	double B = pCell->nearest_density_vector()[nB];
	double C = pCell->nearest_density_vector()[nC];
	double R = pCell->nearest_density_vector()[nR];
	double p = pCell->state.simple_pressure; 

	// necrotic death rate 
	static double base_necrosis_rate = pCD->phenotype.death.rates[nNecrosis];
	static double necrosis_threshold = parameters.doubles("A_necrosis_threshold");
	phenotype.death.rates[nNecrosis] = 0.0;
	if( R < necrosis_threshold )
	{
		phenotype.death.rates[nNecrosis] = base_necrosis_rate; 
		phenotype.death.rates[nNecrosis] *= (1.0 - R / necrosis_threshold);
	}

	// cycle rate 
	static double base_cycle_rate = pCD->phenotype.cycle.data.transition_rate(0,0); 
	phenotype.cycle.data.transition_rate(0,1) = base_cycle_rate;
	phenotype.cycle.data.transition_rate(0,1) *= R; 
	if( p > parameters.doubles( "A_cycle_pressure_threshold") )
	{ phenotype.cycle.data.transition_rate(0,1) = 0.0; }

	up_down_signal sig; 

	// A 
	sig.add_effect( A , parameters.strings("A_cycle_A") );
	// B
	sig.add_effect( B , parameters.strings("A_cycle_B") );
	// C 
	sig.add_effect( C , parameters.strings("A_cycle_C") );
	
	phenotype.cycle.data.transition_rate(0,1) *= sig.compute_effect(); 

	// apoptotic rate 

	double base_apoptosis_rate = pCD->phenotype.death.rates[nApoptosis]; 
	phenotype.death.rates[nApoptosis] = base_apoptosis_rate; 

	sig.reset();

	// A 
	sig.add_effect( A , parameters.strings("A_death_A") );
	// B
	sig.add_effect( B , parameters.strings("A_death_B") );
	// C 
	sig.add_effect( C , parameters.strings("A_death_C") );	
	// R 
	sig.add_effect( C , parameters.strings("A_death_R") );	
	
	phenotype.death.rates[nApoptosis] *= sig.compute_effect(); 
	if( p > parameters.doubles("A_apoptosis_pressure_threshold") )
	{
		phenotype.death.rates[nApoptosis] = 10; 
	}

	// speed 
	phenotype.motility.migration_speed = pCD->phenotype.motility.migration_speed; 

	sig.reset(); 

	// A 
	sig.add_effect( A , parameters.strings("A_speed_A") );
	// B
	sig.add_effect( B , parameters.strings("A_speed_B") );
	// C 
	sig.add_effect( C , parameters.strings("A_speed_C") );	
	// R 
	sig.add_effect( C , parameters.strings("A_speed_R") );	

	phenotype.motility.migration_speed *= sig.compute_effect();

	// secretion 

	phenotype.secretion.secretion_rates[nA] = 
		pCD->phenotype.secretion.secretion_rates[nA]; 

	sig.reset(); 
	// A 
	sig.add_effect( A , parameters.strings("A_signal_A") );
	// B
	sig.add_effect( B , parameters.strings("A_signal_B") );
	// C 
	sig.add_effect( C , parameters.strings("A_signal_C") );	
	// R 
	sig.add_effect( C , parameters.strings("A_signal_R") );	

	phenotype.secretion.secretion_rates[nA] *= sig.compute_effect();

	return; 
}






void B_phenotype( Cell* pCell, Phenotype& phenotype, double dt )
{
	if( phenotype.death.dead == true )
	{
		phenotype.secretion.set_all_secretion_to_zero(); 
		phenotype.secretion.set_all_uptake_to_zero(); 
		phenotype.motility.is_motile = false; 

		pCell->functions.update_phenotype = NULL; 
		return; 
	}
	
	// housekeeping 
	static int nApoptosis = cell_defaults.phenotype.death.find_death_model_index( "apoptosis"); 
	static int nNecrosis  = cell_defaults.phenotype.death.find_death_model_index( "necrosis"); 
	static Cell_Definition* pCD  = find_cell_definition("B");

	// sample A, B, C, resource;
	static int nA = microenvironment.find_density_index( "signal A" ); 
	static int nB = microenvironment.find_density_index( "signal B" ); 
	static int nC = microenvironment.find_density_index( "signal C" ); 
	static int nR = microenvironment.find_density_index( "resource" ); 

	double A = pCell->nearest_density_vector()[nA];
	double B = pCell->nearest_density_vector()[nB];
	double C = pCell->nearest_density_vector()[nC];
	double R = pCell->nearest_density_vector()[nR];
	double p = pCell->state.simple_pressure; 

	// necrotic death rate 
	static double base_necrosis_rate = pCD->phenotype.death.rates[nNecrosis];
	phenotype.death.rates[nNecrosis] *= (1.0-R);

	// cycle rate 
	static double base_cycle_rate = pCD->phenotype.cycle.data.transition_rate(0,0); 
	phenotype.cycle.data.transition_rate(0,1) = base_cycle_rate;
	phenotype.cycle.data.transition_rate(0,1) *= R; 
	if( p > parameters.doubles( "B_cycle_pressure_threshold") )
	{ phenotype.cycle.data.transition_rate(0,1) = 0.0; }

	double factor = 1.0; 
	char temp = parameters.strings("B_cycle_A" )[0]; 
	if( temp == 'p' || temp == 'P' ) // promotes cycling 
	{ factor *= A; }
	if( temp == 'i' || temp == 'I' ) // inhibits cycling 
	{ factor *= (1-A); }
	phenotype.cycle.data.transition_rate(0,1) *= factor;

	factor = 1.0; 
	temp = parameters.strings("B_cycle_B" )[0]; 
	if( temp == 'p' || temp == 'P' ) // promotes cycling 
	{ factor *= B; }
	if( temp == 'i' || temp == 'I' ) // inhibits cycling 
	{ factor *= (1-B); }
	phenotype.cycle.data.transition_rate(0,1) *= factor;

	factor = 1.0; 
	temp = parameters.strings("B_cycle_C" )[0]; 
	if( temp == 'p' || temp == 'P' ) // promotes cycling 
	{ factor *= C; }
	if( temp == 'i' || temp == 'I' ) // inhibits cycling 
	{ factor *= (1-C); }
	phenotype.cycle.data.transition_rate(0,1) *= factor;

	// apoptotic rate 

	double base_apoptosis_rate = pCD->phenotype.death.rates[nApoptosis]; 
	phenotype.death.rates[nApoptosis] = base_apoptosis_rate; 

	factor = 1.0; 
	temp = parameters.strings("B_death_A" )[0]; 
	if( temp == 'p' || temp == 'P' ) // promotes death 
	{ factor *= A; }
	if( temp == 'i' || temp == 'I' ) // inhibits death 
	{ factor *= (1-A); }
	phenotype.death.rates[nApoptosis] *= factor;

	factor = 1.0; 
	temp = parameters.strings("B_death_B" )[0]; 
	if( temp == 'p' || temp == 'P' ) // promotes death 
	{ factor *= B; }
	if( temp == 'i' || temp == 'I' ) // inhibits death 
	{ factor *= (1-B); }
	phenotype.death.rates[nApoptosis] *= factor;

	factor = 1.0; 
	temp = parameters.strings("B_death_C" )[0]; 
	if( temp == 'p' || temp == 'P' ) // promotes death 
	{ factor *= C; }
	if( temp == 'i' || temp == 'I' ) // inhibits death 
	{ factor *= (1-C); }
	phenotype.death.rates[nApoptosis] *= factor;

	factor = 1.0; 
	temp = parameters.strings("B_death_R" )[0]; 
	if( temp == 'p' || temp == 'P' ) // promotes death 
	{ factor *= R; }
	if( temp == 'i' || temp == 'I' ) // inhibits death 
	{ factor *= (1-R); }
	phenotype.death.rates[nApoptosis] *= factor;

	// speed 
	phenotype.motility.migration_speed = pCD->phenotype.motility.migration_speed; 

	factor = 1.0; 
	temp = parameters.strings("B_speed_A" )[0]; 
	if( temp == 'p' || temp == 'P' ) // promotes motility 
	{ factor *= A; }
	if( temp == 'i' || temp == 'I' ) // inhibits motility 
	{ factor *= (1-A); }
	phenotype.motility.migration_speed *= factor;

	factor = 1.0; 
	temp = parameters.strings("B_speed_B" )[0]; 
	if( temp == 'p' || temp == 'P' ) // promotes motility 
	{ factor *= B; }
	if( temp == 'i' || temp == 'I' ) // inhibits motility 
	{ factor *= (1-B); }
	phenotype.motility.migration_speed *= factor;

	factor = 1.0; 
	temp = parameters.strings("B_speed_C" )[0]; 
	if( temp == 'p' || temp == 'P' ) // promotes motility 
	{ factor *= C; }
	if( temp == 'i' || temp == 'I' ) // inhibits motility 
	{ factor *= (1-C); }
	phenotype.motility.migration_speed *= factor;

	factor = 1.0; 
	temp = parameters.strings("B_speed_R" )[0]; 
	if( temp == 'p' || temp == 'P' ) // promotes motility 
	{ factor *= R; }
	if( temp == 'i' || temp == 'I' ) // inhibits motility 
	{ factor *= (1-R); }
	phenotype.motility.migration_speed *= factor;

	return; 
}

void C_phenotype( Cell* pCell, Phenotype& phenotype, double dt )
{
	if( phenotype.death.dead == true )
	{
		phenotype.secretion.set_all_secretion_to_zero(); 
		phenotype.secretion.set_all_uptake_to_zero(); 
		phenotype.motility.is_motile = false; 

		pCell->functions.update_phenotype = NULL; 
		return; 
	}
	
	// housekeeping 
	static int nApoptosis = cell_defaults.phenotype.death.find_death_model_index( "apoptosis"); 
	static int nNecrosis  = cell_defaults.phenotype.death.find_death_model_index( "necrosis"); 
	static Cell_Definition* pCD  = find_cell_definition("C");

	// sample A, B, C, resource;
	static int nA = microenvironment.find_density_index( "signal A" ); 
	static int nB = microenvironment.find_density_index( "signal B" ); 
	static int nC = microenvironment.find_density_index( "signal C" ); 
	static int nR = microenvironment.find_density_index( "resource" ); 

	double A = pCell->nearest_density_vector()[nA];
	double B = pCell->nearest_density_vector()[nB];
	double C = pCell->nearest_density_vector()[nC];
	double R = pCell->nearest_density_vector()[nR];
	double p = pCell->state.simple_pressure; 

	// necrotic death rate 
	static double base_necrosis_rate = pCD->phenotype.death.rates[nNecrosis];
	phenotype.death.rates[nNecrosis] *= (1.0-R);

	// cycle rate 
	static double base_cycle_rate = pCD->phenotype.cycle.data.transition_rate(0,0); 
	phenotype.cycle.data.transition_rate(0,1) = base_cycle_rate;
	phenotype.cycle.data.transition_rate(0,1) *= R; 
	if( p > parameters.doubles( "B_cycle_pressure_threshold") )
	{ phenotype.cycle.data.transition_rate(0,1) = 0.0; }

	double factor = 1.0; 
	char temp = parameters.strings("C_cycle_A" )[0]; 
	if( temp == 'p' || temp == 'P' ) // promotes cycling 
	{ factor *= A; }
	if( temp == 'i' || temp == 'I' ) // inhibits cycling 
	{ factor *= (1-A); }
	phenotype.cycle.data.transition_rate(0,1) *= factor;

	factor = 1.0; 
	temp = parameters.strings("C_cycle_B" )[0]; 
	if( temp == 'p' || temp == 'P' ) // promotes cycling 
	{ factor *= B; }
	if( temp == 'i' || temp == 'I' ) // inhibits cycling 
	{ factor *= (1-B); }
	phenotype.cycle.data.transition_rate(0,1) *= factor;

	factor = 1.0; 
	temp = parameters.strings("C_cycle_C" )[0]; 
	if( temp == 'p' || temp == 'P' ) // promotes cycling 
	{ factor *= C; }
	if( temp == 'i' || temp == 'I' ) // inhibits cycling 
	{ factor *= (1-C); }
	phenotype.cycle.data.transition_rate(0,1) *= factor;

	// apoptotic rate 

	double base_apoptosis_rate = pCD->phenotype.death.rates[nApoptosis]; 
	phenotype.death.rates[nApoptosis] = base_apoptosis_rate; 

	factor = 1.0; 
	temp = parameters.strings("C_death_A" )[0]; 
	if( temp == 'p' || temp == 'P' ) // promotes death 
	{ factor *= A; }
	if( temp == 'i' || temp == 'I' ) // inhibits death 
	{ factor *= (1-A); }
	phenotype.death.rates[nApoptosis] *= factor;

	factor = 1.0; 
	temp = parameters.strings("C_death_B" )[0]; 
	if( temp == 'p' || temp == 'P' ) // promotes death 
	{ factor *= B; }
	if( temp == 'i' || temp == 'I' ) // inhibits death 
	{ factor *= (1-B); }
	phenotype.death.rates[nApoptosis] *= factor;

	factor = 1.0; 
	temp = parameters.strings("C_death_C" )[0]; 
	if( temp == 'p' || temp == 'P' ) // promotes death 
	{ factor *= C; }
	if( temp == 'i' || temp == 'I' ) // inhibits death 
	{ factor *= (1-C); }
	phenotype.death.rates[nApoptosis] *= factor;

	factor = 1.0; 
	temp = parameters.strings("C_death_R" )[0]; 
	if( temp == 'p' || temp == 'P' ) // promotes death 
	{ factor *= R; }
	if( temp == 'i' || temp == 'I' ) // inhibits death 
	{ factor *= (1-R); }
	phenotype.death.rates[nApoptosis] *= factor;

	// speed 
	phenotype.motility.migration_speed = pCD->phenotype.motility.migration_speed; 

	factor = 1.0; 
	temp = parameters.strings("C_speed_A" )[0]; 
	if( temp == 'p' || temp == 'P' ) // promotes motility 
	{ factor *= A; }
	if( temp == 'i' || temp == 'I' ) // inhibits motility 
	{ factor *= (1-A); }
	phenotype.motility.migration_speed *= factor;

	factor = 1.0; 
	temp = parameters.strings("C_speed_B" )[0]; 
	if( temp == 'p' || temp == 'P' ) // promotes motility 
	{ factor *= B; }
	if( temp == 'i' || temp == 'I' ) // inhibits motility 
	{ factor *= (1-B); }
	phenotype.motility.migration_speed *= factor;

	factor = 1.0; 
	temp = parameters.strings("C_speed_C" )[0]; 
	if( temp == 'p' || temp == 'P' ) // promotes motility 
	{ factor *= C; }
	if( temp == 'i' || temp == 'I' ) // inhibits motility 
	{ factor *= (1-C); }
	phenotype.motility.migration_speed *= factor;

	factor = 1.0; 
	temp = parameters.strings("C_speed_R" )[0]; 
	if( temp == 'p' || temp == 'P' ) // promotes motility 
	{ factor *= R; }
	if( temp == 'i' || temp == 'I' ) // inhibits motility 
	{ factor *= (1-R); }
	phenotype.motility.migration_speed *= factor;

	return; 
}

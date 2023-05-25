# -*- coding: utf-8 -*-
"""
Created on Mon Jan 16 16:20:38 2023

@author: Kdeiri
"""

import utilities as util 
import pathlib
import simscale_sdk as sim_sdk
import pandas as pd 
##################################

cht = util.ConjugateHeatTransfer()

"""Setup the API connection"""
cht.set_api_connection()

"""Create Project"""
cht.create_project( name = "CoolingPlate_CHTv2-API",
                   description = "Running the cooling plate example using an automated workflow via the API")

"""Upload Geometry"""
#Provide the name of the files to upload, if it is a directory simply give the name,
#if it is a file then add the file extension to the name ex: example.stl
name_of_files_to_upload = ["Cooling_plate_V1.x_t"]
geometry_path = []
# for i, cad in enumerate(name_of_files_to_upload):  # use in case you're looping over mutliple geo
geometry_path.append(pathlib.Path().cwd() / "Geometries" / name_of_files_to_upload[0])
cht.upload_geometry(name_of_files_to_upload[0], geometry_path[0])

"""Define geometry mappings""" 
#(cooling plate example)
volume_cad  = {"flow_volume" : "NS_Fluid", "tool_body" : "NS_Al"}
surface_cad = {"inlet": "NS_Inlet", "outlet": "NS_Outlet", "battery_load": "NS_Battery"}

#Get the entity ID associated with each CAD part 
for key, value in volume_cad.items(): 
    cht.get_single_entity_name(cht.project_id, cht.geometry_id, key = key,
                               attributes=["SDL/TYSA_NAME"], values=[value])
    
#Get the entity ID associated with each predefined surface
for key, value in surface_cad.items():
    cht.get_single_entity_name(cht.project_id, cht.geometry_id, key = key ,
                               _class = 'face' ,attributes=["SDL/TYSA_NAME"],
                               values=[value])
    
for key, value in cht.single_entity.items(): 
    print("{k} : {v}".format(k = key, v = value))

"""Simulation Setup"""
#Global Settings
cht.set_compressible(state = False)
cht.set_turbulence_model(model = "KOMEGASST") #KOMEGASST ; NONE (laminar)
cht.set_gravity_direction(g_x = 0, g_y = 0, g_z = -9.81) 

#-----------------------------
#Define Boundary Conditions
cht.flow_rate_velocity_inlet(flowrate = 0.13 , temp = 18, type = "mass",
                             name = 'Velocity Inlet', key = 'inlet')

cht.pressure_outlet_bc(value = 0 ,name = "Pressure Outlet", unit = 'bar', key = "outlet" )

cht.external_wall_heat_flux_bc(method = 'FIXED_POWER', power = 2000, key_list = ['battery_load'])

cht.set_boundary_conditions()

#-----------------------------
#Define Simulation Control Settings
cht.set_simulation_end_time(time = 1000)
cht.set_simulation_time_step(time_step = 1)
cht.set_simulation_write_controls(write_interval = 1000)
cht.set_simulation_max_run_time(max_run_time = 40000)
cht.set_simulation_control()

#-----------------------------
#Define Result Controls 
cht.set_area_averages(name = 'inlet-outlet', write_interval = 10,
                      key_list = ["inlet", 'outlet'])

cht.set_area_volumes(name = 'inlet-outlet_vol', write_interval = 10,
                     key_list = ["inlet", 'outlet'])

cht.set_result_control_items()

#-----------------------------
#create simulation
cht.set_simulation_spec( simulation_name = "chtv2_api")
#-----------------------------

#Update the spec with the materials 
cht.add_fluid_material(fluid_name = 'Water'  , key_list = ['flow_volume'])
cht.add_solid_material(solid_name = 'Aluminium', key_list = ['tool_body'])

#-----------------------------    
#Mesh settings
cht.set_mesh_layer_settings(num_of_layers = 3, total_rel_thickness = 0.4, growth_rate = 1.5)
cht.set_advanced_mesh_settings(small_feature_tolerance = 5E-5, gap_ref_factor = 0.05, gradation_rate = 1.22)
cht.complete_mesh_settings(mesh_name = "Mesh_test", fineness = 0.1, physics_based_meshing = True)
cht.estimate_mesh_operation()
cht.start_meshing_operation(run_state = False)

#-----------------------------    
#Sanity checks
cht.check_simulation_and_mesh_settings()

#-----------------------------
#Start Simulation
cht.estimate_simulation(maximum_cpu_consumption_limit = 200)
cht.create_simulation(sim_name = "Run 1")
cht.start_simulation_run(wait_for_results = False)


#############################################################
"""example: boundary conditions"""
# cht.pressure_inlet_bc(value = 1, temp = 95 ,name = "Pressure Inlet", unit = 'bar', key = "inlet" )

# cht.pressure_outlet_bc(value = 0 ,name = "Pressure Outlet", unit = 'bar', key = "outlet" )

# cht.constant_velocity_inlet_bc(speed_x = 1, speed_y = 0, speed_z = 2, temp = 10, key = 'inlet')

# cht.flow_rate_velocity_inlet(flowrate = 0.1, temp = 10, type = "mass", key = 'inlet') #volumetric, mass


# cht.external_wall_heat_flux_bc( amb_temp = 20, htc = 10 ,
#                                method = "DERIVED"  ,  # DERIVED ; FIXED ; FIXED_POWER
#                                name = 'ExternalWalls', 
#                                key_list = ["tool_body_top", "tool_body_bottom",
#                                            "tool_body_side1", "tool_body_side2",
#                                            "tool_body_side3", "tool_body_side4"])

# cht.natural_convection_inlet_outlet_bc(pressure = 0, temp = 24 , key_list = ['inlet', 'outlet'])

"example: geometry primitives"
#single point
# cht.set_single_geometry_primitive_point(name = "test_point", pos_x = 0.0420, pos_y= 0.0438, pos_z = 0)   

#list of points
# probe_points_path = pathlib.Path().cwd() / "probe_points" / "probe_list.txt"                              
# cht.set_multiple_geometry_primitive_points(path_to_csv = probe_points_path)
 
#cylinder 
# cylinder1_uuid = cht.set_cylinderical_geometry_primitive(r = 0.02,      
#                                                          len_x= 0, len_y= 0, len_z= 0.1,
#                                                          vec_x= 0, vec_y= 0, vec_z = 1,
#                                                          name = 'cylinder')

#cartesian box
# box_uuid = cht.set_cartesian_box_geometry_primitive(min_x = 0, min_y = 0, min_z = 0,
#                                                     max_x = 0.1, max_y = 0.1, max_z = 0.1,
#                                                     name= "box")

"""example: advanced concepts"""
#power source using a part 
# cht.set_power_sources(power = 150, name = "power source abs", method = "ABSOLUTE",
#                       key_list = ["tool_body"])

#momentum source using a geometry primitive
# cht.set_momentum_sources(speed_x = 1, speed_y = 0, speed_z = 0, name = "Momentum Source 1",
#                          key_list = [], geo_primitive_ids = cylinder1_uuid )

# cht.set_advanced_concepts()

"""example: download results"""
# cht.find_simulation("chtv2_api")
# cht.find_run("Design4_test")
# cht.get_simulation_results()
# cht.get_probe_point_results(name = "test_probe", field = 'T')
# cht.get_probe_point_results(name = "multi_test_probe", field = 'T')
# cht.get_surface_data_results(data_type = 'average', name = 'inlet-outlet', field = "T")
# cht.get_surface_data_results(data_type = 'average', name = 'inlet-outlet', field = "p")
# cht.get_surface_data_results(data_type = 'integral', name = 'inlet-outlet_vol', field = "Uy")
# cht.get_simulation_case_files()

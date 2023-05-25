# -*- coding: utf-8 -*-
"""
Created on Mon Jan 16 16:21:14 2023

@author: Kdeiri
"""

import os
import time
import zipfile
import pathlib
import isodate
import urllib3
import shutil
import pandas as pd

from simscale_sdk import Configuration, ApiClient, ProjectsApi, StorageApi, GeometryImportsApi, GeometriesApi, \
    MeshOperationsApi, SimulationsApi, SimulationRunsApi, ReportsApi, Project, GeometryImportRequest, ApiException, \
    MaterialsApi, MaterialGroupType, MaterialUpdateRequest, MaterialUpdateOperation, MaterialUpdateOperationReference
from simscale_sdk import ConvectiveHeatTransfer, CoupledConjugateHeatTransfer, FluidModel, DimensionalVectorAcceleration, FluidInitialConditions, \
    AdvancedConcepts, ConvectiveHeatTransferMaterials, CoupledConjugateHeatTransferMaterials, TopologicalReference, \
    FluidNumerics, RelaxationFactor, DimensionalPressure, ResidualControls, Tolerance, \
    FluidSolvers, Schemes, TimeDifferentiationSchemes, GradientSchemes, DivergenceSchemes, LaplacianSchemes, \
    InterpolationSchemes, SurfaceNormalGradientSchemes, VelocityInletBC, FixedValueVBC, DimensionalVectorFunctionSpeed, \
    ComponentVectorFunction, ConstantFunction, FixedValueTBC, DimensionalFunctionTemperature, PressureOutletBC, \
    FixedValuePBC, DimensionalFunctionPressure, WallBC, NoSlipVBC, FluidSimulationControl, DimensionalTime, \
    TimeStepWriteControl, ScotchDecomposeAlgorithm, FluidResultControls, AreaAverageResultControl, \
    ProbePointsResultControl, AbsolutePowerSource, DimensionalFunctionPower, IncompressibleMaterial, NewtonianViscosityModel, DimensionalKinematicViscosity, \
    DimensionalDensity, DimensionalThermalExpansionRate, DimensionalTemperature, DimensionalSpecificHeat, PBICGSolver, DILUPreconditioner, ILUCpPreconditioner, \
    SolidCompressibleMaterial, ConstIsoTransport, IsotropicConductivity, DimensionalFunctionThermalConductivity, \
    DimensionalSpecificHeat, HConstThermo, RhoConstEquationOfState, AutomaticMeshSizingSimmetrix, CrossPlaneOrthotropicConductivity, ConstCrossPlaneOrthotropicTransport

    
from simscale_sdk import GeometryImportRequestLocation, GeometryImportRequestOptions, Point, DimensionalVectorLength, \
    DecimalVector
from simscale_sdk import SimulationSpec, MeshOperation, SimmetrixMeshingFluid, AutomaticLayerOn, SimulationRun
from simscale_sdk import UserInputCameraSettings, ProjectionType, Vector3D, ModelSettings, Part, ScalarField, \
    ScreenshotOutputSettings, Color, ResolutionInfo, ScreenshotReportProperties, ReportRequest
    
import simscale_sdk as sim_sdk


class ConjugateHeatTransfer(): 
    
    def __init__(self): 
        
        #API Variables
        self.api_key        = ""
        self.api_url        = ""
        self.api_key_header = "X-API-KEY"
        self.version        = "/v0"
        self.host           = ""
        self.server         = "prod"
        
        #Client Variables
        self.api_client  = None
        self.project_api = None
        self.storage_api = None
        self.geometry_import_api = None
        self.geometry_api = None
        self.materials_api = None 
        self.mesh_operation_api = None 
        self.simulation_api = None
        self.simulation_run_api = None
        self.table_import_api =None
        self.reports_api = None
        
        #Project Variables 
        self.project_name = ""
        self.project_id   = ""
        
        #Geometry Variables
        self.geometry_name = ""
        self.geometry_id   = ""
        self.geometry_path = ""
        
        #Geometry Mapping 
        self.single_entity     = {} #for later: separate the faces from volumes and add a validation rule against duplicate assignments for the same entity (materials)
        self.multiple_entities = {}
        
        #Global Simulation Settings
        self.compressible = False
        self.turbulence_model = None
        self.fluid_model_gravity = None 
        self.initial_conditions = None 
        self.radiation = False 
        
        #Material Variables
        self.fluid_material = []
        self.solid_material = []
        
        #Numerics 
        self.fluid_numerics = None
        
        #Boundary Condition Variables 
        self.boundary_conditions = []
        self.velocity_inlet      = []
        self.velocity_outlet     = []
        self.pressure_inlet      = []
        self.pressure_outlet     = []
        self.wall_bc             = []
        self.natural_convection  = []
        
        #Advanced Concept Variables
        self.advanced_concepts = sim_sdk.AdvancedConcepts()
        self.power_sources    = []
        self.porous_media     = []
        self.momentum_sources  = []
        self.thermal_resistance_networks = []
        
        #Simulation Control Variables
        self.simulation_control = None
        self.end_time = None
        self.delta_t  = None
        self.write_control = None 
        self.max_run_time = None 
        
        #Result Control Variables
        self.result_control = sim_sdk.FluidResultControls()
        self.surface_data = []
        self.volume_data  = []
        self.probe_points  = []
        self.field_calculations = []
        
        #Geometry Primitive Variables 
        self.geometry_primitive_uuid = None
        self.geometry_primitive_uuid_list = []
        # self.geometry_primitive_body_uuid = []
        
        #Contact Definition
        self.contact_detection = "AUTO"
        
        #Mesh Variables
        self.mesh_operation = None 
        self.mesh_operation_id = None
        self.automatic_layer_settings = None
        self.advanced_mesh_settings = None
        self.mesh_refinement = [] 
        self.mesh_id = None
        self.mesh_max_runtime = None

        #Simulation Creation Variables
        self.model = None 
        self.sim_max_run_time = None
        self.simulation_spec = None
        self.simulation_id   = None
        self.simulation_run  = None 
        self.run_id = None
        
        #Simulation Results Download Variables 
        self.simulation_results = None 
        self.probe_point_plot_info = None 
        self.probe_point_plot_data_response = None
        self.solution_info = None
        
        
    """"Functions to setup the API connection """
    
    def _get_variables_from_env(self):
        
        '''
        looks in your environment and reads the API variables
        
        SimScale API key and URL are read if they are set in the 
        environment as:
            
            SIMSCALE_API_KEY
            SIMSCALE_API_URL

        Returns
        -------
        None.

        '''
        try:
            self.api_key = os.getenv('SIMSCALE_API_KEY')
            self.api_url = os.getenv('SIMSCALE_API_URL')
            self.host = self.api_url + self.version
        except:
            raise Exception("Cannot get Keys from Environment Variables")
        
    def check_api(self):
        '''
        Check API key is set, returns boolean True if not set.
    
        Raises
        ------
        Exception
            If the API key and URL is not set, rasie an exception.
    
        Returns
        -------
        is_not_existent : boolean
            True if not set.
    
        '''
        is_not_existent = not os.getenv("SIMSCALE_API_KEY") or not os.getenv("SIMSCALE_API_URL")
        if is_not_existent:
            raise Exception(
                "Either `SIMSCALE_API_KEY` or `SIMSCALE_API_URL`",
                " environment variable is missing.")
            return is_not_existent
        else:
            print("SimScale API Key and URL found in environment variables.")

    def set_api_connection(self, version=0, server='prod'):
        '''
        Reads API key and URL and returns API clients required.
        
        ----------
        version : int
            Version of SimScale API, at time of writing, only 0 is valid.
            Default is 0.
        
        Returns
        -------
        api_client : object
            An API client that represents the user, and their login 
            credentials.
        api_key_header : object
        api_key : string
            A string that is your API key, read from the environment 
            variables.
        credential : SimscaleCredentials object
            An object contain api keys and credential information
    
        '''
        #Get the API url and key variables from env variables and do a sanity check 
        self._get_variables_from_env()
        self.check_api()
        
        #Setup the API configuration (define host and link the associated key)
        configuration = sim_sdk.Configuration()
        configuration.host = self.host
        configuration.api_key = {self.api_key_header: self.api_key}
        
        #Setup the API client connection 
        self.api_client = sim_sdk.ApiClient(configuration)
        retry_policy = urllib3.Retry(connect=5, read=5, redirect=0, status=5, backoff_factor=0.2)
        self.api_client.rest_client.pool_manager.connection_pool_kw["retries"] = retry_policy
       
        #Define the required API clients for the simulation 
        self.project_api = sim_sdk.ProjectsApi(self.api_client)
        self.storage_api = sim_sdk.StorageApi(self.api_client)
        self.geometry_import_api = sim_sdk.GeometryImportsApi(self.api_client)
        self.geometry_api = sim_sdk.GeometriesApi(self.api_client)
        self.mesh_operation_api = sim_sdk.MeshOperationsApi(self.api_client)
        self.materials_api = sim_sdk.MaterialsApi(self.api_client)
        self.simulation_api = sim_sdk.SimulationsApi(self.api_client)
        self.simulation_run_api = sim_sdk.SimulationRunsApi(self.api_client)
        self.table_import_api = sim_sdk.TableImportsApi(self.api_client)
        self.reports_api = sim_sdk.ReportsApi(self.api_client) 
        self.wind_api = sim_sdk.WindApi(self.api_client)
        
    def create_project(self, name, description, measurement_system = "SI"):
        '''
        Take a name and description and create a new workbench project

        Parameters
        ----------
        name : str
            A string with the exact name for the new project.
            
        description : str
            A string with the exact description for the new project.

        Returns
        -------
        None.

        '''
        
        try:
            #Check if the project already exists
            projects = self.project_api.get_projects(limit=1000).to_dict()['embedded']
            found = None
            for project in projects:
                if project['name'] == name:
                    found = project
                    print('Project found: \n' + str(found['name']))
                    break
            if found is None:
                raise Exception('could not find project with name: ' + name)
            
            self.project_id = found['project_id']
            self.project_name = name
            print("Cannot create project with the same name, using existing project")
        except:
            #If not then create a new project
            project = sim_sdk.Project(name=name, description=description,
                                      measurement_system = measurement_system)
            project = self.project_api.create_project(project)
            self.project_id = project.project_id
            self.project_name = name        
            
    def zip_cad_for_upload(self, file_name, base_path): 
     
         '''
         Take a list of the CAD file names and their associated path then zips 
         each CAD file separately and submits it for upload 
         
         Note: Have all the CAD files stored in the same directory
         
         Parameters
         ----------
         file_name : list
             A list with the exact names of the CAD files to upload
             
         base_path : Pathlib Path
             path to the directory that contains the CAD files 
     
         Returns
         -------
         geometry_path : path of the zipped file
     
         '''
         geometry_path = []
         
         #Loop over the CAD files needed for upload
         for cad in file_name:
             
             #Get the path of each CAD file
             path = base_path / cad
             
             # The output_filename variable saves the zip file at a desired path; 
             # in this case it is same directory
             output_filename = path
              
             #Retruns a zip file(s) path of the associated CAD, 
             geometry_path.append(shutil.make_archive(output_filename, 'zip', path)) 
     
         return geometry_path
 
    def upload_geometry(self, name, path=None, units="m", _format="PARASOLID"):
        
        '''
        Upload a geometry to the SimScale platform to a preassigned project.
        
        Parameters
        ----------
        name : str
            The name given to the geometry.
            
        path : pathlib.Path, optional
            The path to a geometry to upload. 
            
        units : str, optional
            the unit in which to upload the geometry to SimScale.
            
            The default is "m".
            
        _format : str, optional
            The file format. 
            
            The default is "STL".
            
        facet_split : bool, optional
            Decide on weather to split facet geometry (such as .stl file 
            types). We prefer not to do this for API use.
            
            The default is False.
    
        Raises
        ------
        TimeoutError
            DESCRIPTION.
    
        Returns
        -------
        None.
    
        '''
        self.geometry_name = name
        
        #Check if the geometry already exists
        try:
            project_id = self.project_id
            geometry_api = self.geometry_api
        
            geometries = geometry_api.get_geometries(project_id).to_dict()['embedded']
            found = None
            for geometry in geometries:
                if geometry['name'] == name:
                    found = geometry
                    print('Geometry found: \n' + str(found['name']))
                    break
                        
            if found is None:
                raise Exception('could not find geometry with id: ' + name)
                
            self.geometry_name = found
            self.geometry_id = found["geometry_id"]
            print("Cannot upload geometry with the same name, using existing geometry")
    
        except:
            
            self.geometry_path = path
    
            storage = self.storage_api.create_storage()
            with open(self.geometry_path, 'rb') as file:
                self.api_client.rest_client.PUT(url=storage.url, headers={'Content-Type': 'application/octet-stream'},
                                                body=file.read())
            self.storage_id = storage.storage_id
    
            geometry_import = sim_sdk.GeometryImportRequest(
                name=name,
                location=sim_sdk.GeometryImportRequestLocation(self.storage_id),
                format=_format,
                input_unit=units,
                options=sim_sdk.GeometryImportRequestOptions(facet_split= False, sewing=False, improve=True,
                                                         optimize_for_lbm_solver=False),
            )
    
            geometry_import = self.geometry_import_api.import_geometry(self.project_id, geometry_import)
            geometry_import_id = geometry_import.geometry_import_id
    
            geometry_import_start = time.time()
            while geometry_import.status not in ('FINISHED', 'CANCELED', 'FAILED'):
                # adjust timeout for larger geometries
                if time.time() > geometry_import_start + 900:
                    raise TimeoutError()
                time.sleep(10)
                geometry_import = self.geometry_import_api.get_geometry_import(self.project_id, geometry_import_id)
                print(f'Geometry import status: {geometry_import.status}')
            self.geometry_id = geometry_import.geometry_id

    """Functions to get geometry mappings for topological entity assignment"""
    
    def get_single_entity_name(self, project_id, geometry_id, key ,**kwargs):
        
        '''
        Retrieve the entity ID of a face/part 
        
        Parameters
        ----------
        project_id : str
            The id of the project being worked on.
            
        geometry_id : str
            The geometry that some entities are being retrieved from.  
            
        key : str
            the desired faces/parts for later topological assignment are stored 
            in a dict. each face/part has a specific user defined key, that is 
            passed to this function. In the end, this key would point to the 
            desired entitiy id of a face/part
            
        Returns
        -------
        None.
    
        '''        
        entity = self.geometry_api.get_geometry_mappings(project_id, geometry_id, **kwargs)._embedded
        if len(entity) == 1:
            # print(entity[0].name)
            self.single_entity[key] = entity[0].name
            return self.single_entity[key]
        else:
            raise Exception(f"Found {len(self.single_entity[key])} entities instead of 1: {self.single_entity[key]}")

    def get_entity_names(self, project_id, geometry_id, key, number = None ,**kwargs):
        
        '''
        Retrieve multiple entity IDs of faces/parts. Assume you have multiple parts 
        with the same name. This function will help you out in that case 
        
        Parameters
        ----------
        project_id : str
            The id of the project being worked on.
            
        geometry_id : str
            The geometry that some entities are being retrieved from.  
            
        key : str
            the desired faces/parts for later topological assignment are stored 
            in a dict. each face/part has a specific user defined key, that is 
            passed to this function. In the end, this key would point to the 
            desired entitiy id of a face/part
        
        number: number of faces/parts with the same name that you would like to 
        get the entity ids of 
        
        Returns
        -------
        None.
    
        '''                
        entities = self.geometry_api.get_geometry_mappings(project_id, geometry_id, **kwargs)._embedded
    
        if number is None or len(entities) == number:
            # print(len(entities))
            self.multiple_entities[key] = entities[0].name
            return [self.multiple_entities[key] for e in entities]
        else:
            raise Exception(f"Found {len(self.multiple_entities[key])} entities instead of {number}: {self.multiple_entities[key]}")
    
    """Functions to define the global and initial settings of the simulation"""

    def set_compressible(self, state = False):
        
        '''
        A control function that defines whether a simulation would be set 
        as compressible.
        
        As of now, the compressible simulation settings are not completed. 
        For example, boundary conditions that requires a change of the pressure 
        type when compressible is on, is not yet fixed. 
        
        So until this gets updated keep this toggle False
        
        Parameters
        ----------
        state : boolean
            if true - a compressible setting is passed to spec
             
        Returns
        -------
        None.
    
        '''            
        self.compressible = state
    
    def set_radiation(self, state = False): 
        
        '''
        A control function that defines whether a simulation would include 
        radiation.
        
        As of now, the radiation settings are not completed. 
        For example, boundary conditions that requires additional input 
        when radiation is on, is not yet integrated in this code. 
        
        So until this gets updated keep this toggle False
        
        Parameters
        ----------
        state : boolean
            if true - radiation setting is passed to spec
             
        Returns
        -------
        None.
    
        '''                  
        self.radiation = state
    
    def set_turbulence_model(self, model = "KOMEGASST"):
        
        '''
        choose a turbluence model type 
        Not many choices honestly, it is either "KOMEGASST" or "None"
        Parameters
        ----------
        model : str
            type of turbulencec model
             
        Returns
        -------
        None.
    
        '''      

        self.turbulence_model = model
        
    def set_gravity_direction(self, g_x, g_y, g_z):
        
        '''
        set the gravity direction 
        
        Parameters
        ----------
        g_x : float
            gravity vector in x-direction 
            
        g_y : float
            gravity vector in y-direction 
            
        g_z : float
            gravity vector in z-direction 
                                       
        Returns
        -------
        None.
    
        '''      
        
        self.fluid_model_gravity = sim_sdk.FluidModel(
                gravity=sim_sdk.DimensionalVectorAcceleration(
                    value=sim_sdk.DecimalVector(x= g_x, y= g_y, z= g_z),
                    unit="m"))

            
    def _set_initial_conditions(self, p  = 0 , ux = 0 , uy = 0, uz = 0, t = 19.85, k = 0.00375, omega = 3.375):
        
        '''
        set the initial conditions of the simulation
        
        Parameters
        ----------
        p : float or int 
            initial value of pressure in the domain
            
        ux : float or int 
            initial value of velocity x in the domain
            
        uy : float or int 
            initial value of velocity y in the domain

        uz : float or int 
            initial value of velocity z in the domain
            
        t : float or int 
            initial value of temperature in the domain
            
        k : float or int 
            initial value of kinetic energy in the domain
            
        omega : float or int 
            initial value of specific dissipation rate in the domain
                                                     
        Returns
        -------
        None.
    
        '''             
        self.initial_conditions = sim_sdk.FluidInitialConditions(
    		gauge_pressure_rgh=sim_sdk.DimensionalInitialConditionDomainsPressure(
    			_global= sim_sdk.DimensionalPressure(
    				value= p,
    				unit="Pa",
    			),
    		),
    		velocity= sim_sdk.DimensionalVectorInitialConditionDomainsSpeed(
    			_global= sim_sdk.DimensionalVectorSpeed(
    				value= sim_sdk.DecimalVector(
    					x= ux,
    					y= uy,
    					z= uz,
    				),
    				unit="m/s",
    			),
    			subdomains=[],
    		),
    		temperature= sim_sdk.DimensionalInitialConditionDomainsTemperature(
    			_global= sim_sdk.DimensionalTemperature(
    				value= t,
    				unit="°C",
    			),
    			subdomains=[],
    		),
    		turbulent_kinetic_energy= sim_sdk.DimensionalInitialConditionDomainsTurbulenceKineticEnergy(
    			_global= sim_sdk.DimensionalTurbulenceKineticEnergy(
    				value= k ,
    				unit="m²/s²",
    			),
    		),
    		omega_dissipation_rate= sim_sdk.DimensionalInitialConditionDomainsSpecificTurbulenceDissipationRate(
    			_global= sim_sdk.DimensionalSpecificTurbulenceDissipationRate(
    				value= omega,
    				unit="1/s",
    			),
    		),
    	)
        
        return self.initial_conditions
        
    """Functions to set materials for the simulations by retrieving them from the material library"""

    def add_fluid_material(self, fluid_name = '' ,key_list = []):

        '''
        Add a fluid material to the simulation by retrieving it from the 
        SimScale default material library 
        
        The name provided must match the name of the material in the workbench 
        
        Parameters
        ----------
        fluid_name : str 
            name of the fluid material to be imported  

        key_list : list
        
            the desired faces/parts for later topological assignment are stored 
            in a dict. each face/part has a specific user defined key, that is 
            passed to this function. In the end, this key would point to the 
            desired entitiy id of a face/part.
            
            in this case, you are supplying a list of the keys that points 
            to the parts that are expected to be assigned a certain material 
                           
        Returns
        -------
        None.
    
        '''            
        bodies_to_assign = []
        for key in key_list: 
            bodies_to_assign.append(self.single_entity[key])
            
        material_groups = self.materials_api.get_material_groups().embedded
        default_material_group = next((group for group in material_groups if group.group_type == MaterialGroupType.SIMSCALE_DEFAULT), None)
        if not default_material_group:
            raise Exception(f"Couldn't find default material group in {material_groups}")
        
        default_materials = self.materials_api.get_materials(material_group_id=default_material_group.material_group_id).embedded
        material_of_interest = next((material for material in default_materials if material.name == fluid_name), None)
        if not material_of_interest:
            raise Exception("Couldn't find default {m} material in {dm}".format(m = fluid_name, dm = default_materials ))
            
        material_data = self.materials_api.get_material_data(
            material_group_id=default_material_group.material_group_id,
            material_id=material_of_interest.id
        )
        material_update_request = MaterialUpdateRequest(
            operations=[
                MaterialUpdateOperation(
                    path="/materials/fluids",
                    material_data=material_data,
                    reference=MaterialUpdateOperationReference(
                        material_group_id=default_material_group.material_group_id,
                        material_id=material_of_interest.id
                    )
                )
            ]
        )
        self.fluid_material.append( self.simulation_api.update_simulation_materials(self.project_id, self.simulation_id, material_update_request))
        
        # Add assignments to the new material
        self.simulation_spec = self.simulation_api.get_simulation(self.project_id, self.simulation_id)
        self.simulation_spec.model.materials.fluids[0].topological_reference = sim_sdk.TopologicalReference(entities= bodies_to_assign)
        self.simulation_api.update_simulation(self.project_id, self.simulation_id, self.simulation_spec)
            
    def add_solid_material(self, solid_name = '' ,key_list = []):
        
        '''
        Add a solid material to the simulation by retrieving it from the 
        SimScale default material library 
        
        The name provided must match the name of the material in the workbench 
        
        Parameters
        ----------
        solid_name : str 
            name of the solid material to be imported  

        key_list : list
        
            the desired faces/parts for later topological assignment are stored 
            in a dict. each face/part has a specific user defined key, that is 
            passed to this function. In the end, this key would point to the 
            desired entitiy id of a face/part.
            
            in this case, you are supplying a list of the keys that points 
            to the parts that are expected to be assigned a certain material 
                           
        Returns
        -------
        None.
    
        ''' 
        bodies_to_assign = []
        for key in key_list: 
            bodies_to_assign.append(self.single_entity[key])
            
        material_groups = self.materials_api.get_material_groups().embedded
        default_material_group = next((group for group in material_groups if group.group_type == MaterialGroupType.SIMSCALE_DEFAULT), None)
        if not default_material_group:
            raise Exception(f"Couldn't find default material group in {material_groups}")
        
        default_materials = self.materials_api.get_materials(material_group_id=default_material_group.material_group_id).embedded
        material_of_interest = next((material for material in default_materials if material.name == solid_name), None)
        if not material_of_interest:
            raise Exception("Couldn't find default {m} material in {dm}".format(m = solid_name, dm = default_materials ))
        
        material_data = self.materials_api.get_material_data(
            material_group_id=default_material_group.material_group_id,
            material_id=material_of_interest.id
        )
        material_update_request = MaterialUpdateRequest(
            operations=[
                MaterialUpdateOperation(
                    path="/materials/solids",
                    material_data=material_data,
                    reference=MaterialUpdateOperationReference(
                        material_group_id=default_material_group.material_group_id,
                        material_id=material_of_interest.id
                    )
                )
            ]
        )
        self.solid_material.append(self.simulation_api.update_simulation_materials(
            self.project_id, self.simulation_id, material_update_request))
        
        # Add assignments to the new material
        self.simulation_spec = self.simulation_api.get_simulation(self.project_id, self.simulation_id)
        self.simulation_spec.model.materials.solids[0].topological_reference = sim_sdk.TopologicalReference(entities= bodies_to_assign)
        self.simulation_api.update_simulation(self.project_id, self.simulation_id, self.simulation_spec)
        
    """Function to define simulation numerics"""

    def _set_simulation_numerics(self):
    
        '''
        function that allows changing the simulation numerics. 
        
        As of now the settings are not exposed. Hence, any numerics modifications 
        would need to take place directly in this function and the main script. 
        
        Parameters
        ----------
     
        Returns
        -------
        None.
    
        '''         
        self.fluid_numerics = sim_sdk.FluidNumerics(
                relaxation_type="MANUAL",
                relaxation_factor=sim_sdk.RelaxationFactor(
                        pressure_rgh_field=0.7,
                        velocity_equation=0.3,
                        temperature_equation=0.8,
                ),
                num_non_orthogonal_correctors=1,
                solvers=sim_sdk.FluidSolvers(
                        velocity_solver=sim_sdk.PBICGSolver(
                                type="PBICG",
                                absolute_tolerance=1.0E-15,
                                relative_tolerance=0.01,
                                preconditioner=sim_sdk.DILUPreconditioner(
                                        type="DILU",
                                ),
                        ),
                        temperature_solver=sim_sdk.PBICGSolver(
                                type="PBICG",
                                absolute_tolerance=1.0E-15,
                                relative_tolerance=0.01,
                                preconditioner=sim_sdk.ILUCpPreconditioner(
                                        type="ILUCP",
                                        fill_in_level=1,
                                ),
                        ),
                        pressure_rgh_solver=sim_sdk.ILUCpPreconditioner(
                                type="GAMG",
                                #absolute_tolerance=1.0E-15,
                                #relative_tolerance=0.01,
                                # smoother="GAUSSSEIDEL",
                                # num_pre_sweeps=1,
                                # num_post_sweeps=1,
                                # cache_agglomeration_on=True,
                                # num_cells_coarsest_level=100,
                                # num_merge_levels=1,
                        ),
                ),
                schemes=sim_sdk.Schemes(
                        second_order_convection=False,
                ),

            )
        
        return self.fluid_numerics
        
    """Functions to define geometry primitives"""
    
    def set_single_geometry_primitive_point(self, name, pos_x, pos_y, pos_z):
        
        '''
        define a point geometry primitive 
        
        Parameters
        ----------
        name : str 
            name of the defined point geometry primitive 

        pos_x : float or int
            
            x-coordinate of the point        

        pos_y : float or int
            
            y-coordinate of the point      

        pos_z : float or int
            
            z-coordinate of the point                  
        Returns
        -------
        None.
    
        ''' 
        geometry_primitive_point = sim_sdk.Point(
            name= name,
            center=sim_sdk.DimensionalVectorLength(value=sim_sdk.DecimalVector(x=pos_x, y=pos_y, z=pos_z), unit="m"),
        )
        self.geometry_primitive_uuid = self.simulation_api.create_geometry_primitive(
            self.project_id, geometry_primitive_point).geometry_primitive_id
        
        print(f"geometry_primitive_uuid: {self.geometry_primitive_uuid}")
    
    def set_multiple_geometry_primitive_points(self, path_to_csv):
        
        '''
        define multiple point geometry primitives at once using a .csv file
        
        csv format: Label X Y Z (delimiter is a space)

        Parameters
        ----------
        path_to_csv : Pathlib path 
        
            name of the defined point geometry primitive 
            
            ex: pathlib.Path().cwd() / "probe_points" / "probe_list.txt"             
    
        Returns
        -------
        None.
    
        '''         
        
        probe_pd = pd.read_csv( path_to_csv , header = 0, delimiter = ' ') 
        probe_pd.set_index("Label")
        for row in probe_pd.index:
            pos_x = probe_pd['X'][row]
            pos_y = probe_pd['Y'][row]
            pos_z = probe_pd['Z'][row]
            print(pos_x)
            print(pos_y)
            print(pos_z)
        
            geometry_primitive_point = sim_sdk.Point(
                name= "point{}".format(row),
                center=sim_sdk.DimensionalVectorLength(value=sim_sdk.DecimalVector(x=float(pos_x), y=float(pos_y), z=float(pos_z)), unit="m"),
            )
            self.geometry_primitive_uuid_list.append(self.simulation_api.create_geometry_primitive(
                self.project_id, geometry_primitive_point).geometry_primitive_id)
            
        print(f"geometry_primitive_uuid: {self.geometry_primitive_uuid}")
        
        
    def set_cylinderical_geometry_primitive(self, r, pos_x, pos_y, pos_z, len_x, len_y, len_z, name = 'cylinder' ,unit = 'm') : 
        
        '''
        define a cylinderical geometry primitive 
        
        Parameters
        ----------
        r : float or int 
        
            radius of the cylinder 
            
        pos_x: float or int 
            
            The position of the cylinder in x-direction in reference to origin (0,0,0) 
            
        pos_y: float or int 
            
            The position of the cylinder in y-direction in reference to origin (0,0,0) 
            
        pos_z: float or int 
            
            The position of the cylinder in z-direction in reference to origin (0,0,0) 
            
        len_x: float or int 
            
            length of the cylinder in the x-direction 
            
        len_y: float or int 
            
            length of the cylinder in the y-direction 
            
        len_z: float or int 
            
            length of the cylinder in the z-direction 
            
        name: str
            name of the cylinderical geometry primitive
                        
        Returns
        -------
        geometry_primitive_body_uuid: str 
            
            A pointer to the id of the geometry primitive, 
            To be used when assigining it to a topological entity
        '''                   
        geometry_primitive_body = sim_sdk.GeometryPrimitive(type = 'CYLINDER',
                                      name = name,
                                      reference = sim_sdk.DimensionalVectorLength(value= sim_sdk.DecimalVector(x=float(pos_x),
                                                                                                               y=float(pos_y), 
                                                                                                               z=float(pos_z)),
                                        unit="m"), 
                                      axis = sim_sdk.DimensionalVectorLength(value= sim_sdk.DecimalVector(x=float(len_x),
                                                                                                          y=float(len_y),
                                                                                                          z=float(len_z)),
                                        unit="m"), 
                                      radius = sim_sdk.DimensionalLength(value = r, unit = 'm')
                                      
                ) 
        
        geometry_primitive_body_uuid = self.simulation_api.create_geometry_primitive(
            self.project_id, geometry_primitive_body).geometry_primitive_id
        
        print(f"geometry_primitive_uuid: {geometry_primitive_body_uuid}")
    
        return geometry_primitive_body_uuid
    
    
    
    def set_cartesian_box_geometry_primitive(self, min_x, min_y, min_z, max_x, max_y, max_z, name = "box", unit = 'm'):
        
        '''
        define a cartesian box geometry primitive 
        
        Parameters
        ----------
        
        min_x: float or int 
            
            The minimum x-coordinate of the box in reference to origin (0,0,0) 
            
        min_y: float or int 
            
            The minimum y-coordinate of the box in reference to origin (0,0,0) 
            
        min_z: float or int 
            
            The minimum z-coordinate of the box in reference to origin (0,0,0) 
            
        max_x: float or int 
            
            The maximum x-coordinate of the box in reference to origin (0,0,0) 
            
        max_y: float or int 
            
            The maximum y-coordinate of the box in reference to origin (0,0,0) 
            
        max_z: float or int 
            
            The maximum z-coordinate of the box in reference to origin (0,0,0) 
            
        name: str
            name of the cylinderical geometry primitive
                        
        Returns
        -------
        geometry_primitive_body_uuid: str 
            
            A pointer to the id of the geometry primitive, 
            To be used when assigining it to a topological entity
    
        '''  
        if min_x > max_x or min_y > max_y or min_z > max_z : 
            raise Exception("The minumum values of the box must be smaller than the maximums")

        geometry_primitive_body = sim_sdk.GeometryPrimitive(type = 'CARTESIAN_BOX',
                                      name = name,
                                      min = sim_sdk.DimensionalVectorLength(value= sim_sdk.DecimalVector(x=float(min_x),
                                                                                                               y=float(min_y),
                                                                                                               z=float(min_z)),
                                                                                                                unit="m"), 
                                      max = sim_sdk.DimensionalVectorLength(value= sim_sdk.DecimalVector(x=float(max_x),
                                                                                                          y=float(max_y),
                                                                                                          z=float(max_z)),
                                                                                                                unit = 'm')
                                      ) 
        
        geometry_primitive_body_uuid = self.simulation_api.create_geometry_primitive(
            self.project_id, geometry_primitive_body).geometry_primitive_id
        
        print(f"geometry_primitive_uuid: {geometry_primitive_body_uuid}")
    
        return geometry_primitive_body_uuid
            
    """Functions to define the boundary conditions of the simulation"""

    def constant_velocity_inlet_bc(self, speed_x, speed_y, speed_z, temp,  name = "velocity inlet", key = ' '):
        
        '''
        define a constant velocity inlet boundary condition 
        
        Parameters
        ----------
        
        speed_x: float or int 
            
            The speed of the flow in x-direction 
            
        speed_y: float or int 
            
            The speed of the flow in y-direction 
            
        speed_z: float or int 
            
            The speed of the flow in z-direction 
            
        temp: float or int 
            
            The temperature of the fluid in celsius degrees 
            
        name: str
            name of boundary condition
            
        key : str
            the desired faces/parts for later topological assignment are stored 
            in a dict. each face/part has a specific user defined key, that is 
            passed to this function. In the end, this key would point to the 
            desired entitiy id of a face/part
         
                        
        Returns
        -------
        None. 
            
        '''          
        self.velocity_inlet.append(
                sim_sdk.VelocityInletBC(
                name= name,
                velocity= sim_sdk.FixedValueVBC(
                value= sim_sdk.DimensionalVectorFunctionSpeed(
                    value= sim_sdk.ComponentVectorFunction(
                        x= sim_sdk.ConstantFunction(value= speed_x), 
                        y= sim_sdk.ConstantFunction(value=speed_y),
                        z= sim_sdk.ConstantFunction(value= speed_z)
                        )
                    )
                ),                
            temperature= sim_sdk.FixedValueTBC(
                value=DimensionalFunctionTemperature(value=sim_sdk.ConstantFunction(value= temp), unit="°C")
            ),        
            topological_reference=TopologicalReference(entities=[self.single_entity[key]]),
            ))


    def flow_rate_velocity_inlet(self, flowrate, temp, type, name = "Velocity Inlet", key = ' '):
        
        '''
        define a flow rate velocity inlet boundary condition 
        
        depending on the type specified: either a mass flow or a volumetric flow rate is defined 
        
        Parameters
        ----------
        
        flowrate: float or int 
            
            The value of the flowrate
            
        temp: float or int 
            
            The temperature of the fluid in celsius degrees 

        type: str 
            
            to decide on whether a mass flowrate of a volumetric flowrate to be defined 
            
            if you would like to setup a volumetric flowrate then pass "volumetric" to type
                        
        name: str
            name of boundary condition
            
        key : str
            the desired faces/parts for later topological assignment are stored 
            in a dict. each face/part has a specific user defined key, that is 
            passed to this function. In the end, this key would point to the 
            desired entitiy id of a face/part
         
                        
        Returns
        -------
        None. 
            
        '''             
        if type == 'volumetric':
            
            self.velocity_inlet.append(
                sim_sdk.VelocityInletBC(
                    name=name,
                    velocity=sim_sdk.FlowRateInletVBC(flow_rate=sim_sdk.VolumetricFlow(
                        type="VOLUMETRIC",
                        value=sim_sdk.DimensionalFunctionVolumetricFlowRate(
                            value=sim_sdk.ConstantFunction(
                                type="CONSTANT",
                                value=flowrate,
                            ),
                            unit="m³/s",
                        ),
                    ),),
                    temperature=sim_sdk.FixedValueTBC(
                        value=DimensionalFunctionTemperature(
                            value=sim_sdk.ConstantFunction(value=temp), unit="°C")
                    ),
                    topological_reference=TopologicalReference(
                        entities=[self.single_entity[key]]),
                ))

        else:
        
            self.velocity_inlet.append(
                sim_sdk.VelocityInletBC(
                    name=name,
                    velocity=sim_sdk.FlowRateInletVBC(flow_rate=sim_sdk.MassFlow(
                        type="MASS",
                        value=sim_sdk.DimensionalFunctionMassFlowRate(
                            value=sim_sdk.ConstantFunction(
                                type="CONSTANT",
                                value=flowrate,
                            ),
                            unit="kg/s",
                        ),
                    ),),
                    temperature=sim_sdk.FixedValueTBC(
                        value=DimensionalFunctionTemperature(
                            value=sim_sdk.ConstantFunction(value=temp), unit="°C")
                    ),
                    topological_reference=TopologicalReference(
                        entities=[self.single_entity[key]]),
                ))


    def natural_convection_inlet_outlet_bc(self, pressure, temp, name = "Natural Convection Inlet Outlet", key_list = []):
        
        '''
        define a natural convection inlet outlet boundary condition 
                
        Parameters
        ----------
        
        pressure: float or int 
            
            Depending on whether compressible is True or Flase. Pass either, a
            gauge pressure value or absolute pressure value 
            
        temp: float or int 
            
            The temperature of the fluid in celsius degrees 

        name: str
            name of boundary condition
            
        key_list : list
        
            the desired faces/parts for later topological assignment are stored 
            in a dict. each face/part has a specific user defined key, that is 
            passed to this function. In the end, this key would point to the 
            desired entitiy id of a face/part.
            
            in this case, you are supplying a list of the keys that points 
            to the parts that are expected to be assigned a certain material 
           
        Returns
        -------
        None. 
            
        '''      
        faces_to_assign = []
        for key in key_list: 
            faces_to_assign.append(self.single_entity[key])
        
        if self.compressible: 
                    
            self.natural_convection.append(
                sim_sdk.NaturalConvectionInletOutletBC(
                    name= name,
                    type="NATURAL_CONVECTION_INLET_OUTLET",
        			pressure_rgh=sim_sdk.AmbientPBC(
        				type="AMBIENT_PRESSURE",
        				ambient_pressure=sim_sdk.DimensionalPressure(
        					value=pressure,
        					unit="Pa",
        				),
        			),
        			temperature=sim_sdk.AmbientTBC(
        				type="AMBIENT_TEMPERATURE",
        				ambient_temperature=sim_sdk.DimensionalTemperature(
        					value=19.85,
        					unit="°C",
        				),
        			),
                    topological_reference=sim_sdk.TopologicalReference(
                        entities= faces_to_assign
                    ),
                ))
            
        else: 
            
            self.natural_convection.append(
                sim_sdk.NaturalConvectionInletOutletBC(
                    name= name,
                    type="NATURAL_CONVECTION_INLET_OUTLET",
                    gauge_pressure_rgh=sim_sdk.AmbientPBC(
                        type="AMBIENT_PRESSURE",
                        ambient_pressure=sim_sdk.DimensionalPressure(
                            value= pressure,
                            unit="Pa",
                        ),
                    ),
                    temperature=sim_sdk.AmbientTBC(
                        type="AMBIENT_TEMPERATURE",
                        ambient_temperature=sim_sdk.DimensionalTemperature(
                            value= temp,
                            unit="°C",
                        ),
                    ),
                    topological_reference=sim_sdk.TopologicalReference(
                        entities= faces_to_assign
                    ),
                ))
                       
    
    def pressure_inlet_bc(self, value, temp ,name = "Pressure inlet" , unit = 'pa', key = ''):
        
        '''
        define a pressure inlet boundary condition 
        
        Parameters
        ----------
        
        value: float or int 
            
            The value of the inlet pressure depending on the unit
            
        temp: float or int 
            
            The temperature of the fluid in celsius degrees 

        name: str
            name of boundary condition
            
        key : str
            the desired faces/parts for later topological assignment are stored 
            in a dict. each face/part has a specific user defined key, that is 
            passed to this function. In the end, this key would point to the 
            desired entitiy id of a face/part
         
                        
        Returns
        -------
        None. 
            
        '''              
        self.pressure_inlet.append(sim_sdk.PressureInletBC(
            name= name,
            gauge_pressure_rgh= sim_sdk.TotalPBC(
                type = "TOTAL_PRESSURE", total_pressure = sim_sdk.DimensionalFunctionPressure(value = sim_sdk.ConstantFunction(value= value), unit = unit))
            ,            
            temperature= sim_sdk.FixedValueTBC(
                value=DimensionalFunctionTemperature(value=sim_sdk.ConstantFunction(value= temp), unit="°C")
            ),
            topological_reference=TopologicalReference(entities=[self.single_entity[key]]),
        ))   
        
    def pressure_outlet_bc(self, value = 0 ,name = "Pressure outlet" , unit = 'pa', key = ''):
        
        '''
        define a pressure outlet boundary condition 
        
        Parameters
        ----------
        
        value: float or int 
            
            The value of the inlet pressure depending on the unit
            
        temp: float or int 
            
            The temperature of the fluid in celsius degrees 

        name: str
            name of boundary condition
            
        key : str
            the desired faces/parts for later topological assignment are stored 
            in a dict. each face/part has a specific user defined key, that is 
            passed to this function. In the end, this key would point to the 
            desired entitiy id of a face/part
         
                        
        Returns
        -------
        None. 
            
        '''              
        self.pressure_outlet.append(sim_sdk.PressureOutletBC(
            name= name,
            gauge_pressure_rgh= sim_sdk.FixedValuePBC(
                value= sim_sdk.DimensionalFunctionPressure(value= sim_sdk.ConstantFunction(value= value ), unit= unit)
            ),
            topological_reference=sim_sdk.TopologicalReference(entities=[self.single_entity[key]]),
        ))
                
    def no_slip_fixed_temp_wall_bc(self, temp, name = "Fixed Temp" ,key_list = []):
        
        '''
        define a no-slip fixed wall temperature boundary condition 
        
        Parameters
        ----------
        
        temp: float or int 
            
            The temperature of the wall in celsius degrees 

        name: str
            name of boundary condition
            
        key_list : list
        
            the desired faces/parts for later topological assignment are stored 
            in a dict. each face/part has a specific user defined key, that is 
            passed to this function. In the end, this key would point to the 
            desired entitiy id of a face/part.
            
            in this case, you are supplying a list of the keys that points 
            to the parts that are expected to be assigned a certain material 
        
        Returns
        -------
        None. 
            
        '''               
        faces_to_assign = []
        for key in key_list: 
            # print(key)
            # print(self.single_entity[key])
            faces_to_assign.append(self.single_entity[key])
            # print(faces_to_assign)

        self.wall_bc.append(           
            sim_sdk.WallBC(
                    name= name,
                    velocity=sim_sdk.NoSlipVBC(),
                    temperature=sim_sdk.FixedValueTBC(
                        value=sim_sdk.DimensionalFunctionTemperature(value=sim_sdk.ConstantFunction(value=temp), unit="°C")
                    ),
                    topological_reference=sim_sdk.TopologicalReference(entities = faces_to_assign),
                ))

        
    def external_wall_heat_flux_bc(self, method = "DERIVED", amb_temp = None, 
                                   power = None , htc = None , heat_flux = None,
                                   name = "stuff_of_the_heat_flux" ,
                                   key_list = []):
        
        #for later: expose the option to add an wall thermal layer and an 
        #additional heat flux 
        
        '''
        define an external wall heat flux boundary condition 
        choose from three different types: DERIVED, FIXED, FIXED_POWER
        For DERIVED     --> htc and amb_temp 
        for FIXED       --> heat_flux
        for FIXED_POWER --> power
        
        Parameters
        ----------
        
        method: str 
            
            to decide on what type of this boundary condition to be used 
            DERIVED, FIXED, FIXED_POWER
            
        amb_temp: float or int 
            
            The temperature of the wall in celsius degrees 

        power: float or int 
            
            The value of dissipated power in watts  to be assigned 
            
        htc: float or int 
        
            the heat transfer coefficient value
            
        heat_flux: float or int 
        
            the value of the heat flux on the assigned surface
        
        name: str
            name of boundary condition        
            
        key_list : list
        
            the desired faces/parts for later topological assignment are stored 
            in a dict. each face/part has a specific user defined key, that is 
            passed to this function. In the end, this key would point to the 
            desired entitiy id of a face/part.
            
            in this case, you are supplying a list of the keys that points 
            to the parts that are expected to be assigned a certain material 
        
        Returns
        -------
        None. 
            
        '''                
        faces_to_assign = []
        for key in key_list:
            # print(key)
            # print(self.single_entity[key])
            faces_to_assign.append(self.single_entity[key])
            # print(faces_to_assign)

        if method == "DERIVED":
            #for later: add code that allows incorportaing wall thermals into the simulation
            if amb_temp and htc != None:
                self.wall_bc.append(
                    sim_sdk.WallBC(
                        name=name,
                        velocity=NoSlipVBC(),
                        temperature=sim_sdk.ExternalWallHeatFluxTBC(
                            heat_flux=sim_sdk.DerivedHeatFlux(
                                type="DERIVED",
                                heat_transfer_coefficient=sim_sdk.DimensionalThermalTransmittance(
                                    value=htc, unit='W/(K·m²)'),
                                ambient_temperature=sim_sdk.DimensionalTemperature(
                                    value=amb_temp, unit="°C"),
                                additional_heat_flux=None,
                                wall_thermal=None)
                        ),
                        topological_reference=TopologicalReference(
                            entities=faces_to_assign),
                    ))
            else:
               raise Exception(
                   "Provide the ambient temperature and heat transfer coefficient values")

        elif method == "FIXED":
            if heat_flux != None:
                self.wall_bc.append(
                    sim_sdk.WallBC(
                        name=name,
                        velocity=NoSlipVBC(),
                        temperature=sim_sdk.ExternalWallHeatFluxTBC(
                            heat_flux=sim_sdk.FixedHeatFlux(
                                type="FIXED",
                                function=sim_sdk.DimensionalFunctionHeatFlux(value=sim_sdk.ConstantFunction(
                                    type="CONSTANT",
                                    value=heat_flux,
                                ), unit="W/m²"))
                        ),
                        topological_reference=TopologicalReference(
                            entities=faces_to_assign),
                    ))
            else:
                raise Exception("Provide the heat flux value")

        elif method == "FIXED_POWER":
            if power != None:
                self.wall_bc.append(
                    sim_sdk.WallBC(
                        name=name,
                        velocity=NoSlipVBC(),
                        temperature=sim_sdk.ExternalWallHeatFluxTBC(
                            heat_flux=sim_sdk.FixedPowerHeatFlux(
                                type="FIXED_POWER",
                                function=sim_sdk.DimensionalFunctionPower(
                                    value=sim_sdk.ConstantFunction(
                                        type="CONSTANT",
                                        value=power,
                                    ),
                                    unit="W",
                                ),
                            )
                        ),
                        topological_reference=TopologicalReference(
                            entities=faces_to_assign),
                    ))
            else:
                raise Exception("Provide the power value")

    def set_boundary_conditions(self):
        #for later: add a validation against multiple face assignments for a BC
        #for later: think of how you can make this code better looking (one liners?)
        
        '''
        store all of the pre-defined boundary condtions in one variable to pass
        that variable in a clean way to the simulation spec 
        
        just an organizational function
        
        Parameters
        ----------
        None. 
        
        Returns
        -------
        None. 
            
        '''                      
        if len(self.velocity_inlet) >= 1 : 
            for v_inlet in self.velocity_inlet :
                if v_inlet in self.boundary_conditions:
                    continue
                else:
                    self.boundary_conditions.append(v_inlet)

        if len(self.velocity_outlet) >= 1 : 
            for v_outlet in self.velocity_outlet: 
                if v_outlet in self.boundary_conditions:
                    continue
                else:                
                    self.boundary_conditions.append(v_outlet)
                
        if len(self.pressure_inlet) >= 1 : 
            for p_inlet in self.pressure_inlet: 
                if p_inlet in self.boundary_conditions:
                    continue
                else:                
                    self.boundary_conditions.append(p_inlet)
                
        if len(self.pressure_outlet) >= 1 : 
            for p_outlet in self.pressure_outlet:
                if p_outlet in self.boundary_conditions:
                    continue
                else:                
                    self.boundary_conditions.append(p_outlet)
                    
        if len(self.natural_convection) >= 1 : 
            for natural_inlet_outlet in self.natural_convection:
                if natural_inlet_outlet in self.boundary_conditions:
                    continue
                else:                
                    self.boundary_conditions.append(natural_inlet_outlet)
                              
        if len(self.wall_bc) >= 1 : 
            for wall in self.wall_bc: 
                if wall in self.boundary_conditions:
                    continue
                else:                
                    self.boundary_conditions.append(wall)
      
    """Functions to define advanced concepts"""
    
    def set_power_sources(self ,power , name = "power source", method = "ABSOLUTE" , key_list = [], geo_primitive_ids = []):
        #for later: add code to assign bodies that has multiple parts using get entities(for example one assignment for all batteries )
        #for later: think how you can add a validation rule against exact same power source definitions; is this needed?)

        '''
        define a power source on a body/part or geometry primtive using either an absolute power 
        source [w] or a specific one [w/m3] 
        
        Parameters
        ----------
        
        power: float or int 
            
            The value of dissipated power from a particular part 

        name: str
            name of boundary condition
            
        method : str 
        
            A choice of either "ABSOLUTE" or "SPECIFIC "       
        
        key_list : list
        
            the desired faces/parts for later topological assignment are stored 
            in a dict. each face/part has a specific user defined key, that is 
            passed to this function. In the end, this key would point to the 
            desired entitiy id of a face/part.
            
            in this case, you are supplying a list of the keys that points 
            to the parts that are expected to be assigned a certain material 
        
        geo_primitive_ids: list
            
            a list containing all the pointers to the geometry primitives to be used for assignment 
            (a list with all geometry primitive uuids is what you need to pass here)
                        
        Returns
        -------
        None. 
            
        '''     
        #assign parts         
        bodies_to_assign = []
        for key in key_list:
            # print(key)
            # print(self.single_entity[key])
            bodies_to_assign.append(self.single_entity[key])
            # print(bodies_to_assign)

        #assign geometry primitives
        geometry_primitives_to_assign = []
        for primitive_uuid in geo_primitive_ids:
            
            geometry_primitives_to_assign.append(primitive_uuid)

        if method == "ABSOLUTE":
            self.power_sources.append(
                sim_sdk.AbsolutePowerSource(
                    name=name,
                    heat_flux=DimensionalFunctionPower(
                        value=sim_sdk.ConstantFunction(value=power), unit="W"),
                    topological_reference=TopologicalReference(
                        entities=bodies_to_assign),
                    geometry_primitive_uuids= geometry_primitives_to_assign),)

        elif method == "SPECIFIC":
            self.power_sources.append(
                sim_sdk.SpecificPowerSource(
                    name=name,
                    heat_flux=sim_sdk.DimensionalFunctionVolumetricPower(
                        value=sim_sdk.ConstantFunction(value=power), unit="W/m³"),
                    topological_reference=TopologicalReference(
                        entities=bodies_to_assign),

                    geometry_primitive_uuids= geometry_primitives_to_assign))

    def set_porous_media(self):
        
        '''
        Not yet implemented... 
        
        define a porous media on a body/part or a geometry primitive
        
        Parameters
        ----------
                   
        Returns
        -------
        None. 
            
        '''        
        pass
    
    def set_momentum_sources(self, speed_x, speed_y, speed_z, name = "momentum source", key_list = [], geo_primitive_ids = [] ):
        
        '''
        define a momentum source on a body/part using either an absolute power 
        source [w] or a specific one [w/m3] 
        
        Parameters
        ----------
        
        speed_x: float or int 
            
            The speed of the flow in x-direction 

        speed_y: float or int 
            
            The speed of the flow in y-direction 

        speed_z: float or int 
            
            The speed of the flow in z-direction 

        name: str
            name of boundary condition
            
        method : str 
        
            A choice of either "ABSOLUTE" or "SPECIFIC "       
        
        key_list : list
        
            the desired faces/parts for later topological assignment are stored 
            in a dict. each face/part has a specific user defined key, that is 
            passed to this function. In the end, this key would point to the 
            desired entitiy id of a face/part.
            
            in this case, you are supplying a list of the keys that points 
            to the parts that are expected to be assigned a certain material 
            
        geo_primitive_ids: str 
            
            a pointer to the geometry primitive to be used for assignment 
            (geometry primitive uuid is what you need to pass here)
                    
                        
        Returns
        -------
        None. 
            
        '''         
        #assign parts
        bodies_to_assign = []
        for key in key_list:
            bodies_to_assign.append(self.single_entity[key])
            
        #assign geometry primitives
        geometry_primitives_to_assign = []
        for primitive_uuid in geo_primitive_ids:
            
            geometry_primitives_to_assign.append(primitive_uuid)

        self.momentum_sources.append(
            sim_sdk.AverageVelocityMomentumSource(
                type="AVERAGE_VELOCITY",
                name=name,
                average_velocity=sim_sdk.DimensionalVectorSpeed(
                    value=sim_sdk.DecimalVector(
                        x=speed_x,
                        y=speed_y,
                        z=speed_z,
                    ),
                    unit="m/s",
                ),

                topological_reference=TopologicalReference(
                    entities=bodies_to_assign),

            				geometry_primitive_uuids= geometry_primitives_to_assign
            ))

    def set_thermal_resistance_network(self):
        '''
        Not yet implemented... 
        
        define a porous media on a body/part or a geometry primitive
        
        Parameters
        ----------
                   
        Returns
        -------
        None. 
            
        '''        
        pass

    def set_advanced_concepts(self):
        #for later: add a validation against multiple face assignments for a BC
        # (this validation needs to be in the attribute functions e.g set_power_sources)
        
        '''
        store all of the pre-defined advanced concepts in one variable to pass
        that variable in a clean way to the simulation spec 
        
        just an organizational function
        
        Parameters
        ----------
        None. 
        
        Returns
        -------
        None. 
            
        '''           
        self.advanced_concepts = sim_sdk.AdvancedConcepts(
            power_sources    = self.power_sources, 
            porous_mediums   = self.porous_media, 
            momentum_sources = self.momentum_sources,
            thermal_resistance_networks = self.thermal_resistance_networks)

    """Functions to define simulation control settings of the simulation"""

    def set_simulation_end_time(self, time = 1000):
        
        '''
        define the end time of the simulation - this is steady state so the 
        total number of iterations 
        
        (transient would be a future implementation step)
        
        Parameters
        ----------
        
        time: int 
            
            The total number of iterations in the simulation
           
        Returns
        -------
        None. 
            
        '''             
        self.end_time = sim_sdk.DimensionalTime(value= time, unit="s")
        
    def set_simulation_time_step(self, time_step = 1): 
        
        '''
        define the timestep of the simulation - this is steady state so the 
        interval between each iteration
        
        (transient would be a future implementation step)
        
        Parameters
        ----------
        
        timestep: int 
            
            The iteration interval - default is 1 
           
        Returns
        -------
        None. 
            
        '''           
        self.delta_t = sim_sdk.DimensionalTime(value= time_step, unit="s")
        
    def set_simulation_write_controls(self, write_interval = 1000): 
        
        '''
        define the write interval of the results - this is steady state so 
        we are only interested in the final results (iteration)
        
        Have this write_interval = end_time of the simulation 
        
        Parameters
        ----------
        
        write_interval: int 
            
            how frequently are results written to be visualized in the post-processor
            keep this the same as the simulation end time for steady state 
            
        Returns
        -------
        None. 
            
        '''           
        self.write_control = sim_sdk.TimeStepWriteControl(write_interval= write_interval)
    
    def set_simulation_max_run_time(self, max_run_time = 40000):
        
        '''
        define the maximum time allowed for the simulation to run 
                
        Parameters
        ----------
        
        max_run_time: int 
            
            don't make this number extremely large nor way too small
            
        Returns
        -------
        None. 
            
        '''             
        self.max_run_time = sim_sdk.DimensionalTime(value= max_run_time, unit="s")
    
    def set_simulation_control(self):
        
        '''
        store all of the pre-defined simulation control settings in one variable to pass
        that variable in a clean way to the simulation spec 
        
        just an organizational function
        
        Parameters
        ----------
        None. 
        
        Returns
        -------
        None. 
            
        '''            
        self.simulation_control = sim_sdk.FluidSimulationControl(
            
                end_time= self.end_time,
                delta_t=  self.delta_t,
                write_control= self.write_control,
                max_run_time= self.max_run_time,
                decompose_algorithm=sim_sdk.ScotchDecomposeAlgorithm(),)
    
    """Functions to define additional result export"""
    
    def set_area_averages(self, name = "Area average", write_interval = 10, key_list = []):
        
        '''
        set an area average on a face. useful to monitor convergence of a simulation  
        
        Parameters
        ----------
        name: str 
            
            name of the area average result control 
            
        write_interval : int 
        
            how often are the results updated through out the simulation 
            
        key_list : list
        
            the desired faces/parts for later topological assignment are stored 
            in a dict. each face/part has a specific user defined key, that is 
            passed to this function. In the end, this key would point to the 
            desired entitiy id of a face/part.
            
            in this case, you are supplying a list of the keys that points 
            to the parts that are expected to be assigned a certain material 
                    
        Returns
        -------
        None. 
            
        '''   
         
        faces_to_assign = []
        for key in key_list: 
            faces_to_assign.append(self.single_entity[key])

        self.surface_data.append(
            sim_sdk.AreaAverageResultControl(
                name = name, 
                write_control = sim_sdk.TimeStepWriteControl(write_interval = write_interval),
                topological_reference = sim_sdk.TopologicalReference(entities = faces_to_assign),))


    def set_area_volumes(self, name = "Volume average", write_interval = 10 , key_list = []):
        
        '''
        set an area volume on a face. useful to monitor convergence of a simulation  
        
        Parameters
        ----------
        name: str 
            
            name of the area volume result control 
            
        write_interval : int 
        
            how often are the results updated through out the simulation 
            
        key_list : list
        
            the desired faces/parts for later topological assignment are stored 
            in a dict. each face/part has a specific user defined key, that is 
            passed to this function. In the end, this key would point to the 
            desired entitiy id of a face/part.
            
            in this case, you are supplying a list of the keys that points 
            to the parts that are expected to be assigned a certain material 
                    
        Returns
        -------
        None. 
            
        '''           
        faces_to_assign = []
        for key in key_list: 
            faces_to_assign.append(self.single_entity[key])

        self.surface_data.append(
            sim_sdk.AreaIntegralResultControl(
                name = name, 
                write_control = sim_sdk.TimeStepWriteControl(write_interval = write_interval),
                topological_reference = sim_sdk.TopologicalReference(entities = faces_to_assign),))

    
    def set_probe_points(self, name, source = "single"):
        
        '''
        set a probe point(s) result control item
        
        Parameters
        ----------
        name: str 
            
            name of the area volume result control 
            
        source : str 
        
            either 'single' or 'multiple'. just to tell the code if you are trying 
            to define a single probe point item or multiple probe points at once 
             
        Returns
        -------
        None. 
            
        '''            
        if source == "single":
            
            self.probe_points.append(
                sim_sdk.ProbePointsResultControl(
                    name= name,
                    write_control=sim_sdk.TimeStepWriteControl(write_interval=10),
                    geometry_primitive_uuids=[self.geometry_primitive_uuid],
                ))
        
        elif source == "multiple": 
            
            self.probe_points.append(
                sim_sdk.ProbePointsResultControl(
                    name= name,
                    write_control=sim_sdk.TimeStepWriteControl(write_interval=10),
                    geometry_primitive_uuids= self.geometry_primitive_uuid_list,
                ))
    

    def set_field_calculations(self):
        
        '''
        Not yet implemented... 
        
        define a porous media on a body/part or a geometry primitive
        
        Parameters
        ----------
                   
        Returns
        -------
        None. 
            
        '''            
        pass
    
    def set_result_control_items(self):
        
        '''
        store all of the pre-defined result contorl settings in one variable to pass
        that variable in a clean way to the simulation spec 
        
        just an organizational function
        
        Parameters
        ----------
        None. 
        
        Returns
        -------
        None. 
            
        '''              
        self.result_control = sim_sdk.FluidResultControls(
            surface_data = self.surface_data, 
            probe_points = self.probe_points, 
            field_calculations = self.field_calculations)
        
    def _set_contact_detection(self, method = "AUTO"):
        
        '''
        define the contact definition type. Default is AUTO and there is no other 
        choices.. 
        
        in CHTv2 the requirments to define contacts is not required. Interfaces between 
        bodies are handled by the meshing algorithim 
        
        contacts are only needed to be defined manually if one would like to define 
        an thermal resistance layer in an interface 
        
        Parameters
        ----------
        None. 
        
        Returns
        -------
        None. 
            
        '''              
        self.contact_detection = method
        
        return self.contact_detection
    
    """Functions to define the simulation spec"""
    
    def set_simulation_spec(self, simulation_name):
        
        '''
        Add everything that was defined before to the simulation spec 
        
        Parameters
        ----------
        None. 
        
        Returns
        -------
        None. 
            
        '''      
        # create the simulation spec first to pass as a reference to a mesh operation for physics based meshing
        
        self.simulation_spec = None 
        
        model = CoupledConjugateHeatTransfer(
            is_compressible= self.compressible,
            enable_radiation=self.radiation, 
            turbulence_model= self.turbulence_model,
            model= self.fluid_model_gravity,
            
            # initial_conditions = self.initial_conditions,
            initial_conditions = self._set_initial_conditions(),
            
            materials= sim_sdk.CoupledConjugateHeatTransferMaterials(),
            
            # numerics= self.fluid_numerics,
            numerics = self._set_simulation_numerics(),
            
            boundary_conditions= self.boundary_conditions
            ,
            advanced_concepts= self.advanced_concepts
            ,
            simulation_control= self.simulation_control
            ,
            result_control= self.result_control 
            ,
            contact_handling_mode= self._set_contact_detection(method = "AUTO")
        )
        
        self.simulation_spec = sim_sdk.SimulationSpec(name= simulation_name, 
                                                      geometry_id= self.geometry_id,
                                                      model=model)
        self.simulation_id = self.simulation_api.create_simulation(self.project_id, 
                                                                   self.simulation_spec).simulation_id
        
        print(f"simulation_id: {self.simulation_id}")
        
    def reset_simulation_spec_components(self):
        
        '''
        Reset the spec to an intital state to avoid name conflict with the setup of simulations of different CAD
        
        Parameters
        ----------
        None. 
        
        Returns
        -------
        None. 
        
            
        '''          
        # To avoid duplicate assignemnts between multiple simulations 
        self.fluid_material      = []
        self.solid_material      = []
        self.boundary_conditions = []
        self.velocity_inlet      = []
        self.velocity_outlet     = []
        self.pressure_inlet      = []
        self.pressure_outlet     = []
        self.wall_bc             = []
        self.natural_convection  = []
        self.power_sources       = []
        self.porous_media        = []
        self.momentum_sources    = []
        self.thermal_resistance_networks = []    
        self.surface_data        = []
        self.probe_points        = []
        self.field_calculations  = []
        self.mesh_refinement     = []        
        
    """Functions to define mesh settings"""
    
    def set_mesh_layer_settings(self, num_of_layers = 3, total_rel_thickness = 0.4, growth_rate = 1.5):
        
        '''
        define the mesh boundary layer settings         
        Parameters
        ----------
        
        num_of_layers: int
            
            number of boundary layers to be added 
            
        total_rel_thickness: float 
        
            Specifies the overall thickness of all layers combined, relative to the surface mesh element
            adjacent to the boundary layers. It is not recommended to set this value below 0.25.
        
        growth_rate: float
        
            controls the relative sizing between elements of adjacent layers. 1.5 means a 50% difference between layers.
            
        Returns
        -------
        None. 
        
        '''           
        self.automatic_layer_settings = sim_sdk.AutomaticLayerOn(type="AUTOMATIC_LAYER_ON",
                                                                  		number_of_layers= num_of_layers,
                                                                  		total_relative_thickness= total_rel_thickness,
                                                                  		layer_type= sim_sdk.FractionalHeight2(
                                                                  			type="FRACTIONAL_HEIGHT_2",
                                                                  			growth_rate= growth_rate,)
                                                                          )
        
    def set_advanced_mesh_settings(self, small_feature_tolerance = 5E-5, gap_ref_factor = 0.05, gradation_rate = 1.22): 
        
        '''
        define some advanced mesh settings layer settings         
        Parameters
        ----------
        
        small_feature_tolerance: float
            
            This can be used to ignore small surfaces during the meshing process.
            It essentially merges surfaces to avoid unnecessary levels of refinement around small features. 
            
        gap_ref_factor: float 
        
            allows to better capture small gaps in the model
            For values > 1, the value (integer) is also the number of elements across the gap.
            For values < 1, the gap refinement factor is the inverse of the aspect ratio of the elements in the gap.   
            
        gradation_rate: float
        
            ratio between the size of two adjacent cells in the computational domain.
            This setting controls how quickly the cells transition from small to large.   
            
        Returns
        -------
        None. 
        
        '''             
        self.advanced_mesh_settings = sim_sdk.AdvancedSimmetrixFluidSettings(
                                    small_feature_tolerance= sim_sdk.DimensionalLength(value= small_feature_tolerance,unit="m",),
                                    gap_elements = gap_ref_factor, global_gradation_rate = gradation_rate)
    
    
    def complete_mesh_settings(self, mesh_name ,fineness = 5, physics_based_meshing = True):
        
        '''
        combine the mesh settings in one place and define the global mesh fineness  
        
        Parameters
        ----------
        
        mesh_name: str
            
            the name of the mesh being submitted 
            
        fineness: int  
        
            The global fineness level of the mesh. define values between 1 and 10 
            1  being very coarse
            10 being very fine 
            
            
        physics_based_meshing: boolean 
        
            when enabled, builds the mesh taking into account the information entered as a part of the simulation setup.
            This important information ranges from the material properties, boundary conditions,
            additional source terms (e.g., momentum and power sources), etc. 
            Essentially, the physics involved in the fluid simulation is given priority while sizing the mesh elements.
            
        Returns
        -------
        None. 
        
        '''      
        # Start of mesh operation
        
        self.mesh_operation = self.mesh_operation_api.create_mesh_operation(
            self.project_id,
            sim_sdk.MeshOperation(
                name= mesh_name,
                geometry_id= self.geometry_id,
                model= sim_sdk.SimmetrixMeshingFluid(
                    physics_based_meshing= physics_based_meshing, hex_core = True, 
                    sizing=AutomaticMeshSizingSimmetrix(type="AUTOMATIC_V9",fineness= fineness),
                    refinements = self.mesh_refinement,
                    automatic_layer_settings= self.automatic_layer_settings,
                    advanced_simmetrix_settings = self.advanced_mesh_settings
                 ),
                ),
            )
        
        self.mesh_operation_id = self.mesh_operation.mesh_operation_id
        self.mesh_operation_api.update_mesh_operation( self.project_id, self.mesh_operation_id, self.mesh_operation)

    def set_local_element_size_refinement(self, max_element_size ,name = '',  key_list = []):
        
        '''
        add a local element size mesh refinement
        
        Parameters
        ----------
        
        max_element_size: int or float 
            
            maximum edge length for the mesh elements.
            decide on a value based on the size of the edges you are trying to refine.
            
        name: str  
        
            name of the local mesh refinement 
                
        key_list : list
        
            the desired faces/parts for later topological assignment are stored 
            in a dict. each face/part has a specific user defined key, that is 
            passed to this function. In the end, this key would point to the 
            desired entitiy id of a face/part.
            
            in this case, you are supplying a list of the keys that points 
            to the parts that are expected to be assigned a certain material 
              
        Returns
        -------
        None. 
        
        '''    
        #for later: add a method to add multiple entity selections automatically using get_entities
        
        faces_to_assign = []
        for key in key_list: 
            faces_to_assign.append(self.single_entity[key])
        
        self.mesh_refinement.append(   
                        sim_sdk.SimmetrixLocalSizingRefinement(
            			type="SIMMETRIX_LOCAL_SIZING_V10",
            			name= name,
            			max_element_size= sim_sdk.DimensionalLength(
            				value= max_element_size,
            				unit="m",
            			),
            			topological_reference= sim_sdk.TopologicalReference(
            				entities= faces_to_assign,
            				sets=[],
            			),
            		))
        
    def estimate_mesh_operation(self):
        
        '''
        Estimate how much time and computing resources a mesh operation is going to take 
        
        Parameters
        ----------
        None. 
        
        Returns
        -------
        None. 
        
        '''            
        # Estimate Mesh operation
        try:
            mesh_estimation = self.mesh_operation_api.estimate_mesh_operation(self.project_id, self.mesh_operation_id)
            # print(f"Mesh operation estimation: {mesh_estimation}")
            print("*"*10)
            print("Mesh operation estimation:")
            print("Number of cells: {lower} - {upper}. Expected is {avg}\n".format(
                                                               lower =  mesh_estimation.cell_count.interval_min, 
                                                               upper =  mesh_estimation.cell_count.interval_max,
                                                               avg   =  mesh_estimation.cell_count.value))
            
            print("CPU consumption: {lower} - {upper}. Expected is {avg}\n".format(
                                                               lower =  mesh_estimation.compute_resource.interval_min, 
                                                               upper =  mesh_estimation.compute_resource.interval_max,
                                                               avg   =  mesh_estimation.compute_resource.value))
            
            print("Duration: {lower} - {upper}. Expected is {avg}\n".format(
                                                               lower =  mesh_estimation.duration.interval_min.replace('PT',''), 
                                                               upper =  mesh_estimation.duration.interval_max.replace('PT', ''),
                                                               avg   =  mesh_estimation.duration.value.replace('PT','')))          
            print("*"*10)
            
            if mesh_estimation.compute_resource is not None and mesh_estimation.compute_resource.value > 150.0:
                raise Exception("Too expensive", mesh_estimation)
        
            if mesh_estimation.duration is not None:
                self.mesh_max_runtime = isodate.parse_duration(mesh_estimation.duration.interval_max).total_seconds()
                self.mesh_max_runtime = max(3600, self.mesh_max_runtime * 2)
            else:
                self.mesh_max_runtime = 36000
                print(f"Mesh operation estimated duration not available, assuming max runtime of {self.mesh_max_runtime} seconds")
        except ApiException as ae:
            if ae.status == 422:
                self.mesh_max_runtime = 36000
                print(f"Mesh operation estimation not available, assuming max runtime of {self.mesh_max_runtime} seconds")
            else:
                raise ae
                
    def check_simulation_and_mesh_settings(self):
        
        '''
        do a sanity check on whether there are new warnings/errors in either 
        the simulation setup or the mesh
        
        Parameters
        ----------
        None. 
        
        Returns
        -------
        None. 
        
        '''            
        mesh_check = self.mesh_operation_api.check_mesh_operation_setup(self.project_id, self.mesh_operation_id, simulation_id= self.simulation_id)
        warnings = [entry for entry in mesh_check.entries if entry.severity == "WARNING"]
        if warnings: 
            # print(f"Meshing check warnings: {warnings}")
            print("Meshing check warning: {}\n".format(warnings[0].message))
        if len(warnings) > 1 :
            print("*"*10)
            # print(warnings[1].message)
            print("\nSimulation setup check warnings: {}".format(warnings[1].message))
        errors = [entry for entry in mesh_check.entries if entry.severity == "ERROR"]
        if errors:
            raise Exception("Simulation check failed", mesh_check)

    def start_meshing_operation(self, run_state = False):
        
        '''
        submit the mesh operation to start 
        
        Parameters
        ----------
        run_state: boolean 
            
            a boolean to decide on whether the mesh is put on queue or not 
            
            (if you are still reading this line, then please keep run_state False, 
             I would still need to test this thorougly )
        
        Returns
        -------
        None. 
        
        '''    
        #for later: Figure out how to run the mesh operations in parallel without having to wait for each one to finish
        # The problem is related with retreiving the mesh id. It is possible that it is only attainable after the mesh operation is completed - investiagte...
        
        
        if run_state :
            self.mesh_operation_api.start_mesh_operation(self.project_id, self.mesh_operation_id, simulation_id= self.simulation_id)
        
        else: 
            
            self.mesh_operation_api.start_mesh_operation(self.project_id, self.mesh_operation_id, simulation_id= self.simulation_id)            
            # Wait until the meshing operation is complete
            self.mesh_operation = self.mesh_operation_api.get_mesh_operation(self.project_id, self.mesh_operation_id)
            mesh_operation_start = time.time()
            while self.mesh_operation.status not in ("FINISHED", "CANCELED", "FAILED"):
                if time.time() > mesh_operation_start + self.mesh_max_runtime:
                    raise TimeoutError()
                time.sleep(30)
                self.mesh_operation = self.mesh_operation_api.get_mesh_operation(self.project_id, self.mesh_operation_id)
                print(f"Meshing run status: {self.mesh_operation.status} - {self.mesh_operation.progress}")
            
            self.mesh_operation = self.mesh_operation_api.get_mesh_operation(self.project_id, self.mesh_operation_id)
            # print(f"final mesh_operation: {self.mesh_operation}")
            
            # Get the simulation spec and update it with mesh_id from the previous mesh operation
            self.simulation_spec = self.simulation_api.get_simulation(self.project_id, self.simulation_id)
            self.simulation_spec.mesh_id = self.mesh_operation.mesh_id
            self.simulation_api.update_simulation(self.project_id, self.simulation_id, self.simulation_spec)

    """Functions to estimate the simulation resources, create and run the simulation"""
    
    def estimate_simulation(self, maximum_cpu_consumption_limit = 400):
        
        '''
        Estimate how much time and computing resources a mesh operation is going to take 
        
        Parameters
        ----------
        maximum_cpu_consumption_limit: int 
        
            define a maximum limit that you wouldn't want the simulation to consume more than 
            if the simulation estimation finds out that there a possiblity to exceed the 
            maximum defined cpu consumption limit, then an exception is raised 
            
        Returns
        -------
        None. 
        
        '''               
        try:
            estimation = self.simulation_api.estimate_simulation_setup(self.project_id, self.simulation_id)
            # print(f"Simulation estimation: {estimation}")
            
            print("*"*10)
            print("CPU consumption: {lower} - {upper}. Expected is {avg}\n".format(
                                                               lower =  estimation.compute_resource.interval_min, 
                                                               upper =  estimation.compute_resource.interval_max,
                                                               avg   =  estimation.compute_resource.value))
            
            print("Duration: {lower} - {upper}. Expected is {avg}\n".format(
                                                               lower =  estimation.duration.interval_min.replace('PT',''), 
                                                               upper =  estimation.duration.interval_max.replace('PT', ''),
                                                               avg   =  estimation.duration.value.replace('PT','')))          
            print("*"*10)
            
            if estimation.compute_resource is not None and estimation.compute_resource.value > maximum_cpu_consumption_limit:
                raise Exception("Too expensive", estimation)
        
            if estimation.duration is not None:
                self.sim_max_run_time = isodate.parse_duration(estimation.duration.interval_max).total_seconds()
                self.sim_max_run_time = max(3600, self.sim_max_run_time * 2)
            else:
                self.sim_max_run_time = 36000
                print(f"Simulation estimated duration not available, assuming max runtime of {self.sim_max_run_time} seconds")
        except ApiException as ae:
            if ae.status == 422:
                self.sim_max_run_time = 36000
                print(f"Simulation estimation not available, assuming max runtime of {self.sim_max_run_time} seconds")
            else:
                raise ae
        
    def create_simulation(self, sim_name ):
        
        '''
        create a simulation and prepare it for submitting to run 
        
        Parameters
        ----------
        
        sim_name: str  
        
            name of the simulation run
        
        Returns
        -------
        None. 
        
        '''             
        self.simulation_run = sim_sdk.SimulationRun(name= sim_name)
        self.simulation_run = self.simulation_run_api.create_simulation_run(self.project_id, self.simulation_id, self.simulation_run)
        self.run_id = self.simulation_run.run_id
        print(f"runId: {self.run_id}")
        
        # Read simulation run and update with the deserialized model
        self.simulation_run = self.simulation_run_api.get_simulation_run(self.project_id, self.simulation_id, self.run_id)
        self.simulation_run_api.update_simulation_run(self.project_id, self.simulation_id, self.run_id, self.simulation_run)

    def start_simulation_run(self, wait_for_results = True): 
        
        '''
        submit the simulation to solve
        
        Parameters
        ----------
        
        wait_for_results: boolean  
        
            if True - then the simulation run is put on queue and the process 
            waits until the simulation results are obtained before moving on 
            
            this can be helpful if you want to download the results or process them 
            automatically in the same script, which means that you would need 
            to get the results first 
        
        Returns
        -------
        None. 
        
        '''            
        #Run simulation and don't wait until results are finished 
        if not wait_for_results : 
            # Start simulation run 
            self.simulation_run_api.start_simulation_run(self.project_id, self.simulation_id, self.run_id)
            self.simulation_run = self.simulation_run_api.get_simulation_run(self.project_id, self.simulation_id, self.run_id)
            
        #Run simulation and wait until results are finished 
        else: 
            self.simulation_run_api.start_simulation_run(self.project_id, self.simulation_id, self.run_id)
            self.simulation_run = self.simulation_run_api.get_simulation_run(self.project_id, self.simulation_id, self.run_id)
            simulation_run_start = time.time()
            while self.simulation_run.status not in ("FINISHED", "CANCELED", "FAILED"):
                if time.time() > simulation_run_start + self.sim_max_run_time:
                    raise TimeoutError()
                time.sleep(30)
                self.simulation_run = self.simulation_run_api.get_simulation_run(self.project_id, self.simulation_id, self.run_id)
                print(f"Simulation run status: {self.simulation_run.status} - {self.simulation_run.progress}")
            
    """Functions to download the results"""
    
    def get_simulation_results(self): 
        
        '''
        get information on the simulation results that can be used for further processing 
        
        Parameters
        ----------
        None. 
        
        Returns
        -------
        None. 
        
        '''
        #for later: add functions that allows retreiving the project, simulation and run id of any project
        self.simulation_results = self.simulation_run_api.get_simulation_run_results(self.project_id, self.simulation_id, self.run_id)
        print(f"results: {self.simulation_results}")

    def get_probe_point_results(self, name, field = 'T', dir_name = "probe_point_results/"):
        
        '''
        get the results from a probe point result control item 
        you can specifiy which field in particular you would like to get, and then it saves
        those results to a directory you specifiy in a .csv format 
        
        Parameters
        ----------
        
        name : str 
        
            name of the result control item defined in the simulation setup process
            
        field : str  
            
            the name of the field result of interest ex : 'Ux', 'Uy', 'Uz', 'p', 'k', 'T'
        
        dir_name: str 
        
            a string that has the name of the directory where the results are to be saved 
            as of now the directory specified should be in the current working directory 
            
        Returns
        -------
        None. 
        
        '''        
        #for later: add code that will run a validation if the probe points don't all have a field 
        # catch the error and move on. Do not terminate the script
        # and r.name == name
        self.probe_point_plot_info = [r for r in self.simulation_results._embedded if (r.category == "PROBE_POINT_PLOT" and r.name == name and r.quantity == field)][0]
        print(self.probe_point_plot_info)
        
        self.probe_point_plot_data_response = self.api_client.rest_client.GET(
            url= self.probe_point_plot_info.download.url, headers={self.api_key_header: self.api_key}, _preload_content=False
        )
        probe_point_plot_data_csv = self.probe_point_plot_data_response.data.decode("utf-8")
        print(f"Probe point plot data as CSV: {probe_point_plot_data_csv}")
        
        # Write probe points to CSV file
        probe_results_path = pathlib.Path(dir_name)
        try:
            #check if the directory already exists if not create a new one and store in it 
            probe_results_path.mkdir(parents = True, exist_ok = False)
            print("D")
            write_to = probe_results_path / "{}.csv".format(name)
            with open(write_to, "w") as file:
                file.write(probe_point_plot_data_csv)
        except: 
            #write to the already existing directory 
            probe_results_path.mkdir(parents = True, exist_ok = True)
            write_to = probe_results_path / "{}.csv".format(name)
            with open(write_to, "w") as file:
                file.write(probe_point_plot_data_csv)
    
    def get_surface_data_results(self, name, data_type = "average" , field = 'T', dir_name = "surface_data_results"): 
        
        '''
        get the results from a surface data result control item 
        you can specifiy which field in particular you would like to get, and then it saves
        those results to a directory you specifiy in a .csv format 
        
        Parameters
        ----------
        
        name : str 
        
            name of the result control item defined in the simulation setup process
            
        data_type : str 
        
            either an 'average' or 'volume' surface data control item 
                        
        field : str  
            
            the name of the field result of interest ex : 'Ux', 'Uy', 'Uz', 'p', 'k', 'T'
        
        dir_name: str 
        
            a string that has the name of the directory where the results are to be saved 
            as of now the directory specified should be in the current working directory 
            
        Returns
        -------
        None. 
        
        '''                    
                
        if data_type == 'average': 
            area_average_result = [r for r in self.simulation_results._embedded if (r.category == "AREA_AVERAGE" and r.name == name and r.quantity == field)][0]
        
        else: 
            area_average_result = [r for r in self.simulation_results._embedded if (r.category == "AREA_INTEGRAL" and r.name == name and r.quantity == field)][0]
   
        area_average_result_response = self.api_client.rest_client.GET(
            url= area_average_result.download.url, headers={self.api_key_header: self.api_key}, _preload_content=False
        )
        area_average_results_csv = area_average_result_response.data.decode("utf-8")
        print(f"Area average result as CSV: {area_average_results_csv}")

        # Write area averages to CSV file
        area_average_result_path = pathlib.Path(dir_name)
        try:
            #check if the directory already exists if not create a new one and store in it 
            area_average_result_path.mkdir(parents = True, exist_ok = False)
            print("D")
            write_to = area_average_result_path / "{n}_{f}.csv".format(n = name, f = field)
            with open(write_to, "w") as file:
                file.write(area_average_results_csv)
        except: 
            #write to the already existing directory 
            area_average_result_path.mkdir(parents = True, exist_ok = True)
            write_to = area_average_result_path / "{n}_{f}.csv".format(n = name, f = field)
            with open(write_to, "w") as file:
                file.write(area_average_results_csv)
                
    def get_simulation_case_files(self): 
        
        '''
        download the simulation case files 
        
        Parameters
        ----------
        None. 
        
        Returns
        -------
        None. 
        
        '''                
        self.solution_info = [r for r in self.simulation_results._embedded if r.category == "SOLUTION"][0]
        solution_response = self.api_client.rest_client.GET(
            url=self.solution_info.download.url, headers={self.api_key_header: self.api_key}, _preload_content=False)
        with open("case_file_solution.zip", "wb") as file:
               file.write(solution_response.data)
        zip = zipfile.ZipFile("case_file_solution.zip")
        print(f"Averaged solution ZIP file content: {zip.namelist()}")
    
    def get_simulation_report(self, part_id = "region1"):
        
        '''
        generate a simulation report  
        
        Parameters
        ----------
        part_id: str 
        
            the name of the region the results are stored in. 
            
            for example, if you go to the post-processor or in paraview 
            you would find the name of regions where results are stored 
            
            (anyway, this should be enhanced, it is not robust)
        
        Returns
        -------
        None. 
        
        '''          
        
        # Generating simulation run report
        camera_settings = sim_sdk.UserInputCameraSettings(
            projection_type= sim_sdk.ProjectionType.ORTHOGONAL,
            up= sim_sdk.Vector3D(0.5, 0.3, 0.2),
            eye= sim_sdk.Vector3D(0.0, 5.0, 10.0),
            center= sim_sdk.Vector3D(10.0, 12.0, 1.0),
            front_plane_frustum_height=0.5,
        )
        # "Temperature", component="X", data_type="CELL"
        model_settings = sim_sdk.ModelSettings(
            parts=[Part(part_identifier= part_id, solid_color=Color(0.8, 0.2, 0.4))],
            scalar_field=sim_sdk.ScalarField(field_name= "Velocity", component="X", data_type="CELL"),
        )
        output_settings = sim_sdk.ScreenshotOutputSettings(name="Output 1", format="PNG", resolution=ResolutionInfo(800, 800),
                                                   frame_index=0)
        report_properties = sim_sdk.ScreenshotReportProperties(
            model_settings=model_settings,
            filters=None,
            camera_settings=camera_settings,
            output_settings=output_settings,
        )
        report_request = sim_sdk.ReportRequest(name="Report 1", description="Simulation report", result_ids=[self.solution_info.result_id],
                                       report_properties=report_properties)
        
        create_report_response = self.reports_api.create_report(self.project_id, report_request)
        report_id = create_report_response.report_id
        
        # Start report job
        print(f"Starting report with ID {report_id}")
        report_job = self.reports_api.start_report_job(self.project_id, report_id)
        
        report = self.reports_api.get_report(self.project_id, report_id)
        
        while report.status not in ("FINISHED", "CANCELED", "FAILED"):
            time.sleep(30)
            report = self.reports_api.get_report(self.project_id, report_id)
        
        print(f"Report finished with status {report.status}")
        
        if report.status == "FINISHED":
            # Download the report
            print("Downloading report result")
            report_response = self.api_client.rest_client.GET(
                url=report.download.url,
                headers={self.api_key_header: self.api_key},
                _preload_content=False,
            )
        
            file_name = f"report.{report.download.format}"
            with open(file_name, "wb") as file:
                file.write(report_response.data)
                print(f"Finished downloading report with name {file_name}")
        elif report.status == "FAILED":
            raise Exception("Report generation failed", report.failure_reason)


    """Functions to find existing projects, simulations, runs, geometries"""    

    def find_simulation(self, name):
        
        '''
        Take a Simulation Name, return a simulation
    
        Parameters
        ----------
        name : string
            The exact name of the simulation, best copied from the SimScale 
            UI.
        project_id : object
            the ID of the project that you are searching, Simulations are 
            child objects of projects.
        simulation_api : object
            An API object that can be used for querying and creating 
            SimScale simulations..
    
        Raises
        ------
        Exception
            Raise exception if the name matches no simulation in the 
            project.
    
        Returns
        -------
        found : object
            A simulation object that was matched by the provided name.
    
        '''
        project_id = self.project_id
        simulation_api = self.simulation_api
    
        simulations = simulation_api.get_simulations(project_id).to_dict()['embedded']
        found = None
        for simulation in simulations:
            if simulation['name'] == name:
                found = simulation
                print('Simulation found: \n' + str(found['name']))
                break
        if found is None:
            raise Exception('could not find simulation with id: ' + name)
        self.simulation = found
        self.simulation_id = found["simulation_id"]
    
    def find_run(self, name):
        '''
        Take, name, parent simulation and parent project, return run.
    
        Parameters
        ----------
        name : string
            The exact name of the simulation run, best copied from the 
            SimScale UI.
        project_id : object
        simulation_id : object
        simulation_run_api : object
            An API object that can be used for querying and creating and 
            downloading SimScale simulation runs.
    
        Raises
        ------
        Exception
            Raise exception if the name matches no run in the parent 
            simulation.
    
        Returns
        -------
        found : TYPE
            DESCRIPTION.
    
        '''
        project_id = self.project_id
        simulation_id = self.simulation_id
        simulation_run_api = self.simulation_run_api
    
        runs = simulation_run_api.get_simulation_runs(
            project_id, simulation_id).to_dict()['embedded']
    
        found = None
        for run in runs:
            if run['name'] == name:
                found = run
                print('Run found: \n' + str(found['name']))
                break
        if found is None:
            raise Exception('could not find simulation with id: ' + name)
        self.run = found
        self.run_id = found["run_id"]
        self.run_id = found["run_id"]
    
    # def find_geometry(self, name):
    #     '''
    #     Take a Simulation Name, return a simulation
    
    #     Parameters
    #     ----------
    #     name : string
    #         The exact name of the simulation, best copied from the SimScale
    #         UI.
    #     project_id : object
    #         the ID of the project that you are searching, Simulations are
    #         child objects of projects.
    #     simulation_api : object
    #         An API object that can be used for querying and creating
    #         SimScale simulations..
    
    #     Raises
    #     ------
    #     Exception
    #         Raise exception if the name matches no simulation in the
    #         project.
    
    #     Returns
    #     -------
    #     found : object
    #         A simulation object that was matched by the provided name.
    
    #     '''
    #     project_id = self.project_id
    #     geometry_api = self.geometry_api
    
    #     geometries = geometry_api.get_geometries(project_id).to_dict()['embedded']
    #     found = None
    #     for geometry in geometries:
    #         if geometry['name'] == name:
    #             found = geometry
    #             print('Geometry found: \n' + str(found['name']))
    #             break
    #     if found is None:
    #         raise Exception('could not find geometry with id: ' + name)
    #     self.geometry = found
    #     self.geometry_id = found["geometry_id"]
    
    # def find_project(self, name):
    #     '''
    #     Take a project Name, return a project
    
    #     Parameters
    #     ----------
    #     name : string
    #         The exact name of the project, best copied from the SimScale
    #         UI.
    #     project_api : object
    #         An API object that can be used for querying and creating
    #         SimScale projects.
    
    #     Raises
    #     ------
    #     Exception
    #         Raise exception if the name matches no project in the clients*
    #         account.
    
    #         *client is created using the users API key, see create_client.
    #     Returns
    #     -------
    #     found : object
    #         A simulation object that was matched by the provided name.
    
    #     '''
    #     project_api = self.project_api
    
    #     projects = project_api.get_projects(limit=1000).to_dict()['embedded']
    #     found = None
    #     for project in projects:
    #         if project['name'] == name:
    #             found = project
    #             print('Project found: \n' + str(found['name']))
    #             break
    #     if found is None:
    #         raise Exception('could not find project with name: ' + name)
    
    #     self.project_id = found['project_id']
    #     self.project = name
    

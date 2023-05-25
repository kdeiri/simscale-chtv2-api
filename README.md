# SimScale_CHTv2
The main goal of this code is to demonstrate how to run conjugate heat transfer (CHTv2) simulations on SimScale via the API. 

The repository already contains a working example that can be used for testing. The example is of a cold cooling plate for a battery system. 

This code should serve the purpose of showing: 
1. How to setup the API connection. 
2. How to create a project and upload geometries 
3. How to retrieve entity IDs of parts and surfaces, so that they can be used for assignmnet of boundary conditions, advanced concepts, etc.
4. How to define a material directly from the default SimScale material library via the API 
5. How to apply boundary conditions via the API (almost all boundary conditions are coded in the utilities file for reference) 
6. How to define advanced concpets (power sources, momentum sources, porous media: would be added in the future)
7. How to define result control items (area averages, area integrals, probe points) 
8. How to define geometry primtives( points, cartesian box, cylinders) 
9. How to create a mesh and adding mesh refinements 
10. How to check for the computing resources of a mesh and a simulation 
11. Run the simulation 
12. Download results 

To use the code you would need to have a SimScale API key - you can receive one after you contact the SimScale team. 



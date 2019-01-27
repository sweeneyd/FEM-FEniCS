// Gmsh project created on Sat Jan 26 17:45:16 2019
Merge "cylinder_flow.step";
//+
Physical Surface(1) = {2};
//+
Physical Surface(2) = {3};
//+
Physical Surface(3) = {1};
//+
Physical Volume(4) = {1};
Mesh.CharacteristicLengthMax = 0.35; //0.2
Mesh 3;
OptimizeMesh "GMSH";
//RefineMesh;
//OptimizeMesh "GMSH";
//RefineMesh;
//OptimizeMesh "GMSH";
//RefineMesh;
//OptimizeMesh "GMSH";
Coherence Mesh;
Save 'cylinder_flow.msh';

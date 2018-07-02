# CodendOOP
Simulation of trawl cod-end shapes during the towing using an object-oriented C# version of 3D finite element model 
based on triangular elements developed by D.Priour (see [link](http://www.ifremer.fr/web-com/dpriour/web/papers/2005_lisbon_FEM.pdf)
and thesis). Suggested to be used in combination with [Axis-symmetric model](https://github.com/mihsamusev/AxiCodend) 
for generation of initial shape and better robustness of the used calculation algorighm. Input happens throught the "input.txt" file, and should follow the format:

*Block of paths (essential input) contains the paths where the results are stored, 
the input is by default input.txt and has to be in the same directory as .exe file. 
In addition if Axis-symmetric model is used, its .exe and input.txt paths need to be provided.*
```
PATHS
AxiExe		C:\Users\MyFavoriteFolder\AxiCodend\AxiCodend\bin\Debug\AxiCodend.exe
InputAxi	C:\Users\MyFavoriteFolderAxiCodend\AxiCodend\bin\Debug\input.txt
OutputAxiShapes C:\Users\MyFavoriteFolder\precalcShapes.txt
Output3dShapes	C:\Users\MyFavoriteFolder\demoShape.txt
Output3dResults	C:\Users\MyFavoriteFolder\demoResults.txt
```
*Block of cod-end assemply (essential input). Here the stcructure of the cod-end is specified.* 
```
CODEND ASSEMBLY
EntranceRadius 0.3
```
*Describing the panels, it is important to have the same amout of panel rows as* `PanelCount`. 
*User can choose between* `DiamondPanel` *and* `SquarePanel`. *For* `DiamondPaned` *orientation can be* `0` *or* `90`.
```
PANELS
PanelCount 2
ID	Type		MaterialID	MeshesAlong	MeshesAcross	Orientation
1	DiamondPanel 	1 		50 		25		0
1	DiamondPanel 	1 		50 		25		0
```
*There can be many panel materials, but the program will only use those referenced in* `MaterialID`. 
```
PANEL MATERIALS
PanelMaterialCount 2
ID	Type	Density	MeshSide	IsDoubleTwine	TwineThickness	InitialOpeningAngle	EA	EI	OpenningStifness
1	Diamond	1025	0.1		0		0.0014		45			1000	0.002	0.2	
2	Square	1025	0.1		0		0.0014		45			1000	0.002	0.2400
```
*Selvedges can be included by putting the flag of* `1` *next to* `IncludeSelvedges`. *Their amount depends on the panel count.*
```
SELVEDGES
IncludeSelvedges 0
Length	Density	Diameter	EA	EI
4.5	1025	0.012		5000	0
```
*Round straps can be included by putting the flag of* `1` *next to* `IncludeRoundStraps`.
```
ROUND STRAPS
IncludeRoundStraps 1
RoundStrapCount 1
ID	Length	Position	Density	 Diameter	EA	EI
1	2.5	0.3		1025	0.012		5000	0
```
*Block of Towing (essential input). Netting drag can be included using corresponding flag.*
```
TOWING
IncludeNettingDrag 0
CX 0
CY 0.8
CZ 0
```
*Block of catches (essential input). Can be applied by two different* `ApplyMethod` *that are:*
`ByBlockedMeshes` *specifying the amount of meshes blocked along the cod-end (relevant for cod-ends made of diamond mesh panels only)*
`ByBlockingRatio` *specifying the fraction of the length from the rear end 0-1 that is blocked by catch.
Make sure that Count matches the amount of elements in the following catches column.*
```
CATCH
CatchCount 5
ApplyMethod ByBlockedMeshes
20
25
30
35
40
```
*Block of panel meshing (optional input). Specifies meshing settings for each panel. The* `MeshingMethod` *can be:*
`ByElement` *specifying how many meshes should be inside an element along or across.*
`ByMesh` *specifying how many elements there are along or across the panel.*
*By default the settings are as following:*
```
MESHING
MeshingMethod ByElement
ElemAcrossPanel 10
ElemAlongPanel 10
```
*If another method is chosen, then the parameters are set using keywords* `MeshPerElemAlong` *and* `MeshPerElemAcross`


*Block of Solver settings (optional), if nothing is input the solver uses following default settings:*
```
SOLVER PARAMETERS
IterMax             5000
ResidualTol         1e-3
DisplacementTol     1e-4
ResidualMax         1e20
StiffnessTol        1
DiagStiffness       1
ReduceStiffnessBy   0.1
IncreaseStiffnessBy 2
ShowLineSearchSteps 0
LineSearchMax       6
ReductionCondition  1e-4
MinCatchBlock       5
MinTowingSpeed      1.0
IncludeLineSearch   0
InitialShape        cylinder
```

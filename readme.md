# As Rigid As Possible Surface Modelling
#### _In Python_
#### Tana Tanoi
#### CGRA409

### Command Line Arguments

`python3 arap.py $1 $2 $3 $4`

##### Argument 1
A `.off` file to be deformed

##### Argument 2
A `.sel` file that selects the handles and fixed points of the mesh.
Each line corresponds to a vertex in the `.off` file, where a **2** indicates a handle, a **1** indicates a flexible/deformable point, and **0** is a fixed/constrained point.

##### Argument 3
A `.def` file that represents the deformation to apply to the handles (from argument 2). This is a 4 by 4 matrix, where each number is separated by white space.

##### Argument 4
An integer specifying the max number of iterations to run. _It may run fewer than this if the energy difference between iterations is less than 0.01._

Performs acceptably on smaller meshes (sub 1k faces), but has trouble with larger meshes.
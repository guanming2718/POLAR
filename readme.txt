
1 Run the code by typing
python3 test.py -in <inputfile> -out <outputdirectory> -v [1D or 2D]

2 The input parameters
[Configuration]
# Configuration parameters
nx             :  number of points in x direction
ny             :  number of points in y direction
Lx             :  size of the simulation box in x dimension
BC_density     :  boundary condition for the density = 'pbc','no flux','zero'
BC_polar       :  boundary condition for the polar field = 'pbc','free-slip','zero'
BC_velocity    :  boundary condition for the velocity field  = 'pbc','free-slip','zero'
init-config    :  initial configuration = (confluent,cluster,one-side-wound or two-side-wound)
[Model]
#parameters for the model
Ac             :  free energy of the density 0.5*Ac*[rho(rho-1)]^2
Kc             :  elasticity of the density Kc*[grad(rho)]^2
#The free energy of the polar terms
Ap             :  free energy of the polar field 0.5*Ap*p^2
Bp             :  free energy of the polar field 0.25*Bp*p^4
Kp             :  elastic term of polar free energy  Kp*(\nabla p)^2
w              :  contact inhibition coupling parameter     
gamma          :  the relaxatoin mobility
nu             :  tumbling parameter
zeta           :  the dipolar activity
alpha          :  self-propulsion activity
Ks             :  shear modulus
Kv             :  compressibility
xi             :  cell-substrate friction



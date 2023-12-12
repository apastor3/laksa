function NP  =  Load_Numerical_Parameters

% Numerical parameters
N_fe         = 2;         % Number of finite elements
N_Int        = 5;         % Number of point to carry out the integrals inside an element

Mode         = 5;         % 0 -> Test Global functions
                          % 1 -> Test the code with the catenary solution
                          % 3 -> Study Catenary solution 
                          % 4 -> Study periodic orbits 
                          % 5 -> Integrate the equations of motion
                          % 6 -> Bifurcation Diagram
% Integrator
RelTol       = 1e-10;
AbsTol       = 1e-10;
Tol          = 1e-8;       % Integrator Tolerance
dh           = 1e-5;

Newton_Tol   = 1e-8;
Max_It       = 100;


% Derived Numerical Parameters
N_nc         = N_Int*N_fe+1;     % Number of points for Newton-Cotes integration
L_fe         = 1/N_fe;           % Length of finite element
L_nc         = L_fe/N_Int;       % Length of Newton-Cotes segment
Np           = N_fe+1;           % Number of points
N_Dim        = 3*(N_fe-1)+Np;    % Dimension of the state vector
 
NP           = [N_fe N_Int Mode N_nc L_fe L_nc Np N_Dim Tol RelTol AbsTol dh Newton_Tol Max_It ]; 


end

mex -v '-I/usr/local/Cellar/eigen/3.3.4/include/eigen3' ...
     '-I/usr/local/Cellar/suite-sparse/5.1.0/include' ...
     SIMengine/Classes/stdafx.cpp SIMengine/Classes/mersenne.cpp ...
     SIMengine/Classes/SIMcore.cpp SIMengine/Classes/SIMparameters.cpp ...
     SIMengine/Classes/ChemotaxisMapSuiteS.cpp SIMengine/Classes/NecrosisMapSuiteS.cpp...
     SIMengine/Classes/Environment.cpp SIMengine/Classes/TumorCells.cpp ...
     SIMengine/Classes/Lymphocytes.cpp SIMengine/Classes/Macrophages.cpp ...
     SIMengine_interface_mex.cpp ...
 -lamd -lcamd -lccolamd -lcholmod -lcolamd -lmetis -lsuitesparseconfig -lcxsparse -lblas -llapack

 copyfile('stdafx.mexmaci64','SIMengine_interface_mex.mexmaci64')
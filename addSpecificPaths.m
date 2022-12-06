exampleRepo = '../CertifiablyRobustPerception/RotationSearch/';
cvxpath     = strcat(exampleRepo, '../../cvx');
mosekpath   = strcat(exampleRepo,'../../mosek');
utilspath   = strcat(exampleRepo,'../utils');
manoptpath  = strcat(exampleRepo,'../manopt');
sdpnalpath  = strcat(exampleRepo,'../../SDPNAL+v1.0');
stridepath  = strcat(exampleRepo,'../STRIDE');
sostools    = strcat(exampleRepo,'../../SOSTOOLS');

addpath(genpath(exampleRepo))
addpath(genpath(cvxpath))
addpath(genpath(mosekpath))
addpath(genpath(sostools))
addpath(genpath(utilspath))
addpath(genpath('./lib'))
addsostools
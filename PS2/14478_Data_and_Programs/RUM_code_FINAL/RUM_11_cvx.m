%% Code Description: Setup cvx
% This function installs and set paths for cvx.
% Make sure licence.dat file is put in cvx/cvx_licence folder.

%% Check version of MatLab
if mexext == 'mexw32'
    cvx_path = 'cvx/cvx_win32/';
elseif mexext == 'mexw64'
    cvx_path = 'cvx/cvx_win64/';
elseif mexext == 'mexglx'
    cvx_path = 'cvx/cvx_linux32/';
elseif mexext == 'mexa64'
    cvx_path = 'cvx/cvx_linux64/';
elseif mexext == 'mexmaci'
    cvx_path = 'cvx/cvx_mac32/';
elseif mexext == 'mexw64'
    cvx_path = 'cvx/cvx_mac64/';
elseif mexext == 'mexmaci64'
    error('Non-intel based version of Mac not supported by CVX; terminating program') 
else
    error('Cannot detect operating system; terminating program') 
end

%% Check to see if CVX license exists
if exist('cvx/cvx_license/cvx_license.dat', 'file') ~= 2
   error(sprintf('Cannot find cvx_license.dat.   Please request one from http://cvxr.com/cvx/academic/.  Save in the folder %s/cvx/cvx_license',cd))
end

%% Install cvx 
cd(cvx_path)
cvx_setup ../cvx_license/cvx_license.dat
cvx_startup
cd ../..
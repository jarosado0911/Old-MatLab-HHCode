function  [nX,nT, X, T,mthd]=set_sim_params()
    nX = 400;       % number of spatial steps
    nT = 3e5;     % number of time steps
    X = 1000*1e-6;  % length of dendrite segment [m]
    T = 0.035;       % endtime                    [s]
    mthd=0; % set to 0 for FE on ODEs and set to 1 for RK4
end


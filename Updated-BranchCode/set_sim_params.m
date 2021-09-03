function  [nX,nT, X, T,end_loc,brc_ngbr_loc,brc_loc,src_loc,mthd]=set_sim_params()
    nX = 300;
    nT = 20000;
    X = 1;
    T = 25;
    end_loc = floor(nX*2/3)-1;
    brc_ngbr_loc = end_loc+1;
    brc_loc = floor(nX/3);
    src_loc=1;
    mthd=1; % set to 0 for FE on ODEs and set to 1 for RK4
end


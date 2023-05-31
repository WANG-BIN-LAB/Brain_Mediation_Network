pa='H:\Preprocessing\HCP360\REST1\*.mat';
Listn=dir(pa);
parfor k=1:1%j
    A=importdata(Listn(k).name); 
    B = diaomedia(A);
    fname = sprintf( 'B_sub_%d', k );
    parsave( fname, B);
end



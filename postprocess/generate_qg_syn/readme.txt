This folder creates the QG field and synthetic QG field from the simulation output Shafer gave me.
even_extension.m first makes an even extension of the streamfunction. Then it removes high vertical mdoes to satisfiy no-flux boundary condition. 
psi.mat is the original field Shafer gave to me.
psi_filter.mat is the field without high vertical modes.
Scale the streamfunciton in psi_filter.mat to Ro=0.1 leads to the streamfunction in psiN100.mat
generate_randph.m creates the random phase synthetic QG field from psiN100.mat. 

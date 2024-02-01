

function [rho,SH] = exponential_model_data(h)
          [T,~,rho,grav_const,MolecularWeight] = ATM1976(h); 

          % SH 
          R = 8.314 /(10^(-3));  % J / kmol
          
          SH = R*T/(MolecularWeight*grav_const);


end



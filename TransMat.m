function mat = TransMat(phi,beta,d,a)
mat = [cos(phi),-cos(beta)*sin(phi),sin(beta)*sin(phi),a*cos(phi);sin(phi),cos(beta)*cos(phi),-sin(beta)*cos(phi),a*sin(phi);0,sin(beta),cos(beta),d;0,0,0,1];
end
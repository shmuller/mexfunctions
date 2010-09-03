function gpc_tpatch_test

n1 = 10;
n2 = 10;
x = {rand(n1,1),rand(n2,1)};
y = {rand(n1,1),rand(n2,1)};

c = [0,0,1; 0,1,0; 1,0,0];
a = [0.5,0.5,0.5];

sm_tpatch(x,y,c,a);
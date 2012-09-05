v1 = 0;
v2 = 0;
u1 = linspace(-2,2,5);
u2 = linspace(-2,2,5);

z = [0,0,0];

I = dblMaxwmex('one',1,1,z,1,z,int32([1,1]),'MM',v1,v2,u1,u2);


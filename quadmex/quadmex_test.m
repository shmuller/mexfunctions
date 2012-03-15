format long

f = @sin;

y = quadmex(int32(1),int32(8),f,[0,pi])

y_exact = 2


y = quadfunmex(int32([1,2,3]),int32([32,32,32]),'Maxw',[-5,5],[-5,5],[-5,5],0,0,0,1)

y = quadfunmex(int32(1),int32(128),'Maxw_r',[0,10],0,0,0,1)

y_exact = 8*pi/(2*pi)^(3/2)


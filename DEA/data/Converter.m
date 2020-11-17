data=load('C:\Users\Viktor Cheng\Documents\GitHub\OSE\DEA\data\DataAgriculture.mat');
%remeber to change the location
f=fieldnames(data);
 for k=1:size(f,1)
   xlswrite('DataAgriculture.xlsx',data.(f{k}),f{k})
 end
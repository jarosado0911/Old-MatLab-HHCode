function [numPts,Inter_node,End_node]=get_IntrsctNode()

data = readmatrix('Y_shape_1D3D.txt');
shp_data = readmatrix('Y_shape.txt');
shp_data= [shp_data(:,1),shp_data(:,end-1)];
outLst = data;

numPts=length(shp_data(:,1));
for j=1: length(outLst(:,1))
    if outLst(j,2)==0
        Inter_node = outLst(j,1);
        break
    end
end

End_node = Inter_node*2-1;
end
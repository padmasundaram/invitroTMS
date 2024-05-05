
% use even number of layers around the wire core

function [cx_arr, cy_arr] = get_wirecenters(numloops,numlayers,wireradius,coreradius)

Ntotal = round(numloops * (numlayers/2) + (numloops-1) * (numlayers/2));
fprintf('Computing %d centers\n', Ntotal);

cx_arr = zeros(Ntotal,1);
cy_arr = zeros(Ntotal,1);
indx=1;

% figure;
% hold on;

for layer=1:numlayers,
    if mod(layer,2) == 0, % this will have L-1 loops
        fprintf('even layer\n\n');
        for loop = 1:numloops-1,
            cx = coreradius + wireradius + (wireradius*1.73*(layer-1));
            cy = wireradius * 2* loop;
%             fprintf('%f %f\n', cx,cy);
%             plot(cx,cy,'*');
            cx_arr(indx) = cx; cy_arr(indx) = cy; indx = indx+1;
            [xx,yy] = getcirc(cx,cy,wireradius); 
            %plot(xx,yy); %x(i),z(i), a);
        end;
        %lowercorner = [wirediam*(layer-1) + coreradius,
    else
        fprintf('odd layer [%d]\n', layer);
        for loop = 1:numloops,
            cx = coreradius + wireradius + (wireradius*1.73*(layer-1));
            cy = wireradius * (2*loop-1);
%             fprintf('%f %f\n', cx,cy);
%             plot(cx,cy,'*');
            cx_arr(indx) = cx; cy_arr(indx) = cy; indx = indx+1;
            [xx,yy] = getcirc(cx,cy,wireradius); 
            %plot(xx,yy); %x(i),z(i), a);
        end;
    end;
end;


% figure; plot(cx_arr,cy_arr,'*');
% hold on;
% for i=1:Ntotal,
%     [xx,yy] = getcirc(cx_arr(i),cy_arr(i), wireradius);
%     plot(xx,yy);
% end


    


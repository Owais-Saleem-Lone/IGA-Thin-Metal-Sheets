Length=1;
Width=1;
CP(:,:,1) = [0  0
             1  1]*Length;
% y-coordinates
CP(:,:,2) = [-Width/2 Width/2
             -Width/2 Width/2];
% z-coordinates
CP(:,:,3)=[0 0
           0 0];
% Weights
CP(:,:,4) = [1 1
             1 1];
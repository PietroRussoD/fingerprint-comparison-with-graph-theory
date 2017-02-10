%{
File Name:          main.m
Description:        Fingerprint comparison algorithm based on graph theory
Based on the article: Ratha N. K., Pandit V.D., Bolle R. M., Vaish V., “Robust Fingerprint
                      Authentication Using Local Structural Similarity,” in Proc. Workshop on
                      Applications of Computer Vision, pag. 29-34, 2000. 
Author:             Pietro Russo
Start:		          30/03/2010
Last Modify:	      14/04/2010
%}

clear all;

%If other_fingerprint=0, it is selected the other fingerprint, that is similar to the first,
%or it selects the third, that is different to the first.
other_fingerprint=0;

is_print_all=false;

%{
The V set is the set of the minutiae.
The minutiae are represented as nodes, and, every node is represented with the (x,y,theta) triple.
It is builded an Nx3 matrix, where N (rows) are the numbers of the minutiae, and,
3 (columns) are the fields of every node, that contains the positions of the cartesian coordinate x and y,
and the local orientation with theta.
%}

%Values of the pixel's coordinations whre the minutiae is present
%Fingerprint #1
V1=[162 274 1.36 ; 144 246 1.01 ; 135 251 1.1 ; 132 260 1.1 ;
    167 231 1.45 ; 129 308 0.75 ; 108 269 0.85 ; 114 319 0.72 ;
    175 245 0.98 ; 192 221 0.63 ; 165 206 0.68 ; 215 244 0.79 ;
    196 305 1.43 ; 219 228 0.17 ; 169 212 1.1 ; 203 190 0.58];

if other_fingerprint==0
    %Fingerprint #2
    V2=[164 256 1.13 ; 151 227 1.08 ; 139 232 1.18 ; 136 241 1.18 ;
        171 213 0.93 ; 130 288 0.96 ; 114 247 1.06 ; 114 298 0.80 ;
        179 225 1.19 ; 200 200 0.91 ; 172 192 0.84 ; 220 226 0.68 ;
        194 287 1.38 ; 224 211 2.13 ; 176 194 0.82 ; 212 174 0.51];
else
    %Another set of nodes of another images
    %Fingerprint #3
    V2=[199 364 0.87 ; 140 340 1.26 ; 164 226 1.41 ; 170 209 1.43 ;
        217 182 0.54 ; 211 160 0.72 ; 249 178 0.45 ; 244 198 0.40 ;
        243 235 0.73 ; 243 251 1.17 ; 262 294 1.50 ; 263 318 1.06 ;
        271 318 1.43 ; 277 239 2.11 ; 295 277 1.57 ; 308 235 1.79];
end

disp('The V1 matrix (x,y,theta)');
disp(V1);
disp('The V2 matrix (x,y,theta)');
disp(V2);

%The lines matrix is created to keep track of the number of the lines between
%the i-th row and rhe j-th column.
%Number of the papillary crests between the minuties of the fingerprint #1
lines1=[ 0 5 6 6 3 4 8 6 2 2 4 6 7 5 4 4 ;
        5 0 2 2 3 2 4 2 4 5 2 8 9 7 2 4 ;
        6 2 0 2 4 3 3 2 4 6 2 10 11 8 2 4 ;
        6 2 2 0 5 3 3 2 6 7 2 10 12 9 2 5 ;
        3 3 4 5 0 3 6 4 2 3 3 7 7 6 2 2 ;
        4 2 3 3 3 0 4 2 4 5 2 9 10 8 2 3 ;
        8 4 3 3 6 4 0 3 8 9 4 12 14 11 5 7 ;
        6 2 2 2 4 2 3 0 6 6 3 10 10 9 3 6 ;
        2 4 4 6 2 4 8 6 0 2 5 5 6 4 3 2 ;
        2 5 6 7 3 5 9 6 2 0 5 5 6 4 4 2 ;
        4 2 2 2 3 2 4 3 5 5 0 9 10 8 2 3 ;
        6 8 10 10 7 9 12 10 5 5 9 0 3 2 8 7 ;
        7 9 11 12 7 10 14 10 6 6 10 3 0 4 9 8;
        5 7 8 9 6 8 11 9 4 4 8 2 4 0 7 5 ;
        4 2 2 2 2 2 5 3 3 4 2 8 9 7 0 2 ;
        4 4 4 5 2 3 7 6 2 2 3 7 8 5 2 0];

if other_fingerprint==0
    %The number of the papillary crests between the minuties of the fingerprint #2
    %are the same of the fingerprint #1
    lines2=lines1;
else
    %Number of papillary crest between the minuties of the fingerprint #3
    lines2=[ 0 8 10 10 8 10 8 7 4 2 3 4 5 5 7 9 ;
            8 0 3 3 2 3 2 3 5 6 10 11 12 11 14 16 ;
            10 3 0 2 3 2 4 5 8 8 12 13 13 9 16 12 ;
            10 3 2 0 3 2 4 5 8 8 12 13 13 9 16 12 ;
            8 2 3 3 0 3 2 3 6 7 9 10 10 4 8 6 ;
            10 3 2 2 3 0 5 6 8 9 11 13 11 7 10 7 ;
            8 2 4 4 2 5 0 2 6 7 6 7 7 3 3 3 ;
            7 3 5 5 3 6 2 0 4 5 4 6 6 2 4 4 ;
            4 5 8 8 6 8 6 4 0 2 3 3 3 4 6 6 ;
            2 6 8 8 7 9 7 5 2 0 3 5 6 5 8 8 ;
            3 10 12 12 9 11 6 4 3 3 0 2 2 2 5 6 ;
            4 11 13 13 10 13 7 6 3 5 2 0 2 2 5 6;
            5 12 13 13 10 11 7 6 3 6 2 2 0 2 3 4;
            5 11 9 9 4 7 3 2  4 5 2 2 2 0 3 4 ;
            7 14 16 16 8 10 3 4 6 8 2 5 3 3 0 3 ;
            9 16 12 12 6 7 3 4 6 8 6 6 4 4 3 0];
end

%-----Used variables

%Max distance between two nodes
dmax=100;

%Min number of links that the star must to have, to consider it as neighbourhood of the minutiae
Nmin=5;

%Max difference of the consecutive archs
drel=0.19;

%Max difference of the number of rows between two consecutive archs
cdiff=3;

%Max difference of the angle between two coupled archs
Dphi=0.2;

%It is the number of the best couples assosiated between the two fingerprints.
%It is selected from Cmn
Tm=8;

%It is the expected percentage of the total number of the minutiae between two fingerprints
fmin=0.4;

%Min number of the coupled archs of the TOP set. It is the number to guarantee in the consistence phase.
Nm=round(Tm*fmin);

%CONST is a number such that the result is between 0 and 100,
%because the final formula is score=100-CONST*Cstrict*Cconsistency*Cext
CONST=100000;

%-----

%Put the V1 rows into dim1
dim1=length(V1(:,1));
%Initialize the matrix that constitutes the neighbourhood of the indicated node on the row
NEIGHBOURHOOD1 = zeros(1,dim1+1);
%Calculate the neighbourhood matrix, where every row constitutes a star of the i-th row
for i=1:dim1
    for j=1:dim1
        %Calculate the distance between i node and the others.
        %Where V1(i,1) is the x coordinate and V1(i,2) is the y coordinate
        distance = sqrt((V1(i,1)-V1(j,1))^2 + (V1(i,2)-V1(j,2))^2);
        %The obtained distance sould not be bigger than the dmax varible
        if distance~=0 && distance<=dmax
            NEIGHBOURHOOD1(i,j) = distance;
            %In the last row are stored the number of the archs of the node i-th
            NEIGHBOURHOOD1(i,dim1+1) = NEIGHBOURHOOD1(i,dim1+1) + 1;
        end
    end
end
%The i-th row is the start node, the j-th column is the arrival node, with the relative distance.
%In the last column is present the number of the linked archs to the i-th node.
if(is_print_all)
  disp('Neighbourhood 1:');
  disp(NEIGHBOURHOOD1);
end

%Repeat the procedure for V2
dim2=length(V2(:,1));
NEIGHBOURHOOD2 = zeros(1,dim2+1);
for i=1:dim2
    for j=1:dim2
        distance = sqrt((V2(i,1)-V2(j,1))^2 + (V2(i,2)-V2(j,2))^2);
        if distance~=0 && distance<=dmax    
            NEIGHBOURHOOD2(i,j) = distance;
            NEIGHBOURHOOD2(i,dim2+1) = NEIGHBOURHOOD2(i,dim2+1) + 1;
        end
    end
end
if(is_print_all)
  disp('Neighbourhood 2:');
  disp(NEIGHBOURHOOD2);
end

%--------------------
%START OF FIRST PHASE
%--------------------


%{
The E set represents the set of the archs to consider.
It is composed by an Mx5 matrix, where M is the number of the archs and
5 is the number of the fields of every arch (u,e,rad,rc,phi), where u is
the origin node, v the destination node, rad is the distance between the u and the v nodes,
rc is the number of lines between u and v, phi is the subtendend angle from the line with the x-axis.
%}
E1=zeros(1,5,1);
for i=1:length(NEIGHBOURHOOD1(:,1))
    %Select the neighourhood of nodes with at least Nmin archs that composes the star.
    k=0;
    if NEIGHBOURHOOD1(i,end) >= Nmin       
        for j=1:length(NEIGHBOURHOOD1(1,:))-1;
            if (i~=j) && (NEIGHBOURHOOD1(i,j)~=0)
                k=k+1;
                %Store on the first column the start node
                E1(k,1,i)=i;
                %Store on the second column the arrival node
                E1(k,2,i)=j;
                %Store on the third colum the distance between the two nodes
                E1(k,3,i)=NEIGHBOURHOOD1(i,j);
                %Store on the fourth column the number of the lines between two nodes
                E1(k,4,i)=lines1(i,j);
                %Calculate the phi angle.
                %V1(i,1) is the x coordinate of the referred node
                %V1(i,2) is the y coordinate of the referred node
                %V1(j,1) is the x coordinate of the near node
                %V1(j,2) is the y coordinate of the near node
                y=V1(j,2)-V1(i,2);
                x=V1(j,1)-V1(i,1);
                phi=atan(y/x);
                %The angle is set in fucntion of the quadrant, respect to the referred node
                if x>0 && y>0
                    E1(k,5,i)=phi;                
                elseif x<0 && y>0
                    E1(k,5,i)=pi-abs(phi);
                elseif x<0 && y<0
                    E1(k,5,i)=pi+phi;
                elseif x>0 && y<0
                    E1(k,5,i)=2*pi-abs(phi);
                end  
            end
        end
    end
end
%Order clockwise, according to the angle that is formed respect to the reference node,
%and then reordering the first column in function of the fifth column
E1_temp=zeros(length(E1(:,1,1)),5);
for m=1:length(E1(1,1,:))
    E1_temp(:,:,m)=sortrows(E1(:,:,m),[1 5]);    
end
%Put the null elements in the last rows
E1=zeros(length(E1(:,1,1)),5);
for m=1:length(E1_temp(1,1,:))
    k=0;
    for i=1:length(E1_temp(:,1,1))
        if E1_temp(i,:,m)~=zeros(1,length(E1_temp(1,:,1)))
            k=k+1;
            E1(k,:,m)=E1_temp(i,:,m);
        end
    end
end
if(is_print_all)
  disp('Set E1:');
  disp('       m        u_i    distance   n_linee  angle_x');
  disp(E1);
end

%Repeat the process for the other matrix
E2=zeros(1,5,1);
for i=1:length(NEIGHBOURHOOD2(:,1))
    k=0;
    if NEIGHBOURHOOD2(i,end) >= Nmin       
        for j=1:length(NEIGHBOURHOOD2(1,:))-1;
            if (i~=j) && (NEIGHBOURHOOD2(i,j)~=0)
                k=k+1;
                E2(k,1,i)=i;
                E2(k,2,i)=j;
                E2(k,3,i)=NEIGHBOURHOOD2(i,j);
                E2(k,4,i)=lines2(i,j);
                y=V2(j,2)-V2(i,2);
                x=V2(j,1)-V2(i,1);
                phi=atan(y/x);
                if x>0 && y>0
                    E2(k,5,i)=phi;                
                elseif x<0 && y>0
                    E2(k,5,i)=pi-abs(phi);
                elseif x<0 && y<0
                    E2(k,5,i)=pi+phi;
                elseif x>0 && y<0
                    E2(k,5,i)=2*pi-abs(phi);
                end  
            end
        end
    end
end

E2_temp=zeros(length(E2(:,1,1)),5);
for m=1:length(E2(1,1,:))
    E2_temp(:,:,m)=sortrows(E2(:,:,m),[1 5]);    
end

E2=zeros(length(E2(:,1,1)),5);
for m=1:length(E2_temp(1,1,:))
    k=0;
    for i=1:length(E2_temp(:,1,1))
        if E2_temp(i,:,m)~=zeros(1,length(E2_temp(1,:,1)))
            k=k+1;
            E2(k,:,m)=E2_temp(i,:,m);
        end
    end
end
if(is_print_all)
  disp('Set E2:');
  disp('       m        u_i    distance   n_linee  angle_x');
  disp(E2);
end

%Determine the coupled arches
%The B matrix has as columns: the reference node, the two linked nodes to the reference node,
%the angle that is formed between two arches and the radial cost to couple the two arches 
B1=zeros(1,5,1);
%For every arch in E
for m=1:(length(E1(1,1,:)))
    k=0;
    for i=1:(length(E1(:,1,1))-1)
        %If it is not the second last element
        if i~=length(E1(:,1,1))-1 && E1(i,1,m)~=0    
            %If the i-th row of E1 and that one immediately after are the same, the nodes are consecutive.
            if E1(i,1,m)==E1(i+1,1,m) && E1(i,1,m)~=0
                %Variables useful for the comparison
                d1=abs(E1(i,3,m)-E1(i+1,3,m))/min(E1(i,3,m),E1(i+1,3,m));   
                c1=abs(E1(i,4,m)-E1(i+1,4,m));
                %Two arches are coupled if they respect the following condition
                if d1<=drel && c1<=cdiff
                    k=k+1;
                    %Reference node
                    B1(k,1,m)=E1(i,1,m);
                    %First linked node
                    B1(k,2,m)=E1(i,2,m);
                    %Second linked node
                    B1(k,3,m)=E1(i+1,2,m);   
                    %Angle between two arches
                    B1(k,4,m)=E1(i+1,5,m)-E1(i,5,m);
                    %Radial cost to couple two arches
                    B1(k,5,m)=d1/drel;
                end 
            %If they are different, the followed row is full of zeros and
	    %it is necessary to check if the last arch (i-th) and the first are coupleable
            elseif E1(i,1,m)~=E1(i+1,1,m) && E1(i,1,m)~=0 && E1(i+1,1,m)==0
                d1=abs(E1(i,3,m)-E1(1,3,m))/min(E1(i,3,m),E1(1,3,m));   
                c1=abs(E1(i,4,m)-E1(1,4,m));
                if d1<=drel && c1<=cdiff
                    k=k+1;
                    B1(k,1,m)=E1(i,1,m);
                    B1(k,2,m)=E1(i,2,m);
                    B1(k,3,m)=E1(1,2,m);   
                    B1(k,4,m)=2*pi-(E1(i,5,m)-E1(1,5,m));
                    B1(k,5,m)=d1/drel;
                end           
            end
        %If it is the second last element and it is not equal to zero
        elseif i==length(E1(:,1,1))-1 && E1(i,1,m)~=0 && E1(i+1,1,m)~=0
            %Compare it with the last one
            d1=abs(E1(i,3,m)-E1(i+1,3,m))/min(E1(i,3,m),E1(i+1,3,m));
            c1=abs(E1(i,4,m)-E1(i+1,4,m));
            if d1<=drel && c1<=cdiff
                k=k+1;            
                B1(k,1,m)=E1(i,1,m);
                B1(k,2,m)=E1(i,2,m);
                B1(k,3,m)=E1(i+1,2,m);   
                B1(k,4,m)=E1(i+1,5,m)-E1(i,5,m);
                B1(k,5,m)=d1/drel;
            end
            %Compare the last with the first
            d1=abs(E1(i+1,3,m)-E1(1,3,m))/min(E1(i+1,3,m),E1(1,3,m));   
            c1=abs(E1(i+1,4,m)-E1(1,4,m));
            if d1<=drel && c1<=cdiff
                k=k+1;
                B1(k,1,m)=E1(i+1,1,m);
                B1(k,2,m)=E1(i+1,2,m);
                B1(k,3,m)=E1(1,2,m);   
                B1(k,4,m)=2*pi-(E1(i+1,5,m)-E1(1,5,m));
                B1(k,5,m)=d1/drel;
            end
        %If the second last is not equal to zero, but thte last is equal to zero,
	%it is compared the second last with the first
        elseif i==length(E1(:,1,1))-1 && E1(i,1,m)~=0 && E1(i+1,1,m)==0
            d1=abs(E1(i,3,m)-E1(1,3,m))/min(E1(i,3,m),E1(1,3,m));
            c1=abs(E1(i,4,m)-E1(1,4,m));
            if d1<=drel && c1<=cdiff
                k=k+1;            
                B1(k,1,m)=E1(i,1,m);
                B1(k,2,m)=E1(i,2,m);
                B1(k,3,m)=E1(1,2,m);   
                B1(k,4,m)=2*pi-(E1(i,5,m)-E1(1,5,m));
                B1(k,5,m)=d1/drel;
            end            
        end
    end
end
if(is_print_all)
  disp('B1:');
  disp('       m        u_i      u_i+1    angle  rad_cost');
  disp(B1);
end

B2=zeros(1,5,1);
%Determine the coupled arches of the other matrix
for m=1:(length(E2(1,1,:)))
    k=0;
    for i=1:(length(E2(:,1,1))-1)
        if i~=length(E2(:,1,1))-1 && E2(i,1,m)~=0    
            if E2(i,1,m)==E2(i+1,1,m) && E2(i,1,m)~=0
                d2=abs(E2(i,3,m)-E2(i+1,3,m))/min(E2(i,3,m),E2(i+1,3,m));   
                c2=abs(E2(i,4,m)-E2(i+1,4,m));
                 if d2<=drel && c2<=cdiff
                    k=k+1;
                    B2(k,1,m)=E2(i,1,m);
                    B2(k,2,m)=E2(i,2,m);
                    B2(k,3,m)=E2(i+1,2,m);   
                    B2(k,4,m)=E2(i+1,5,m)-E2(i,5,m);
                    B2(k,5,m)=d2/drel;
                end 
           elseif E2(i,1,m)~=E2(i+1,1,m) && E2(i,1,m)~=0 && E2(i+1,1,m)==0
                d2=abs(E2(i,3,m)-E2(1,3,m))/min(E2(i,3,m),E2(1,3,m));   
                c2=abs(E2(i,4,m)-E2(1,4,m));
                if d2<=drel && c2<=cdiff
                    k=k+1;
                    B2(k,1,m)=E2(i,1,m);
                    B2(k,2,m)=E2(i,2,m);
                    B2(k,3,m)=E2(1,2,m);   
                    B2(k,4,m)=2*pi-(E2(i,5,m)-E2(1,5,m));
                    B2(k,5,m)=d2/drel;
                end           
            end
        elseif i==length(E2(:,1,1))-1 && E2(i,1,m)~=0 && E2(i+1,1,m)~=0
            d2=abs(E2(i,3,m)-E2(i+1,3,m))/min(E2(i,3,m),E2(i+1,3,m));
            c2=abs(E2(i,4,m)-E2(i+1,4,m));
            if d2<=drel && c2<=cdiff
                k=k+1;            
                B2(k,1,m)=E2(i,1,m);
                B2(k,2,m)=E2(i,2,m);
                B2(k,3,m)=E2(i+1,2,m);   
                B2(k,4,m)=E2(i+1,5,m)-E2(i,5,m);
                B2(k,5,m)=d2/drel;
            end
            d2=abs(E2(i+1,3,m)-E2(1,3,m))/min(E2(i+1,3,m),E2(1,3,m));   
            c2=abs(E2(i+1,4,m)-E2(1,4,m));
            if d2<=drel && c2<=cdiff
                k=k+1;
                B2(k,1,m)=E2(i+1,1,m);
                B2(k,2,m)=E2(i+1,2,m);
                B2(k,3,m)=E2(1,2,m);   
                B2(k,4,m)=2*pi-(E2(i+1,5,m)-E2(1,5,m));
                B2(k,5,m)=d2/drel;
            end
        elseif i==length(E2(:,1,1))-1 && E2(i,1,m)~=0 && E2(i+1,1,m)==0
            d2=abs(E2(i,3,m)-E2(1,3,m))/min(E2(i,3,m),E2(1,3,m));
            c2=abs(E2(i,4,m)-E2(1,4,m));
            if d2<=drel && c2<=cdiff
                k=k+1;            
                B2(k,1,m)=E2(i,1,m);
                B2(k,2,m)=E2(i,2,m);
                B2(k,3,m)=E2(1,2,m);   
                B2(k,4,m)=2*pi-(E2(i,5,m)-E2(1,5,m));
                B2(k,5,m)=d2/drel;
            end            
        end
    end
end
if(is_print_all)
  disp('B2:');
  disp('       m        u_i      u_i+1    angle  rad_cost');
  disp(B2);
end

%Calculation of angular cost, determined on the basis of all combinations
%of consecutive arches of the two matrices representing different fingerprints
dphi=zeros(1,7,1);
h=0;
%Angular cost to couple two consecutive arches of the m node with those od the n node
for m=1:length(B1(1,1,:))
    for n=1:length(B2(1,1,:))
        if B1(1,1,m)~=0 && B2(1,1,n)~=0
            h=h+1;    
            k=0;
            for i=1:length(B1(:,1,1));
                for j=1:length(B2(:,1,1))                
                    if B1(i,1,m)~=0 && B2(j,1,n)~=0
                        phi=abs(B1(i,4,m)-B2(j,4,n));
                        if phi <= Dphi
                            k=k+1;
                            %Reference node of the first image
                            dphi(k,1,h)=B1(i,1,m);
                            %Node near to the reference node of the first image
                            dphi(k,2,h)=B1(i,2,m);
                            %Another node near to the reference node of the firs image
                            dphi(k,3,h)=B1(i,3,m);
                            %Radial cost to couple two consecutive arches of the first fingerprint
                            dphi(k,4,h)=B1(i,5,m);
                            %Reference node of the second image
                            dphi(k,5,h)=B2(j,1,n);
                            %Node near to the reference node of the second image
                            dphi(k,6,h)=B2(j,2,n);
                            %Another node near to the reference node of the second image
                            dphi(k,7,h)=B2(j,3,n);
                            %Radial cost to couple two consecutive arches of the second fingerprint
                            dphi(k,8,h)=B2(j,5,n);
                            %Angular cost to couple consecutive arches of the two fingerprints
                            dphi(k,9,h)=abs(B1(i,4,m)-B2(j,4,n))/Dphi;                            
                        end
                    end
                end
            end
        end
    end
end
if(is_print_all)
  disp('dphi:');
  disp('       m        ui       ui+1        n        vi       vi+1      dphi     drad');
  disp(dphi);
end

%Cost to couple m with n
Cmn=zeros(1,3);
ne1=zeros(1,3);
ne2=zeros(1,3);
k=0;

for m=1:length(dphi(1,1,:))
    sum_dphi=0;
    rows=0;
    
    %Sum all angular costs dphi of the coupled arches
    for i=1:length(dphi(:,1,1))
        if dphi(i,1,m)~=0
            sum_dphi=sum_dphi+dphi(i,9,m);            
            rows=rows+1;
        end
    end
    
    if dphi(1,1,m)~=0
        %Determine Ne with the selection of the unique coupled arches
        ne1_temp=dphi(1:rows,1:4,m);
        ne2_temp=dphi(1:rows,5:8,m);
    
        ne1_temp=sortrows(ne1_temp,[1 2]);
        ne2_temp=sortrows(ne2_temp,[1 2]);
        
        ne1=zeros(1,4);
        ne2=zeros(1,4);
        prev_ne1=zeros(1,3);
        prev_ne2=zeros(1,3);
        
        h1=0;
        h2=0;
               
        for i=1:rows
            if prev_ne1(1,2:3)~=ne1_temp(i,2:3)
                h1=h1+1;
                ne1(h1,1:4)=ne1_temp(i,1:4); 
            end        
            prev_ne1=ne1_temp(i,1:3); 
            
            if prev_ne2(1,2:3)~=ne2_temp(i,2:3)
                h2=h2+1;
                ne2(h2,1:4)=ne2_temp(i,1:4); 
            end        
            prev_ne2=ne2_temp(i,1:3);
        end
        
        Ne=0;
        sum_dr=0;
        
        for i=1:length(ne1(:,1))
            sum_dr=sum_dr+ne1(i,4);
            Ne=Ne+1;
        end
        
        for i=1:length(ne2(:,1))
            sum_dr=sum_dr+ne2(i,4);
            Ne=Ne+1;
        end
        %Determine the cost of every couple of nodes
        k=k+1;
        Cmn(k,1)=dphi(1,1,m);
        Cmn(k,2)=dphi(1,5,m);
        Cmn(k,3)=(sum_dr+sum_dphi)/(Ne)^2; 
        
    end        
 end
%Sort the matrix in ascending order, on the basis of the cost value, because after we take the lowest values.
Cmn=sortrows(Cmn,3);
if(is_print_all)
  disp(Cmn);
end
%Take the unique couples with the lowest value
Cmn_temp=Cmn;
TOP_temp=Cmn_temp(1,:);
k=1;
fine=0;
while fine==0
    for i=1:length(Cmn_temp)    
        if TOP_temp(k,1)==Cmn_temp(i,1) || TOP_temp(k,2)==Cmn_temp(i,2)
            Cmn_temp(i,3)=999;            
        end
    end
    Cmn_temp=sortrows(Cmn_temp,3);
    if Cmn_temp(1,3)~=999
        k=k+1;
        TOP_temp(k,:)=Cmn_temp(1,:);
    else
        fine=1;
    end    
end
%Take the first Tm elements
TOP=TOP_temp(1:Tm,:);

if(is_print_all)
  disp(TOP);
end

%Determine the cost of the first phase
Cstrict=sum(TOP(1:end,3))/Tm;
if(is_print_all)
  disp(Cstrict);
end

%---------------------
%START OF SECOND PHASE
%---------------------


%Start the phase that checks the consistency
%Select the nodes of the TOP set
Q1=TOP(:,1);
Q1=sortrows(Q1,1);
if(is_print_all)
  disp(Q1);
end

Q2=TOP(:,2);
Q2=sortrows(Q2,1);
if(is_print_all)
  disp(Q2);
end

%Select the arches of the nodes of the TOP set
Ec1=zeros(1,2,length(E1(1,1,:)));
comb1=combnk(Q1(:,1),2);
comb2=[comb1(:,2) comb1(:,1)];
comb=[comb1; comb2];
comb=sortrows(comb,1);
dim_prec=comb(1,1);
k=0;
for i=1:length(comb(:,1))
    if comb(i,1)==dim_prec    
        k=k+1;    
        Ec1(k,1,comb(i,1))=comb(i,1);
        Ec1(k,2,comb(i,1))=comb(i,2);
    else
        k=0;
        k=k+1;    
        Ec1(k,1,comb(i,1))=comb(i,1);
        Ec1(k,2,comb(i,1))=comb(i,2);                   
    end
    dim_prec=comb(i,1);
end
for m=1:(length(Ec1(1,1,:)))
    for i=1:length(Ec1(:,1,1))
        for j=1:length(E1(:,1,1))
            if Ec1(i,2,m)==E1(j,2,m)
                Ec1(i,3,m)=E1(j,3,m);
                Ec1(i,4,m)=E1(j,4,m);
            end
        end
    end
end
if(is_print_all)
  disp('Selected arches 1');
  disp(Ec1);
end

Ec2=zeros(1,2,length(E2(1,1,:)));
comb1=combnk(Q2(:,1),2);
comb2=[comb1(:,2) comb1(:,1)];
comb=[comb1; comb2];
comb=sortrows(comb,1);
dim_prec=comb(1,1);
k=0;
for i=1:length(comb(:,1))
    if comb(i,1)==dim_prec    
        k=k+1;    
        Ec2(k,1,comb(i,1))=comb(i,1);
        Ec2(k,2,comb(i,1))=comb(i,2);
    else
        k=0;
        k=k+1;    
        Ec2(k,1,comb(i,1))=comb(i,1);
        Ec2(k,2,comb(i,1))=comb(i,2);                   
    end
    dim_prec=comb(i,1);
end
for m=1:(length(Ec2(1,1,:)))
    for i=1:length(Ec2(:,1,1))
        for j=1:length(E2(:,1,1))
            if Ec2(i,2,m)==E2(j,2,m)
                Ec2(i,3,m)=E2(j,3,m);
                Ec2(i,4,m)=E2(j,4,m);
            end
        end
    end
end
if(is_print_all)
  disp('ASelected arches 2');
  disp(Ec2);
end

%It starts to couple the arches of Ec1 with the archs of Ec2
h=0;
Bcons=zeros(1,5,1);
for m=1:(length(Ec1(1,1,:)))
    for n=1:(length(Ec2(1,1,:)))
        if Ec1(1,1,m)~=0 && Ec2(1,1,n)~=0
            h=h+1;  
            k=0;
            for i=1:length(Ec1(:,1,1))             
                for j=1:length(Ec2(:,1,1))
                    if Ec1(i,3,m)~=0 && Ec2(j,3,n)~=0
                        d=abs(Ec1(i,3,m)-Ec2(j,3,n))/min(Ec1(i,3,m),Ec2(j,3,n));
                        c=abs(Ec1(i,4,m)-Ec2(j,4,n));
                        if d<=drel && c<=cdiff
                            k=k+1;
                            Bcons(k,1,h)=Ec1(i,1,m);
                            Bcons(k,2,h)=Ec1(i,2,m);
                            Bcons(k,3,h)=Ec2(j,1,n);
                            Bcons(k,4,h)=Ec2(j,2,n);
                            Bcons(k,5,h)=d;
                        end
                    end                
                end
            end
        end
    end
end
if(is_print_all)
  disp('Coupled arches');
  disp(Bcons);
end

%It stores the max number of the coupled arches
N1=length(Bcons(:,1,1));
%Calculate the number of consistent couples, that have at least Nm arches coupled
%For every couple of nodes, it stores how many coupled arches there are
N2=zeros(length(Bcons(1,1,:)),1);
for m=1:length(Bcons(1,1,:))
    for i=1:length(Bcons(:,1,1))
        if Bcons(i,1,m)~=0
            N2(m,1)=N2(m,1)+1;            
        end
    end
end
%Calculate the averages per every coupled node
ccons_pair=zeros(1);
k=0;
for i=1:length(N2(:,1))
    if N2(i)>=Nm
        k=k+1;
        ccons_pair(k,1)=Bcons(1,1,i);
        ccons_pair(k,2)=Bcons(1,3,i);  
        ccons_pair(k,3)=mean(Bcons(1:N2(i),5,i));        
    end    
end
ccons_pair=sortrows(ccons_pair, 3);
%Select the distinct couple of nodes with the best value
ccons_pair_temp=ccons_pair;
TOP_ccons_pair=ccons_pair(1,:);
k=1;
fine=0;
while fine==0
    for i=1:length(ccons_pair_temp)    
        if TOP_ccons_pair(k,1)==ccons_pair_temp(i,1) || TOP_ccons_pair(k,2)==ccons_pair_temp(i,2)
            ccons_pair_temp(i,3)=999;            
        end
    end
    ccons_pair_temp=sortrows(ccons_pair_temp,3);
    if ccons_pair_temp(1,3)~=999
        k=k+1;
        TOP_ccons_pair(k,:)=ccons_pair_temp(1,:);
    else
        fine=1;
    end    
end
ccons_pair=TOP_ccons_pair;
%Check the consistency
if N1<Nm && length(ccons_pair(:,1))<Nm
    disp('It is not consistent');
else
    disp('It is consistent');
    %Cost of the second phase
    Cconsistency=sum(ccons_pair(:,3))/length(ccons_pair(:,3));
    if(is_print_all)
      disp('Cconsistency');
      disp(Cconsistency);
    end
        
%-----------------
%START OF THIRD PHASE
%-----------------


  %Select nodes that are not in the TOP set
  Qe1=zeros(1);
  for i=1:length(V1(:,1))
      Qe1(i,1)=i;
  end
  for i=1:length(Q1(:,1))
      for j=1:length(Qe1(:,1))        
          if Q1(i,1)==Qe1(j,1)
              Qe1(j,1)=0;            
          end
      end
  end
  Qe1_temp=sortrows(Qe1, 1);
  Qe1=zeros(1);
  k=0;
  for i=1:length(Qe1_temp(:,1))
      if Qe1_temp(i,1)~=0
          k=k+1;
          Qe1(k,1)=Qe1_temp(i,1);
      end        
  end
  if(is_print_all)
    disp(Qe1);
  end
  
  Qe2=zeros(1);
  for i=1:length(V2(:,1))
      Qe2(i,1)=i;
  end
  for i=1:length(Q2(:,1))
      for j=1:length(Qe2(:,1))        
          if Q2(i,1)==Qe2(j,1)
              Qe2(j,1)=0;            
          end
      end
  end
  Qe2_temp=sortrows(Qe2, 1);
  Qe2=zeros(1);
  k=0;
  for i=1:length(Qe2_temp(:,1))
      if Qe2_temp(i,1)~=0
          k=k+1;
          Qe2(k,1)=Qe2_temp(i,1);
      end        
  end
  if(is_print_all)
    disp(Qe2);
  end
  
  %Determine the formed arches with the extracted nodes previously
  Ee1=zeros(1,2,length(E1(1,1,:)));
  for i=1:length(Qe1(:,1))
      k=0;
      for j=1:length(Q1(:,1))
          k=k+1;
          Ee1(k,1,Qe1(i,1))=Qe1(i,1);
          Ee1(k,2,Qe1(i,1))=Q1(j,1);            
      end
  end
  %Fill the data of te coupled nodes, that was calculated in the first phase
  for m=1:(length(Ee1(1,1,:)))
      for i=1:length(Ee1(:,1,1))
          for j=1:length(E1(:,1,1))
              if Ee1(i,2,m)==E1(j,2,m)
                  Ee1(i,3,m)=E1(j,3,m);
                  Ee1(i,4,m)=E1(j,4,m);
              end
          end
      end
  end
  if(is_print_all)
    disp(Ee1);
  end
  
  Ee2=zeros(1,2,length(E2(1,1,:)));
  for i=1:length(Qe2(:,1))
      k=0;
      for j=1:length(Q2(:,1))
          k=k+1;
          Ee2(k,1,Qe2(i,1))=Qe2(i,1);
          Ee2(k,2,Qe2(i,1))=Q2(j,1);            
      end
  end
  for m=1:(length(Ee2(1,1,:)))
      for i=1:length(Ee2(:,1,1))
          for j=1:length(E2(:,1,1))
              if Ee2(i,2,m)==E2(j,2,m)
                  Ee2(i,3,m)=E2(j,3,m);
                  Ee2(i,4,m)=E2(j,4,m);
              end
          end
      end
  end
  if(is_print_all)
    disp(Ee2);
  end
  
  %Determina the coupled arches between two sets of nodes
  Bext=zeros(1,5,1);
  h=0;
  for m=1:(length(Ee1(1,1,:)))
      for n=1:(length(Ee2(1,1,:)))
          if (Ee1(1,1,m)~=0 && Ee2(1,1,n)~=0)
              h=h+1;  
              k=0;
              for i=1:length(Ee1(:,1,1))             
                  for j=1:length(Ee2(:,1,1))
                      if (Ee1(i,3,m)~=0 && Ee2(j,3,n)~=0)
                          d=abs(Ee1(i,3,m)-Ee2(j,3,n))/min(Ee1(i,3,m),Ee2(j,3,n));
                          c=abs(Ee1(i,4,m)-Ee2(j,4,n));
                          if (d<=drel && c<=cdiff)
                              k=k+1;
                              Bext(k,1,h)=Ee1(i,1,m);
                              Bext(k,2,h)=Ee1(i,2,m);
                              Bext(k,3,h)=Ee2(j,1,n);
                              Bext(k,4,h)=Ee2(j,2,n);
                              Bext(k,5,h)=d;
                          end
                      end                
                  end
              end
          end
      end
  end
  if(is_print_all)
    disp('Coupled arches');
    disp(Bext);
  end

  %Calculate the number of coupled arches per couple
  Nm_ext=zeros(length(Bext(1,1,:)),1);
  for m=1:length(Bext(1,1,:))
      for i=1:length(Bext(:,1,1))
          if Bext(i,1,m)~=0
              Nm_ext(m,1)=Nm_ext(m,1)+1;            
          end
      end
  end
  %Determina the cost to couple two nodes
  cext=zeros(1);
  k=0;
  for i=1:length(Nm_ext(:,1))
      if Nm_ext(i,1)>=Nm
          k=k+1;
          cext(k,1)=Bext(1,1,i);
          cext(k,2)=Bext(1,3,i);        
          cext(k,3)=mean(Bext(1:Nm_ext(i),5,i));        
      end    
  end
  cext=sortrows(cext,3);
  %Select the couple with the best cost
  cext_temp=cext;
  TOP_cext=cext_temp(1,:);
  k=1;
  fine=0;
  while fine==0
      for i=1:length(cext_temp)    
          if TOP_cext(k,1)==cext_temp(i,1) || TOP_cext(k,2)==cext_temp(i,2)
              cext_temp(i,3)=999;            
          end
      end
      cext_temp=sortrows(cext_temp,3);
      if cext_temp(1,3)~=999
          k=k+1;
          TOP_cext(k,:)=cext_temp(1,:);
      else
          fine=1;
      end    
  end

  cext_nmat=TOP_cext;
  if(is_print_all)
    disp(cext)
  end

  %Expected number of nodes coupled in the extension phase.
  Next=round(fmin*min(length(V1(:,1)),length(V2(:,1))));
  if(is_print_all)
    disp(Next);
  end
  
  %Cost of the extension phase
  Cext=(sum(cext_nmat(1:end,3))*Next)/(length(cext_nmat(:,1)))^2;
  if(is_print_all)
    disp(Cext);
  end
  
  %Cost to couple two graphs
  score=100-CONST*Cstrict*Cconsistency*Cext;
  disp(score);

end %% end of the condition of the consistency

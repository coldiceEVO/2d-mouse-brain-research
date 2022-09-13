imgIdx=1;
dirList=dir;  % getting current list of sub directories with images inside
% for each sub directories
for i=3: size(dirList) % 3:ignoring ./..
    % if not a directory, skip
    if dirList(i).isdir == 0 | startsWith(dirList(i).name, "output")
        continue
    end
    
    cd (dirList(i).name) % change to current directories
    % for each image inside subdirectories
    imgList=dir('*.tif');
    for j=1:size(imgList)
        BW=opengauss(imgList(j).name); %open blur
        [r,c,d]=size(BW);
        CC = bwconncomp(BW) ;%connected component
        %L=labelmatrix(CC);
        %EuDist = bwdist(~BW);
        %Skele = bwmorph(BW,'skel',inf); %this takes too long
        CH=bwconvhull(BW);
        stats = regionprops(CH,"Centroid","ConvexHull");
        xcenter=stats.Centroid(1);
        [bst1,bst2]=symm(BW,xcenter);
        
        CHpoly=stats.ConvexHull;
        %[xv,yv] = polyxpoly([bst1,bst2],[0 r],CHpoly(:,1),CHpoly(:,2));
        %[xv,yv] = polyxpoly([bst1,bst2],[0 r],CHpoly(:,1),CHpoly(:,2));
        [inV,outV] = intersect(polyshape(CHpoly),[bst1 0;bst2 r]);
        inV=max(1,min(r,inV)); %inV=inV(:,2);% vertical contact of convext hull with vertical sysmmetry
        mid=[(bst1+bst2)/2 (inV(1,2)+inV(2,2))/2]; dx_dy=(bst2-bst1)/(r);
        [inH,outH] = intersect(polyshape(CHpoly),[0 mid(2)*(1+dx_dy);c mid(2)*(1-dx_dy)]);% convex hull contact with midpoint othogonal line
        inH=max(1,min(c,inH)); %inH=inH(:,1);
        
        t=-tan(dx_dy)-pi;%CCW  rotation
        %line creation
        % vec=[0.2 0.7; 0.4 -0.1];vec=[vec;vec*[-1 0;0 1]];vec(:,3)=1;% [x;y;1]
        % bvorig=[0.6 0.5;0.1,0.4;-0.1,0.4;-0.6 0.5]; bvorig(:,3)=1; %%%%
        % border=[1,-0.5;1 1;0 -.5;0 1];border=[border;border*[-1 0;0 1]]; border(:,3)=1;
        % Coordinates [x;y;1] (right sided)
        Vs =[.2 .7 1;.4 -.3 1]; %split line -TopDown
        Vm =[.9,1,1;1,-.5,1];Vm=flip(Vm,1); %outside boundaries of MCA origin -BottomUP
        Vms=[.4,.9,1;.6,-.4,1]; %inside boundaires of MCA origin (near split line) -TopDown
        Va =[-.1,.8,1;-.2,-.4,1]; Va=flip(Va,1);%inside boundaries of A -BottomUp
        Vas=[.1,.8,1;.2,-.4,1]; %outside boundaries of ACA origin --TopDown
        V1={Vs;Vm;Vms;Va;Vas};
        V2=V1;
        for i=1:5 %mirroring v2
            tmp=V2{i};
            tmp(:,1)=-tmp(:,1);
            V2{i}=tmp;
        end
        VC={V1;V2};
        
%         % print hour glass for axis
%         BV=[1 1 1;-1 1 1;0 0 1;-1 -1 1; 1 -1 1;0 0 1];
%         BV2=(tf*BV')';
%         BV2=BV2(:,1:2);
%         imshow(I);hold on; plot(polyshape(BV2(:,1),BV2(:,2))); hold off;
%         figure; plot(polyshape(BV(:,1),BV(:,2)));
       
        % coordinates are transformed to the other side as mirror 

        sx=norm(mid-inV(1,:));sy=max(norm(mid-inH(1,:)),norm(mid-inH(2,:)));
        tf=[cos(t) -sin(t) mid(1);sin(t) cos(t) mid(2);0 0 1]*diag([sx,sy,1]);
        %vec=(tf*vec')'; 
        %bvorigtf=(tf*(bvorig'))'; bvorigtf(:,3)=20;%circle size
        %border=(tf*(border'))';
        %vec=vec'; vec=vec(:,1:2); vec=[vec([1,3],:) vec([2,4],:)];
        %msk=poly2mask([],[],r,c)
        
%         Ts =(tf*Vs')';
%         Tm =(tf*Vm')';
%         Tms=(tf*Vms')';
%         Ta =(tf*Va')';
%         Tas=(tf*Vas')';
        
        I=read(Tiff(imgList(j).name,'r'));
        if size(I,3)>3
            I=I(:,:,1:3); %remove extra dimension, really weird
        end
        K=I;
        for k= 1:2
            V=VC{k};
            for i=1:size(V,1)
                tmp=(tf*V{i}')';%transformation
                T{i}=tmp(:,1:2);%remove last col
            end
            ep=T{2};
            mask{1}=[T{2};T{1};ep(1,:)];%Outer CC bound
            mask{2}=[T{2};T{3};ep(1,:)];%Outer MCA region
            ep=T{4};
            mask{3}=[T{4};T{1};ep(1,:)];%Innter CC bound
            mask{4}=[T{4};T{5};ep(1,:)];%Inner ACA region
            ofst=[-1,1];
                    f=figure;
        %imshow(BW);hold on; % plot in image
        
        imshow(K);hold on;
        plot(mask{1}(:,1),mask{1}(:,2),'r','LineWidth',2);
        plot(mask{2}(:,1),mask{2}(:,2),'m','LineStyle','--','LineWidth',2);
        plot(mask{3}(:,1),mask{3}(:,2),'g','LineStyle','-.','LineWidth',2);
        plot(mask{4}(:,1),mask{4}(:,2),'b','LineStyle','--','LineWidth',2);
        legend('MCA Flow','MCA Source','ACA Flow','ACA Source');
        hold off;
% 
%         %hold on; %plot in normalzied space
%         plot([Vs(:,1);Vm(:,1);Vs(1,1)],[Vs(:,2);Vm(:,2);Vs(1,2)],'r');
%         plot([Vms(:,1);Vm(:,1);Vms(1,1)],[Vms(:,2);Vm(:,2);Vms(1,2)],'m','LineStyle','--');
%         plot([Vs(:,1);Va(:,1);Vs(1,1)],[Vs(:,2);Va(:,2);Vs(1,2)],'g','LineStyle','-.');
%         plot([Vas(:,1);Va(:,1);Vas(1,1)],[Vas(:,2);Va(:,2);Vas(1,2)],'b','LineStyle','--');
%         legend('MCA Flow','MCA Source','ACA Flow','ACA Source');
%         hold off;
        exportgraphics(f,strcat('../output5RegionPoly/',imgList(j).name,int2str(k),'.png'))
        close; 
    %         for i=1: size(mask,2)
    %             figure(i);
    %             m=mask{i};imshow(poly2mask(m(:,1),m(:,2),r,c));
    %         end
            m=mask{1};CC1=bwconncomp(BW.*poly2mask(m(:,1),m(:,2),r,c));
            CH1=regionprops(CC1,"ConvexHull");
            I1=false(r,c);
            for i=1:CC1.NumObjects
                m=mask{2};p1=polyshape(m);
                p2=polyshape(CH1(i).ConvexHull);
                if intersect(p1,p2).NumRegions~=0
                    tmp=CC1.PixelIdxList{i};
                    I1(tmp)=1;
                end
            end
            m=mask{3};CC2=bwconncomp(BW.*poly2mask(m(:,1),m(:,2),r,c));
            CH2=regionprops(CC2,"ConvexHull");
            I2=false(r,c);
            for i=1:CC2.NumObjects
                m=mask{4};p1=polyshape(m);
                p2=polyshape(CH2(i).ConvexHull);
                if intersect(p1,p2).NumRegions~=0
                    tmp=CC2.PixelIdxList{i};
                    I2(tmp)=1;
                end
            end        
            imgAND=I1&imtranslate(I2,[ofst(k),0]);
            CCAND=bwconncomp(imgAND);
            CCANDrp=regionprops(CCAND);
            
            ANDcoo=struct2cell(CCANDrp);
            areaList=ANDcoo(1,:)';
            BorderAreaCount(k)=sum(cell2mat(areaList));
            tmp=ANDcoo(2,:)';
            ANDcoo=[cell2mat(tmp) 100*ones(size(tmp,1),1) ];
            BorderCount(k)=size(ANDcoo,1);
            I=insertShape(I,'Circle',ANDcoo,'LineWidth',8,'Color','white');
            
            
            IB=I1;IB(:,:,2:3)=0;IB=circshift(IB,2,3);
            IG=I2;IG(:,:,2:3)=0;IG=circshift(IG,1,3);
            I(IB)=255; I(IG)=255;
        end
        clear K;% save memeory, ROI output 5
        %I=zeros(size(BW));
        %I=insertShape(I,'Line',vec(:,1:2),'LineWidth',8,'Color','white');
        %I=insertShape(I,'Circle',bvorigtf,'LineWidth',8,'Color','white');
        %doesn'twork %I=insertShape(I,'Polygon',{[Ts(:,1:2);Tm(:,1:2)];[Tms(:,1:2);Tm(:,1:2)];[Ts(:,1:2);Ta(:,1:2)];[Tas(:,1:2);Ta(:,1:2)]},'LineWidth',8,'Color','white');
        %bwmask=poly2mask(,,r,c)


        
        %I=imbinarize(I(:,:,1));
        
        %look for nearest AMA 
        
        % insert text 
        % need connection count, total width, from both side, and
        % percentages between
        [~,BorderCountMinIndex]=min(BorderCount);
        BorderCountPercent=min(BorderCount)/max(BorderCount);
        [~,BorderAreaCountMinIndex]=min(BorderAreaCount);
        BorderAreaCountPercent=min(BorderAreaCount)/max(BorderAreaCount);
        text_str=cell(6,1);
        text_str{1}=['LMA#:' num2str(BorderCount(1))];
        text_str{2}=['LMA#:' num2str(BorderCount(2))];
        text_str{3}=['LMA# ratio:' num2str(BorderCountPercent*100) '%'];
        text_str{4}=['sum LMA width:' num2str(BorderAreaCount(1))];
        text_str{5}=['sum LMA width:' num2str(BorderAreaCount(2))];
        text_str{6}=['LMA width ratio:' num2str(BorderAreaCountPercent*100) '%'];
        rightPos=c-1000;
        position=[0,0;c-1000,0;(BorderCountMinIndex-1)*(c-1000),100;0,200;c-1000,200;(BorderAreaCountMinIndex-1)*(c-1000),300];
        I=insertText(I,position,text_str,'FontSize',70,'BoxColor','white','BoxOpacity',0,'TextColor','white','AnchorPoint','LeftTop');
        
        imgName=strcat('../output4Connect/',imgList(j).name);
        imwrite(I,imgName); clear I;
        
        
        
        imgIdx=imgIdx+1;
    end
    cd .. % revert back to upper directories
end
%save statcell.mat statcell 

function BW=opengauss(imgName)
    t=Tiff(imgName,'r') ;%open image
    imageData = read(t); %read image
    imageData(4600:end,4000:end)=0; % remove label by hand
    imageData1=imageData(:,:,1);
    imageData1=imgaussfilt(imageData1,2); 
    %^^many benefits, filling holes of granny textures, remove fringe noises to
    %be below thresholds, and connect the fault lines of the images
    imageData1(imageData1<40)=0;% threshold
    %imshow(imageData1)
    BW=imbinarize(imageData1,'adaptive','Sensitivity',1);
    BW=bwareaopen(BW, 500); %remove isolated
end

function [bst1,bst2]=symm(BW,xcenter)         
         %xcenter=stats.Centroid(1);
         change=0; minsum=255*5000;  bst1=0;bst2=0;
         [r,c,d]=size(BW);yi=1:1:r;
         intercept=50; width=800; 
         xminS=xcenter-intercept-width/2; xmaxS=xcenter+intercept+width/2;
         xmin=xminS; %xmax=xmaxS;x1=xmin; x2=xmin;
         
         while xmin<xmaxS-width
             xmax=xmin+width;
             x1=xmin;
             
             while x1<=xmax
                x2=xmin+(xmax-x1);
                xi=ceil(linspace(x1,x2,yi(end)));
                %c=improfile(imageData,[x1 x2], [1 5000]);
                %figure("Visible",false);
                ind=sub2ind([r c],yi,xi);
                cursum=sum(BW(ind));

                if cursum<minsum %bingo
                    minsum=cursum;
                    bst1=x1; bst2=x2;
                    change=change+1;
                end
                %x2=x2+5;
                x1=x1+5;
             end     % new end
             
             xmin=xmin+5;
         end
end
function []=borderCC(tf, Vs,Vm,Vms,Va,Vas)
    %input
    %output
    Ts =(tf*Vs')';
    Tm =(tf*Vm')';
    Tms=(tf*Vms')';
    Ta =(tf*Va')';    
    Tas=(tf*Vas')';
    cell={Ts;Tm;Tms;Ta;Tas};
    %prompt
end

%show images in tight
% imds=imageDatastore('output5RegionPoly');
% montage(imds,'Interpolation','bilinear');
% exportgraphics(gca,'RegionPolyNonstroke.png');
% imds=imageDatastore('output4Connect');
% montage(imds,'Interpolation','bilinear');
% exportgraphics(gca,'ConnectStroke.png');
% imds=imageDatastore('output4Connect');
% mont=montage(imds,'Interpolation','bilinear','ThumbnailSize',[500 500]);
% montage_IM=mont.CData;
% imwrite(montage_IM,'ConnectNonStroke.png');%new write to file
opt_AOI= 'true';
statcell={};
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
        t=Tiff(imgList(j).name,'r') ;%open image
        imageData = read(t); %read image
        imageData(4600:end,4000:end)=0; % remove label by hand
        imageDataBinary=imageData(:,:,1)>100;
        
        % find best line to divide
        change=0; minsum=255*5000;  bst1=0;bst2=0;
        [r,c,d]=size(imageData);yi=1:1:r;


        %new start
         imageDataBinary = bwareaopen(imageDataBinary, 50); % remove small pixels groups
         CH=bwconvhull(imageDataBinary);
         stats = regionprops(CH);
         xcenter=stats.Centroid(1);
         intercept=50; width=700; 
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
                cursum=sum(imageData(ind));

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
         

        if opt_AOI=='true'
            figure
            imshow(imageData);
            
            roi=images.roi.Polygon;
           
            
            draw(roi) % interactive drawing polygon if area of interest
            
            line=images.roi.Line;
            line.Position=[bst1 0; bst2 r];
            line.Parent=roi.Parent;% force display without interaction

            coeff=polyfit([bst1;bst2],[0;r],1);
          % ---flip alone line y=mx+b
            roip=roi.Position;
            % translate so that y=mx, -Tb
            Tb=[0;coeff(2)];
            % rotate, so the line is on x axis Rt^-1
            theta=atan(coeff(1));
            Rt=[cos(theta) -sin(theta);sin(theta) cos(theta)];
            % flip along x axis
            Xr=[1 0;0 -1];
            % rotate back Rt: Rt\
            % translate back +Tb
            roip_=(Tb+Rt*Xr*(Rt\(roip'-Tb)))';
            roi_=images.roi.Polygon;
            roi_.Position=roip_;
            roi_.Parent=roi.Parent;
            %roi_.InteractionsAllowed='none';
            AOISum =sum(sum(sum(imageData.*uint8(createMask(roi,imageData)))));
            AOISum_=sum(sum(sum(imageData.*uint8(createMask(roi_,imageData)))));
            %
            bw=uint8(createMask(roi,imageData)|createMask(roi_,imageData));
            CCrgb=imageData.*repmat(bw,[1,1,3]);
            
            % Final image stat insertion
            imgName=strcat('../outputAOI/',imgList(j).name);
            text_str = cell(3,1);
            text_str{1}=['RedSum1: ' num2str(AOISum) ];
            text_str{2}=['RedSum2: ' num2str(AOISum_) ];
            AOI_P=min(AOISum,AOISum_)/max(AOISum,AOISum_);% percentage, lesser over greater
            [cx,cy]=centroid(polyshape(roip(:,1),roip(:,2)));
            text_str{3}=[num2str(AOI_P*100) '%'];
            position = [(cx>c/2)*(c-1000) 0;(cx<c/2)*(c-1000) 0;xor(cx>c/2,AOISum>AOISum_)*(c-1000) 100]; 
            box_color = 'white';
            anchorPoint='LeftTop';
            
            CCrgb=insertText(CCrgb,position,text_str,'FontSize',70,'BoxColor',...
                box_color,'BoxOpacity',0,'TextColor','white','AnchorPoint',anchorPoint);
            
            % Output
            imwrite(CCrgb,imgName); 
        end
        
        
        imgIdx=imgIdx+1;
        
    end
    
    cd .. % revert back to upper directories
end
%save statcell.mat statcell 

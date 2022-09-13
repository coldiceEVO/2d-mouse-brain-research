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

%         xmin=2000;xmax=3000;x1=xmin; x2=xmin;
%         while x1<=xmax
%             x2=xmin;
%             while x2<=xmax
%                 xi=ceil(linspace(x1,x2,yi(end)));
%                 %c=improfile(imageData,[x1 x2], [1 5000]);
%                 %figure("Visible",false);
%                 ind=sub2ind([r c],yi,xi);
%                 cursum=sum(imageData(ind));
% 
%                 if cursum<minsum %bingo
%                     minsum=cursum;
%                     bst1=x1; bst2=x2;
%                     change=change+1;
%                 end
%                 x2=x2+5;
%             end
%             x1=x1+5;
%         end     
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
         
        close 
        figure
            f=imshow(imageData);

            % Left Right boundary of each line
            pl=line([bst1 bst2], [1 c]);
            pl.LineWidth=1.5;
            % disp([get(gcf,'Number') bst1 bst2]);
            imgName=imgList(j).name;
            imgName=strcat('../output/',imgName);
            % exportgraphics()
            saveas(f,imgName)
            
        % section stats L/R pixel count, ratio, 
        % a line cut first, revert potential 
        imageDataBinary=insertShape(imageData,'line',[bst1,0,bst2,r],'color','black','LineWidth',2) ;%sever it with dark line to prevent connection from left and right
        % (1) static luminance, 
        %         imageDataBinary=im2bw(imageDataBinary,0.05);% new binary image with lower threshold, value is lu
        % or
        %         imageDataBinary=imageDataBinary(:,:,1)>50;
        
        % (2) graythresh function for each image, adjust by multiplier if needed
                imageDataBinary=imbinarize(imageDataBinary,graythresh(imageData)*0.7);% new binary image with lower threshold, value is lu
                imageDataBinary=imageDataBinary(:,:,1);
        
        imageDataBinary = bwareaopen(imageDataBinary, 500) ; %reduce small ones
        
        % show effect of filtering
        % imshowpair(imageData,imbinarize(imageData,graythresh(imageData)),'montage')
        % imshowpair(imageData,im2bw(imageData,graythresh(imageData)),'montage')
        
        LMask=poly2mask([0 bst1 bst2 0],[0 0 r r],r,c); % vertical, horizontal
        RMask=poly2mask([c bst1 bst2 c],[0 0 r r],r,c);
        LSum=sum(sum(imageDataBinary.*LMask));
        RSum=sum(sum(imageDataBinary.*RMask));
        PercentLesser=min(LSum,RSum)/max(LSum,RSum);% persentage of lesser compare to greater side
        

        
        CC = bwconncomp(imageDataBinary) ;%connected component
        L=labelmatrix(CC);
        imgName=strcat('../output1/',imgList(j).name);
        CCrgb=label2rgb(L,'hsv','k','shuffle');
        imwrite(CCrgb,imgName); % arg2: color map
      % infograph
        imgName=strcat('../output2/',imgList(j).name);
        % line
        CCrgb=insertShape(CCrgb,'line',[bst1,0,bst2,r],'color','red','LineWidth',20);

        
        stats=regionprops(CC) ;% connect component's stats: Area/Centroid/BoundingBox
        CCCentCord=cat(1,stats.Centroid) ;% separate L/R based on centroid's position to dividing line
        listR=inpolygon(CCCentCord(:,1),CCCentCord(:,2), [c bst1 bst2 c],[0 0 r r]); % 0=L 1=R
        idxR=find(listR);
        idxL=find(~listR);
%         %show image
%         BW_R=L;
%         BW_R(ismember(L,idxR))=1;
%         BW_R(ismember(L,idxL))=2;
%         %imshow(label2rgb(BW_R,'jet','k','shuffle'));
%         CCArea=cat(1,stats.Area);
%         allAreaLeft=CCArea(idxL);
%         avgAreaLeft=mean(allAreaLeft);
%         allAreaRight=CCArea(idxR);
%         avgAreaRight=mean(allAreaRight);
%         %CC.PixelIdxList(idxR)
        
        
        numPixels = cellfun(@numel,CC.PixelIdxList);
        
        statcell(imgIdx,1:8)={imgList(j).name [bst1 bst2] [LSum RSum PercentLesser] CC numPixels {idxL idxR} [size(idxL,1) size(idxR,1)] [avgAreaLeft avgAreaRight]};
        
      % stat text
        
        text_str = cell(6,1);
        text_str{1}=['pixel count: ' num2str(LSum) ];
        text_str{2}=['CC count: ' num2str(size(idxL,1))];
        text_str{3}=['pixel count: ' num2str(RSum) ];
        text_str{4}=['CC count: ' num2str(size(idxR,1))];
        text_str{5}=['x1=' num2str(bst1)];
        text_str{6}=['x2=' num2str(bst2)];
        text_str{7}=[num2str(PercentLesser*100) '%'];
        position = [0 0;0 100;c-1000 0;c-1000 100;bst1 0;bst2 r-100;(LSum>RSum)*(c-1000) 200]; 
        box_color = 'white';
        anchorPoint={'LeftTop','RightTop','CenterTop','CenterBottom'};

        CCrgb=insertText(CCrgb,position,text_str,'FontSize',70,'BoxColor',...
            box_color,'BoxOpacity',0,'TextColor','white','AnchorPoint',anchorPoint{1});
        imwrite(CCrgb,imgName); % arg2: color map
        
        imgIdx=imgIdx+1;
        
        
    end
    
    cd .. % revert back to upper directories
end
save statcell.mat statcell 
% stat=cat(1,statcell)
% stat(3,2)
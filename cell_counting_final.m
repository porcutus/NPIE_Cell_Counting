%run this line before running others 
% close all;clear;clc;
% score_arr_total = zeros(11, 1);
% free_cells_total = 0

close all;clearvars -except score_arr_total free_cells_total;clc;
%--------read image(load all images from directory and montage)

% img=imread('/Users/Porcutus/Desktop/104-1/念祖生醫實驗室/Cell_Counting/922  matlab images/best/img2.tif');
% num_of_pics = 1;

D = dir('*.tif');
num_of_pics = numel(D);
imcell = cell(1,numel(D));
for i = 1:numel(D)
  imcell{i} = D(i).name;
end
handle = montage(imcell);
my_montage = getframe(gca);
img = my_montage.cdata;
%--------split channel
red = img(:,:,1); % Red channel
green = img(:,:,2); % Green channel
blue = img(:,:,3); % Blue channel
%--------subtract channels
onlycells = green - blue;
%figure(4); imshow(onlycells); title('onlycells');
bwCells = im2bw(onlycells, 0.4);
bwOnlyCells = bwareaopen(bwCells, 110);%2nd parameter: cell size
%figure(5); imshow(bwCells); title('Cells in binary');
%---------imfill
BWdfill = imfill(bwOnlyCells, 'holes');
%figure(6), imshow(BWdfill); title('binary image with filled holes');
%---watershed
aI=BWdfill(:,:,1);
D = -bwdist(~aI);
L=watershed(D);
bw2 = aI;
bw2(L == 0) = 0;
mask = imextendedmin(D,1);

D2=imimposemin(D,mask);
Ld2=watershed(D2);
bw3 = aI;
bw3(Ld2 == 0) = 0;
%---------calculate number of cells
[cell_map,num]=bwlabel(bw3);
cell_num=max(num(:))%why max?
%---------get centroid
s = regionprops(bw3, 'centroid');
ctrd_f = cat(1, s.Centroid);
x_f = ctrd_f(:,1);
y_f = ctrd_f(:,2);
% figure(9); imshow(bw3); title('Centroid');
% hold on
% plot(x_f, y_f, 'b*')
% hold off

%-------------------------------pluri-beads--------------------------------

pout_imadjust = imadjust(blue, [0.1, 0.2],[]);%調整第二個參數可以辨識...還內環外？
img_filter = im2bw(pout_imadjust, 0.9);%第二個參數（level）可以修正破掉的環
% se = strel('ball',5,5);
% img_filter = imdilate(blue, se);
%---------remove objects by the border
img_clearborder1 = imclearborder(img_filter, 8);
%---------get area
a = regionprops(img_clearborder1, 'Area');
Ars = cat(1, a.Area);
%---------remove small objects
switch num_of_pics
    case 1
        pluri = bwareaopen(img_clearborder1, 250);
        test = 1;
    case 4
        pluri = bwareaopen(img_clearborder1, 110);
        test = 4;
    case 9
        pluri = bwareaopen(img_clearborder1, 55);
        test = 9;
    otherwise
        pluri = bwareaopen(img_clearborder1, 55);
        test = 'otws';
end
%pluri = bwareaopen(img_clearborder1, 250); %for 1 pic
%pluri = bwareaopen(img_clearborder1, 110); %for 4 pics
%pluri = bwareaopen(img_clearborder1, 55); %for 9 pics
%figure; imshow(pluri);
%---------find radii
%[centers, radii] = imfindcircles(pluri,[7 15],'Sensitivity',0.9);
%因為有些取pluri內圈 有些取pluri整圈
[centers, radii] = imfindcircles(pluri,[15 25],'Sensitivity',0.9);
%[centers, radii] = imfindcircles(pluri,[20 50],'Sensitivity',0.9);
% hold on
% h = viscircles(centers,radii);
% %plot(centers(:,1),centers(:,2), 'b*');
% hold off;
x_p = centers(:,1);
y_p = centers(:,2);
%--------pluri_num
pluri_num = numel(radii);
%--------visualize labels
% imshow(img, 'InitialMag', 'fit')
% hold on
% plot(centers(:,1), centers(:,2), 'b*');
% for k = 1:numel(radii)
%    t = text(centers(k), centers(k+pluri_num), sprintf('%d', k));
% end
% hold off

%-------------------------------distance--------------------------------

%---combine coordinates
all_cells_xy = squeeze(num2cell(permute(cat(3,x_f,y_f),[3,1,2]),1));
%---count how many cells attach to each bead
cntarr = zeros(numel(x_p),1);
for i = 1:numel(x_p)
%for i = 1:1
    %---count distance
    dist = zeros(numel(all_cells_xy),1);
    a = 1;
    
    for j = 1:numel(all_cells_xy)
        dist(a) = sqrt((x_p(i)-all_cells_xy{j}(1))^2 + (y_p(i)-all_cells_xy{j}(2))^2);
        a = a + 1;
    end;
    %---map distance and cell coordinates together
    distance_map = containers.Map(dist, all_cells_xy);
    %---sort distance from near to far
%     sortarr = sortrows(dist);
    %---count cells near beads
%     for x = 1:numel(sortarr)
%         %if sortarr(x)>45 %montage 1
%         if sortarr(x)>22.5 %montage 4
%         %if sortarr(x)>17 %montage 6
%             break;
%         end
%     end 
    switch num_of_pics
        case 1
            out_of_radius = num2cell(dist(dist > 45));
        case 4
            out_of_radius = num2cell(dist(dist > 22.5));
        case 6 
            out_of_radius = num2cell(dist(dist > 17));
        otherwise
            out_of_radius = num2cell(dist(dist > 17));
    end
    remove(distance_map, out_of_radius);
    catched = length(distance_map);
    cntarr(i) = catched;
    
    used_cells_xy = values(distance_map);

    remove_arr = zeros(numel(used_cells_xy),1);
    for k = 1:numel(used_cells_xy)
        for l = 1:numel(all_cells_xy)
            if(isequal(all_cells_xy(l),used_cells_xy(k)))
                remove_arr(k) = l;
                break;
            else
                continue;
            end
        end
    end
    if(remove_arr~=0)
        all_cells_xy(remove_arr) = [];
    end;
    
    if(numel(all_cells_xy)==0)
        break;
    end;
end;

free_cells = numel(all_cells_xy);
free_cells_total = free_cells_total + free_cells;
%---plot coordinates
figure; imshow(img);
hold on
plot(x_f, y_f, 'b*')
for k = 1:numel(radii)
   t = text(centers(k), centers(k+numel(radii)), sprintf('%d', k));
end
hold off
%d = imdistline;

%-------------------------------draw-graph---------------------------------

score_arr = zeros(11, 1);
for i = 1:numel(cntarr)
   score_arr(cntarr(i)+1) = score_arr(cntarr(i)+1) + 1;
   score_arr_total(cntarr(i)+1) = score_arr_total(cntarr(i)+1) + 1;
end
num_arr = [0 1 2 3 4 5 6 7 8 9 10];
figure; bar(num_arr, score_arr);
% for i = numel(score_arr)
%     text(num_arr(i), score_arr(i), num2str(score_arr(i),'%0.2f'),... 
%     'HorizontalAlignment','center',... 
%     'VerticalAlignment','top');
% end
sum = 0;
for i = 1:numel(num_arr)
    sum = sum + num_arr(i)*score_arr(i);
end
avg = sum/pluri_num;
dev = std(cntarr);
  
display(pluri_num);
display(cell_num);
display(avg)
display(dev);
display(free_cells);

%Z1 = all_cells_xy(find(all_cells_xy~=test))


 

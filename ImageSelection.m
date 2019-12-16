clear all

annotations = readtable('D:\Documents\GitHub\idr-metadata\idr0003-breker-plasticity\screenA\idr0003-screenA-annotation.csv');

%Disregard columns for Phenotype 12, 13, 14 in the future, but not now

annotations.Phenotype1 = categorical(annotations.Phenotype1);
annotations.Phenotype2 = categorical(annotations.Phenotype2);
annotations.Phenotype3 = categorical(annotations.Phenotype3);
annotations.Phenotype4 = categorical(annotations.Phenotype4);
annotations.Phenotype5 = categorical(annotations.Phenotype5);
annotations.Phenotype6 = categorical(annotations.Phenotype6);
annotations.Phenotype7 = categorical(annotations.Phenotype7);
annotations.Phenotype8 = categorical(annotations.Phenotype8);
annotations.Phenotype9 = categorical(annotations.Phenotype9);
annotations.Phenotype10 = categorical(annotations.Phenotype10);
annotations.Phenotype11 = categorical(annotations.Phenotype11);
annotations.Phenotype12 = categorical(annotations.Phenotype12);
annotations.Phenotype13 = categorical(annotations.Phenotype13);
annotations.Phenotype14 = categorical(annotations.Phenotype14);

annotations.ControlGFPLocalization = categorical(annotations.ControlGFPLocalization);

summary(annotations)


%Create trimedAnnotations which contains wells with one and only one
%phenotype (excluding Phenotypes 12, 13, 14) and that were not labeled
%'below threshold' in ControlGFPLocalization

newLabelsAnnotations = annotations(:, {'Phenotype1', 'Phenotype2', 'Phenotype3', 'Phenotype4', ...
    'Phenotype5', 'Phenotype6', 'Phenotype7', 'Phenotype8', 'Phenotype9', 'Phenotype10', ...
    'Phenotype11', 'Phenotype12', 'Phenotype13', 'Phenotype14'});
sumHolder = summary(newLabelsAnnotations(:, 1));
tempLabelsAnnotations(:, 1) = newLabelsAnnotations.Phenotype1 == sumHolder.Phenotype1.Categories;
sumHolder = summary(newLabelsAnnotations(:, 2));
tempLabelsAnnotations(:, 2) = newLabelsAnnotations.Phenotype2 == sumHolder.Phenotype2.Categories;
sumHolder = summary(newLabelsAnnotations(:, 3));
tempLabelsAnnotations(:, 3) = newLabelsAnnotations.Phenotype3 == sumHolder.Phenotype3.Categories;
sumHolder = summary(newLabelsAnnotations(:, 4));
tempLabelsAnnotations(:, 4) = newLabelsAnnotations.Phenotype4 == sumHolder.Phenotype4.Categories;
sumHolder = summary(newLabelsAnnotations(:, 5));
tempLabelsAnnotations(:, 5) = newLabelsAnnotations.Phenotype5 == sumHolder.Phenotype5.Categories;
sumHolder = summary(newLabelsAnnotations(:, 6));
tempLabelsAnnotations(:, 6) = newLabelsAnnotations.Phenotype6 == sumHolder.Phenotype6.Categories;
sumHolder = summary(newLabelsAnnotations(:, 7));
tempLabelsAnnotations(:, 7) = newLabelsAnnotations.Phenotype7 == sumHolder.Phenotype7.Categories;
sumHolder = summary(newLabelsAnnotations(:, 8));
tempLabelsAnnotations(:, 8) = newLabelsAnnotations.Phenotype8 == sumHolder.Phenotype8.Categories;
sumHolder = summary(newLabelsAnnotations(:, 9));
tempLabelsAnnotations(:, 9) = newLabelsAnnotations.Phenotype9 == sumHolder.Phenotype9.Categories;
sumHolder = summary(newLabelsAnnotations(:, 10));
tempLabelsAnnotations(:, 10) = newLabelsAnnotations.Phenotype10 == sumHolder.Phenotype10.Categories;
sumHolder = summary(newLabelsAnnotations(:, 11));
tempLabelsAnnotations(:, 11) = newLabelsAnnotations.Phenotype11 == sumHolder.Phenotype11.Categories;
sumHolder = summary(newLabelsAnnotations(:, 12));
tempLabelsAnnotations(:, 12) = newLabelsAnnotations.Phenotype12 == sumHolder.Phenotype12.Categories;
sumHolder = summary(newLabelsAnnotations(:, 13));
tempLabelsAnnotations(:, 13) = newLabelsAnnotations.Phenotype13 == sumHolder.Phenotype13.Categories;
sumHolder = summary(newLabelsAnnotations(:, 14));
tempLabelsAnnotations(:, 14) = newLabelsAnnotations.Phenotype14 == sumHolder.Phenotype14.Categories;

belowThreshold = annotations.ControlGFPLocalization == 'below threshold';

forCullingRef = int8(tempLabelsAnnotations);
toCull = false(length(forCullingRef(:, 1)), 1);
cullingTotals = sum(forCullingRef, 2);
for i = 1:length(cullingTotals)
    if cullingTotals(i) ~= 1 || forCullingRef(i, 12) == 1 || forCullingRef(i, 13) == 1 || forCullingRef(i, 14) == 1 || belowThreshold(i) == 1
        toCull(i) = true;
    end
end
toKeep = logical(ones(length(toCull), 1) - toCull);

trimedAnnotations = annotations(toKeep, :);


%Optimize the data structure
sampleList = trimedAnnotations(:, {'Plate', 'Well'});
tempForCullingRef = forCullingRef(toKeep, :);
for i = 1:length(tempForCullingRef(1, :))
    tempForCullingRef(:, i) = tempForCullingRef(:, i) * i;
end
newColumn = sum(tempForCullingRef, 2);
sampleList = addvars(sampleList, newColumn, 'NewVariableNames', 'Localization');

save('LocalizationReference_v1.mat', 'sampleList');




api = 'https://idr.openmicroscopy.org';


%Create a conversion chart for plate name + well grid position to image ID.
url3 = [api '/webclient/api/plates/?id=51'];
data3 = webread(url3);
[plateRegistry(1:size(data3.plates)).name] = data3.plates.name;
[plateRegistry(1:size(data3.plates)).id] = data3.plates.id;
imageRegistry = struct('plate', {}, 'well', {}, 'imageid', {});
for i = 1:length(plateRegistry)
    url2 = [api '/webgateway/plate/' num2str(plateRegistry(i).id) '/'];
    data2 = webread(url2);
    for j = 1:16
        for k = 1:24
            locus = ((i - 1) * 384) + ((j - 1) * 24) + k;
            imageRegistry(locus).plate = plateRegistry(i).name;
            imageRegistry(locus).well = [char(data2.rowlabels(j)) num2str(data2.collabels(k))];
            imageRegistry(locus).imageid = data2.grid(j, k).id;
        end
    end
end

save('ImageRegistry_v1.mat', 'imageRegistry');



%Scan all images for usable signal
% options = weboptions('Timeout', 10);
% outcome = zeros(length(sampleList.Localization), 1); % 0 = untried, 1 = timeout, 2 = image inadequate, 3 = image adequate
% for j = 1:length(sampleList.Localization)
%     downSelect = imageRegistry(strcmp({imageRegistry.plate}, char(ourSamples.Plate(j))));
%     well = downSelect(strcmp({downSelect.well}, char(ourSamples.Well(j))));
%     url = [api '/webgateway/render_image/' num2str(well.imageid) '/0/0/'];
%     try
%         img = webread(url, options);
%         mCherry = img(:, :, 1);
%         GFP = img(:, :, 2);
% 
%         %perform initial check for dynamic range and signal
%         %resolution in both channels
%         dim = size(img);
%         pixel_num = dim(1) * dim(2);
%         mins_mCherry = zeros(round(pixel_num / 10), 1) - ones(round(pixel_num / 10), 1);
%         mins_GFP = zeros(round(pixel_num / 10), 1) - ones(round(pixel_num / 10), 1);
%         max_mCherry = max(max(mCherry));
%         max_GFP = max(max(GFP));
%         temp_mCherry = mCherry;
%         temp_GFP = GFP;
%         for k = 1:round(pixel_num / 10)
%             [temp_mins, rows] = min(temp_mCherry);
%             [mins_mCherry(k), column] = min(temp_mins);
%             temp_mCherry(rows(column), column) = max_mCherry;
%             [temp_mins, rows] = min(temp_GFP);
%             [mins_GFP(k), column] = min(temp_mins);
%             temp_GFP(rows(column), column) = max_GFP;
%         end
%         floor_mCherry = mean(mins_mCherry);
%         floor_GFP = mean(mins_GFP);
%         above_ceiling_mCherry = sum(sum(mCherry >= (5 * floor_mCherry)));
%         above_ceiling_GFP = sum(sum(GFP >= (5 * floor_GFP)));
%         res_mCherry = sum(sum(mCherry >= 50));
%         res_GFP = sum(sum(GFP >= 50));
% 
%         if res_mCherry >= 100 && res_GFP >= 100 && (above_ceiling_mCherry / pixel_num) >= 0.001 && (above_ceiling_GFP / pixel_num) >= 0.001
%             outcome(j) = 3;
%         else
%             outcome(j) = 2;
%         end
%     catch
%         outcome(j) = 1;
%     end
%     j
% end



%Create a set of cells with composition of choice
%Define ratio of localization types
locRatio = [1, 1, 10, 2, 3, 1, 1, 6, 3, 1, 1];
% locRatio = [0, 1, 0, 2, 0, 0, 0, 0, 0, 0, 1];

%Set minimum number of cells per type
minCells = 2000;

%Set minimum number of images per type
minImages = 20;

%Adjust search constraints
minCellsList = locRatio * minCells;

options = weboptions('Timeout', 10);
firstTime = true;
actualCellCounts = zeros(1, 11);
actualImageCounts = zeros(1, 11);
for i = 1:11
    if locRatio(i) > 0
        ourSamples = sampleList(sampleList.Localization == i, :);
        while actualImageCounts(i) < minImages && ~isempty(ourSamples)
            selection = randi(length(ourSamples.Localization));
            downSelect = imageRegistry(strcmp({imageRegistry.plate}, char(ourSamples.Plate(selection))));
            well = downSelect(strcmp({downSelect.well}, char(ourSamples.Well(selection))));
            ourSamples = [ourSamples(1:(selection - 1), :); ourSamples((selection + 1):end, :)];
            url = [api '/webgateway/render_image/' num2str(well.imageid) '/0/0/'];
            try
                img = webread(url, options);
                mCherry = img(:, :, 1);
                GFP = img(:, :, 2);

                %first perform initial check for dynamic range and signal
                %resolution in both channels
                dim = size(img);
                pixel_num = dim(1) * dim(2);
                mins_mCherry = zeros(round(pixel_num / 10), 1) - ones(round(pixel_num / 10), 1);
                mins_GFP = zeros(round(pixel_num / 10), 1) - ones(round(pixel_num / 10), 1);
                max_mCherry = max(max(mCherry));
                max_GFP = max(max(GFP));
                temp_mCherry = mCherry;
                temp_GFP = GFP;
                for k = 1:round(pixel_num / 10)
                    [temp_mins, rows] = min(temp_mCherry);
                    [mins_mCherry(k), column] = min(temp_mins);
                    temp_mCherry(rows(column), column) = max_mCherry;
                    [temp_mins, rows] = min(temp_GFP);
                    [mins_GFP(k), column] = min(temp_mins);
                    temp_GFP(rows(column), column) = max_GFP;
                end
                floor_mCherry = mean(mins_mCherry);
                floor_GFP = mean(mins_GFP);
                above_ceiling_mCherry = sum(sum(mCherry >= (5 * floor_mCherry)));
                above_ceiling_GFP = sum(sum(GFP >= (5 * floor_GFP)));
                res_mCherry = sum(sum(mCherry >= 50));
                res_GFP = sum(sum(GFP >= 50));

                if res_mCherry >= 100 && res_GFP >= 100 && (above_ceiling_mCherry / pixel_num) >= 0.001 && (above_ceiling_GFP / pixel_num) >= 0.001

                    filter = fspecial('average', [30 30]);
                    temp = imfilter(mCherry, filter);
                    temp2 = mCherry - temp;
                    mask = temp2 > 3*mean(temp2, 'all');

                    %Temp
%                     test = uint8(mask2);
%                     test(mask2 == 1) = mCherry(mask2 == 1);
%                     figure()
%                     imshow(test)

        %             lmask = bwconncomp(mask);
                    cells = regionprops(mask, 'Area', 'BoundingBox', 'PixelIdxList', 'Circularity', 'Eccentricity');
                    area = [cells.Area];
                    cells = cells(area >= 75 & area <= 576);
                    dim1 = zeros(length(cells), 1);
                    dim2 = zeros(length(cells), 1);
                    for k = 1:length(cells)
                        dim1(k) = cells(k).BoundingBox(3);
                        dim2(k) = cells(k).BoundingBox(4);
                    end
                    cells = cells(dim1 <= 35 & dim2 <= 35);
                    circularity = [cells.Circularity];
                    cells = cells(circularity >= 0.75);
                    eccentricity = [cells.Eccentricity];
                    cells = cells(eccentricity <= 0.8);


                    %Make mask out of size-gated cells
                    cellPixels = cells(1).PixelIdxList;
                    for k = 2:length(cells)
                        cellPixels = [cellPixels; cells(k).PixelIdxList];
                    end
                    mask2 = false(512, 672);
                    mask2(cellPixels) = 1;


                    filter2 = fspecial('average', [50 50]);
                    temp3 = imfilter(GFP, filter);
                    new_GFP = GFP - temp3;

                    mask3 = new_GFP > 3*mean(new_GFP, 'all');
%                     figure();
%                     imshow(mask3)
%                     mask_combo = mask2 & mask3;
%                     figure()
%                     imshow(mask_combo)
%                     mask_remain = mask2 & ~mask_combo;
%                     figure()
%                     imshow(mask_remain)

                    mask_GFPtop = mask3 & ~mask2;
%                     figure()
%                     imshow(mask_GFPtop)
                    tells = regionprops(mask_GFPtop, 'Area', 'PixelIdxList');
                    tells_area = [tells.Area];
                    tells = tells(tells_area <= (2 * mean(tells_area)) & tells_area > 2);
                    tellPixels = tells(1).PixelIdxList;
                    for k = 2:length(tells)
                        tellPixels = [tellPixels; tells(k).PixelIdxList];
                    end
                    mask_tells = false(512, 672);
                    mask_tells(tellPixels) = 1;
%                     figure();
%                     imshow(mask_tells)

                    adjust = true;
                    while adjust == true
                        adjust = false;
                        directions = [1, 0, 0; -1, 0, 0; 0, 1, 0; 0, -1, 0];
                        tells_total = sum(sum(mask_tells));
                        tells_objects = zeros(4, 1);
                        for k = 1:4
                            new_mask_GFPtop = imtranslate(mask3, directions(k, 1:2)) & ~mask2;
                            new_tells = regionprops(new_mask_GFPtop, 'Area', 'PixelIdxList');
                            new_tells_area = [new_tells.Area];
                            new_tells = new_tells(new_tells_area <= (2 * mean(new_tells_area)) & new_tells_area > 2);
                            tells_objects(k) = length(new_tells);
                            new_tellPixels = new_tells(1).PixelIdxList;
                            for m = 2:length(new_tells)
                                new_tellPixels = [new_tellPixels; new_tells(m).PixelIdxList];
                            end
                            new_mask_tells = false(512, 672);
                            new_mask_tells(new_tellPixels) = 1;
                            new_tells_total = sum(sum(new_mask_tells));
                            if new_tells_total <= tells_total
                                directions(k, 3) = tells_total - new_tells_total;
                            end
                        end
                        if max(directions(:, 3)) >= (2.5 * mean(tells_objects))
                            adjust = true;
                            [~, I] = max(directions(:, 3));

                            new_mask_GFPtop = imtranslate(mask3, directions(I, 1:2)) & ~mask2;
                            new_tells = regionprops(new_mask_GFPtop, 'Area', 'PixelIdxList');
                            new_tells_area = [new_tells.Area];
                            new_tells = new_tells(new_tells_area <= (2 * mean(new_tells_area)) & new_tells_area > 2);
                            new_tellPixels = new_tells(1).PixelIdxList;
                            for m = 2:length(new_tells)
                                new_tellPixels = [new_tellPixels; new_tells(m).PixelIdxList];
                            end
                            new_mask_tells = false(512, 672);
                            new_mask_tells(new_tellPixels) = 1;
                            mask_tells = new_mask_tells;

                            mask3 = imtranslate(mask3, directions(I, 1:2));
                            new_GFP = imtranslate(new_GFP, directions(I, 1:2));
%                             figure();
%                             imshow(mask_tells);
                        end
                    end


                    %For grabbing images
%                     figure();
%                     imshow(new_GFP);
%                     figure();
%                     imshow(mask2);


                    cells = regionprops(mask2, new_GFP, 'Area', 'BoundingBox', 'Image', 'PixelList', 'PixelIdxList', 'PixelValues', 'Centroid', 'WeightedCentroid', 'MeanIntensity');
                    cellsin = [cells.MeanIntensity];
                    inmean = mean(cellsin);
                    instd = std(cellsin);
                    cells = cells(cellsin >= (inmean - (2 * instd)) & cellsin <= (inmean + (2 * instd)));
                    cellsroster = [zeros(length(cells), 1) + i, zeros(length(cells), 1) + well.imageid];
                    actualCellCounts(i) = actualCellCounts(i) + length(cells);
                    actualImageCounts(i) = actualImageCounts(i) + 1;
                    if firstTime
                        firstTime = false;
                        cellMasterData = cells;
                        cellMasterRoster = cellsroster;
                    else
                        cellMasterData(end + 1 : end + length(cells)) = cells;
                        cellMasterRoster = [cellMasterRoster; cellsroster];
                    end

                    %From presentation code - to explore
%                     newGFPMatrix = zeros(length(cells), 625) - ones(length(cells), 625);
%                     for k = 1:length(cells)
%                         for m = 1:length(cells(k).PixelIdxList)
%                             newGFPMatrix(k, m) = new_GFP(cells(k).PixelIdxList(m));
%                         end
%                     end
%                     GFPMatrix((sum(cellCounts) + 1):(sum(cellCounts) + length(cells)), :) = newGFPMatrix;
%                     cellLabels((sum(cellCounts) + 1):(sum(cellCounts) + length(cells))) = i;
%                     cellCounts(i) = cellCounts(i) + length(cells);
                end
            catch
                disp("Error occured, trying again...");
            end
        end
        if actualImageCounts(i) == minImages
            while actualCellCounts(i) < minCellsList(i) && ~isempty(ourSamples)
                selection = randi(length(ourSamples.Localization));
                downSelect = imageRegistry(strcmp({imageRegistry.plate}, char(ourSamples.Plate(selection))));
                well = downSelect(strcmp({downSelect.well}, char(ourSamples.Well(selection))));
                ourSamples = [ourSamples(1:(selection - 1), :); ourSamples((selection + 1):end, :)];
                url = [api '/webgateway/render_image/' num2str(well.imageid) '/0/0/'];
                try
                    img = webread(url, options);
                    mCherry = img(:, :, 1);
                    GFP = img(:, :, 2);

                    %first perform initial check for dynamic range and signal
                    %resolution in both channels
                    dim = size(img);
                    pixel_num = dim(1) * dim(2);
                    mins_mCherry = zeros(round(pixel_num / 10), 1) - ones(round(pixel_num / 10), 1);
                    mins_GFP = zeros(round(pixel_num / 10), 1) - ones(round(pixel_num / 10), 1);
                    max_mCherry = max(max(mCherry));
                    max_GFP = max(max(GFP));
                    temp_mCherry = mCherry;
                    temp_GFP = GFP;
                    for k = 1:round(pixel_num / 10)
                        [temp_mins, rows] = min(temp_mCherry);
                        [mins_mCherry(k), column] = min(temp_mins);
                        temp_mCherry(rows(column), column) = max_mCherry;
                        [temp_mins, rows] = min(temp_GFP);
                        [mins_GFP(k), column] = min(temp_mins);
                        temp_GFP(rows(column), column) = max_GFP;
                    end
                    floor_mCherry = mean(mins_mCherry);
                    floor_GFP = mean(mins_GFP);
                    above_ceiling_mCherry = sum(sum(mCherry >= (5 * floor_mCherry)));
                    above_ceiling_GFP = sum(sum(GFP >= (5 * floor_GFP)));
                    res_mCherry = sum(sum(mCherry >= 50));
                    res_GFP = sum(sum(GFP >= 50));

                    if res_mCherry >= 100 && res_GFP >= 100 && (above_ceiling_mCherry / pixel_num) >= 0.001 && (above_ceiling_GFP / pixel_num) >= 0.001

                        filter = fspecial('average', [30 30]);
                        temp = imfilter(mCherry, filter);
                        temp2 = mCherry - temp;
                        mask = temp2 > 3*mean(temp2, 'all');

                        %Temp
    %                     test = uint8(mask2);
    %                     test(mask2 == 1) = mCherry(mask2 == 1);
    %                     figure()
    %                     imshow(test)

            %             lmask = bwconncomp(mask);
                        cells = regionprops(mask, 'Area', 'BoundingBox', 'PixelIdxList', 'Circularity', 'Eccentricity');
                        area = [cells.Area];
                        cells = cells(area >= 75 & area <= 576);
                        dim1 = zeros(length(cells), 1);
                        dim2 = zeros(length(cells), 1);
                        for k = 1:length(cells)
                            dim1(k) = cells(k).BoundingBox(3);
                            dim2(k) = cells(k).BoundingBox(4);
                        end
                        cells = cells(dim1 <= 35 & dim2 <= 35);
                        circularity = [cells.Circularity];
                        cells = cells(circularity >= 0.75);
                        eccentricity = [cells.Eccentricity];
                        cells = cells(eccentricity <= 0.8);


                        %Make mask out of size-gated cells
                        cellPixels = cells(1).PixelIdxList;
                        for k = 2:length(cells)
                            cellPixels = [cellPixels; cells(k).PixelIdxList];
                        end
                        mask2 = false(512, 672);
                        mask2(cellPixels) = 1;


                        filter2 = fspecial('average', [50 50]);
                        temp3 = imfilter(GFP, filter);
                        new_GFP = GFP - temp3;

                        mask3 = new_GFP > 3*mean(new_GFP, 'all');
    %                     figure();
    %                     imshow(mask3)
    %                     mask_combo = mask2 & mask3;
    %                     figure()
    %                     imshow(mask_combo)
    %                     mask_remain = mask2 & ~mask_combo;
    %                     figure()
    %                     imshow(mask_remain)

                        mask_GFPtop = mask3 & ~mask2;
    %                     figure()
    %                     imshow(mask_GFPtop)
                        tells = regionprops(mask_GFPtop, 'Area', 'PixelIdxList');
                        tells_area = [tells.Area];
                        tells = tells(tells_area <= (2 * mean(tells_area)) & tells_area > 2);
                        tellPixels = tells(1).PixelIdxList;
                        for k = 2:length(tells)
                            tellPixels = [tellPixels; tells(k).PixelIdxList];
                        end
                        mask_tells = false(512, 672);
                        mask_tells(tellPixels) = 1;
    %                     figure();
    %                     imshow(mask_tells)

                        adjust = true;
                        while adjust == true
                            adjust = false;
                            directions = [1, 0, 0; -1, 0, 0; 0, 1, 0; 0, -1, 0];
                            tells_total = sum(sum(mask_tells));
                            tells_objects = zeros(4, 1);
                            for k = 1:4
                                new_mask_GFPtop = imtranslate(mask3, directions(k, 1:2)) & ~mask2;
                                new_tells = regionprops(new_mask_GFPtop, 'Area', 'PixelIdxList');
                                new_tells_area = [new_tells.Area];
                                new_tells = new_tells(new_tells_area <= (2 * mean(new_tells_area)) & new_tells_area > 2);
                                tells_objects(k) = length(new_tells);
                                new_tellPixels = new_tells(1).PixelIdxList;
                                for m = 2:length(new_tells)
                                    new_tellPixels = [new_tellPixels; new_tells(m).PixelIdxList];
                                end
                                new_mask_tells = false(512, 672);
                                new_mask_tells(new_tellPixels) = 1;
                                new_tells_total = sum(sum(new_mask_tells));
                                if new_tells_total <= tells_total
                                    directions(k, 3) = tells_total - new_tells_total;
                                end
                            end
                            if max(directions(:, 3)) >= (2.5 * mean(tells_objects))
                                adjust = true;
                                [~, I] = max(directions(:, 3));

                                new_mask_GFPtop = imtranslate(mask3, directions(I, 1:2)) & ~mask2;
                                new_tells = regionprops(new_mask_GFPtop, 'Area', 'PixelIdxList');
                                new_tells_area = [new_tells.Area];
                                new_tells = new_tells(new_tells_area <= (2 * mean(new_tells_area)) & new_tells_area > 2);
                                new_tellPixels = new_tells(1).PixelIdxList;
                                for m = 2:length(new_tells)
                                    new_tellPixels = [new_tellPixels; new_tells(m).PixelIdxList];
                                end
                                new_mask_tells = false(512, 672);
                                new_mask_tells(new_tellPixels) = 1;
                                mask_tells = new_mask_tells;

                                mask3 = imtranslate(mask3, directions(I, 1:2));
                                new_GFP = imtranslate(new_GFP, directions(I, 1:2));
    %                             figure();
    %                             imshow(mask_tells);
                            end
                        end


                        %For grabbing images
    %                     figure();
    %                     imshow(new_GFP);
    %                     figure();
    %                     imshow(mask2);


                        cells = regionprops(mask2, new_GFP, 'Area', 'BoundingBox', 'Image', 'PixelList', 'PixelIdxList', 'PixelValues', 'Centroid', 'WeightedCentroid', 'MeanIntensity');
                        cellsin = [cells.MeanIntensity];
                        inmean = mean(cellsin);
                        instd = std(cellsin);
                        cells = cells(cellsin >= (inmean - (2 * instd)) & cellsin <= (inmean + (2 * instd)));
                        cellsroster = [zeros(length(cells), 1) + i, zeros(length(cells), 1) + well.imageid];
                        actualCellCounts(i) = actualCellCounts(i) + length(cells);
                        actualImageCounts(i) = actualImageCounts(i) + 1;
                        if firstTime
                            firstTime = false;
                            cellMasterData = cells;
                            cellMasterRoster = cellsroster;
                        else
                            cellMasterData(end + 1 : end + length(cells)) = cells;
                            cellMasterRoster = [cellMasterRoster; cellsroster];
                        end
                        
                        actualCellCounts
                        actualImageCounts

                        %From presentation code - to explore
    %                     newGFPMatrix = zeros(length(cells), 625) - ones(length(cells), 625);
    %                     for k = 1:length(cells)
    %                         for m = 1:length(cells(k).PixelIdxList)
    %                             newGFPMatrix(k, m) = new_GFP(cells(k).PixelIdxList(m));
    %                         end
    %                     end
    %                     GFPMatrix((sum(cellCounts) + 1):(sum(cellCounts) + length(cells)), :) = newGFPMatrix;
    %                     cellLabels((sum(cellCounts) + 1):(sum(cellCounts) + length(cells))) = i;
    %                     cellCounts(i) = cellCounts(i) + length(cells);

                    else
                        j = j - 1;
                    end
                catch
                    j = j - 1;
                    disp("Error occured, trying again...");
                end
            end
        end
%         for j = 1:minCellsList(i)
%             %pick random cells
%         end
    end
end

save('dataset1.mat', 'actualCellCounts', 'actualImageCounts', 'cellMasterRoster', 'cellMasterData');




% %Visualize the distribution of cell sizes
% cellCount = 0;
% cells = zeros(200000, 1);
% for i = 1:100
%     ourSamples = sampleList(sampleList.Localization == randi(11), :);
%     selection = randi(length(ourSamples.Localization));
%     downSelect = imageRegistry(strcmp({imageRegistry.plate}, char(ourSamples.Plate(selection))));
%     well = downSelect(strcmp({downSelect.well}, char(ourSamples.Well(selection))));
%     url = [api '/webgateway/render_image/' num2str(well.imageid) '/0/0/'];
%     img = webread(url);
%     mCherry = img(:, :, 1);
% 
%     filter = fspecial('average', [30 30]);
%     test = imfilter(mCherry, filter);
%     test2 = mCherry - test;
%     mask = test2 > 3*mean(test2, 'all');
% 
%     cellAreas = regionprops(mask, 'Area');
%     spot = 1;
%     for j = (cellCount + 1):(cellCount + length(cellAreas))
%         cells(j) = cellAreas(spot).Area;
%         spot = spot + 1;
%     end
%     cellCount = cellCount + length(cellAreas);
% end
% cells = cells(cells > 0);
% cells = cells(cells > 75);
% cells = cells(cells < 625);
% figure();
% histogram(cells);








% %Temporary workflow for 12/5 presentation
% locRatio = [0, 0, 10, 0, 0, 0, 0, 6, 0, 0, 0];
% 
% %Set minimum number of cells per type
% minCells = 100;
% 
% %Set minimum number of images per type
% minImages = 10;
% 
% %Adjust search constraints
% minCellsList = locRatio * minCells;
% 
% cellCounts = zeros(1, 11);
% for i = 1:11
%     if locRatio(i) > 0
%         ourSamples = sampleList(sampleList.Localization == i, :);
%         for j = 1:minImages
%             %Collect all useful cells
%             selection = randi(length(ourSamples.Localization));
%             downSelect = imageRegistry(strcmp({imageRegistry.plate}, char(ourSamples.Plate(selection))));
%             well = downSelect(strcmp({downSelect.well}, char(ourSamples.Well(selection))));
%             url = [api '/webgateway/render_image/' num2str(well.imageid) '/0/0/'];
%             img = webread(url);
%             mCherry = img(:, :, 1);
%             GFP = img(:, :, 2);
%             
%             filter = fspecial('average', [30 30]);
%             temp = imfilter(mCherry, filter);
%             temp2 = mCherry - temp;
%             mask = temp2 > 3*mean(temp2, 'all');
%             
% %             lmask = bwconncomp(mask);
%             cells = regionprops(mask, 'Area', 'BoundingBox', 'Image', 'PixelIdxList');
%             area = [cells.Area];
%             cells = cells(area >= 75 & area <= 625);
%             dim1 = zeros(length(cells), 1);
%             dim2 = zeros(length(cells), 1);
%             for k = 1:length(cells)
%                 dim1(k) = cells(k).BoundingBox(3);
%                 dim2(k) = cells(k).BoundingBox(4);
%             end
%             cells = cells(dim1 <= 40 & dim2 <= 40);
% 
%             
%             filter2 = fspecial('average', [50 50]);
%             temp3 = imfilter(GFP, filter);
%             new_GFP = GFP - temp3;
%             
%             newGFPMatrix = zeros(length(cells), 625) - ones(length(cells), 625);
%             for k = 1:length(cells)
%                 for m = 1:length(cells(k).PixelIdxList)
%                     newGFPMatrix(k, m) = new_GFP(cells(k).PixelIdxList(m));
%                 end
%             end
%             GFPMatrix((sum(cellCounts) + 1):(sum(cellCounts) + length(cells)), :) = newGFPMatrix;
%             cellLabels((sum(cellCounts) + 1):(sum(cellCounts) + length(cells))) = i;
%             cellCounts(i) = cellCounts(i) + length(cells);
%             
%         end
%     end
% end

% %Analysis
% ranges = zeros(length(cellLabels), 1);
% for i = 1:length(ranges)
%     temp = GFPMatrix(i, GFPMatrix(i, :) >= 0);
%     ranges(i) = max(temp) - min(temp);
% end
% 
% figure(1);
% histogram(ranges)
% xlabel('Range of GFP intensities in cell');
% ylabel('Cells');
% 
% 
% lowerBound = min(ranges) + 1;
% upperBound = max(ranges) - 1;
% precision = zeros(upperBound - lowerBound, 1);
% recall = zeros(upperBound - lowerBound, 1);
% score = zeros(upperBound - lowerBound, 1);
% threshold = zeros(upperBound - lowerBound, 1);
% for j = lowerBound:upperBound
%     threshold(j - lowerBound + 1) = j;
% 
%     output = zeros(length(ranges), 1);
%     for i = 1:length(output)
%         state = ranges(i) >= threshold(j - lowerBound + 1);
%         if state
%             output(i) = 8;
%         else
%             output(i) = 3;
%         end
%     end
% 
%     outcome = confusionmat(cellLabels', output);
%     %below, consider nuclear localization to be "positive"
%     TP = outcome(2, 2);
%     TN = outcome(1, 1);
%     FP = outcome(1, 2);
%     FN = outcome(2, 1);
%     precision(j - lowerBound + 1) = outcome(2, 2) / (outcome(2, 2) + outcome(1, 2));
%     recall(j - lowerBound + 1) = outcome(2, 2) / (outcome(2, 2) + outcome(2, 1));
%     score(j - lowerBound + 1) = (TP * TN - FP * FN) / sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN));
% end
% 
% %Matthews correlation coefficient
% 
% figure(2);
% plot(threshold, score, 'LineWidth', 1.5);
% xlabel('Threshold for Classification (Pixel Intensity Range)');
% ylabel('Matthews Correlation Coefficient');
% xlim([lowerBound upperBound]);
% figure(3);
% plot(threshold, recall, 'LineWidth', 1.5);
% hold on
% plot(threshold, precision, 'LineWidth', 1.5);
% legend('Recall', 'Precision');
% ylabel('Performance Metric');
% xlabel('Threshold for Classification (Pixel Intensity Range)');
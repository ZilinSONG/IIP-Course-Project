function Task1()

function [faceModel] = getFaceModel1()  % face1 use YUV model
    img = imread('face1.1.jpg'); % read file
    ycbcr = rgb2ycbcr(img); % RGB -> YUV
    cb = ycbcr(:,:,2); % U space
    cr = ycbcr(:,:,3); % V space
    faceModel1=[max(max(cb)), min(min(cb)),max(max(cr)), min(min(cr))]; % get threshold for U and V space
    img = imread('face1.2.jpg'); % get another image
    ycbcr = rgb2ycbcr(img);
    cb = ycbcr(:,:,2);
    cr = ycbcr(:,:,3);
    faceModel2=[max(max(cb)), min(min(cb)),max(max(cr)), min(min(cr))]; % get threshold for U and V space
    faceModel_max = max(faceModel1,faceModel2); % compare two threshold 
    faceModel_min = min(faceModel1,faceModel2);
    faceModel=[faceModel_max(1),faceModel_min(2),faceModel_max(3),faceModel_min(4)]; % upload new threshold
end

function [faceModel] = getFaceModel2()  % face2 use HSV model
    img = imread('face2.1.jpg'); % read file
    hsv = rgb2hsv(img); % RGB -> HSV
    H = hsv(:,:,1); % H space
    S = hsv(:,:,2); % S space
    faceModel=[max(max(H)), min(min(H)),max(max(S)), min(min(S))]; % get threshold for H and S space
end

function binarizedFaceImage = getBinarizedFaceImage1(testImage,faceModel)
    img = imread(testImage); % read test file
    binarizedFaceImage = rgb2gray(img); % convert to grayscale image
    ycbcr = rgb2ycbcr(img); % get YUV
    cb = ycbcr(:,:,2); % get U space
    cr = ycbcr(:,:,3); % get V space
    [m, n] = size(binarizedFaceImage); % size of the grayscale image
    for m1 = 1:m
        for n1 = 1:n
            if((cb(m1,n1) <= faceModel(1) && cb(m1,n1)>= faceModel(2))&&(cr(m1,n1) <= faceModel(3) && cr(m1,n1)>= faceModel(4)))
                binarizedFaceImage(m1,n1) = 255; % set to 255 if the value within the threshold
            else
                binarizedFaceImage(m1,n1) = 0;
            end
        end
    end
end


function binarizedFaceImage = getBinarizedFaceImage2(testImage,faceModel)
    img = imread(testImage); % read test file
    binarizedFaceImage = rgb2gray(img); % convert to grayscale image
    HSV = rgb2hsv(img); % get HSV
    H = HSV(:,:,1); % get H space
    S = HSV(:,:,2); % get S space
    [m, n] = size(binarizedFaceImage); % size of the grayscale image
    for m1 = 1:m
        for n1 = 1:n
            if((H(m1,n1) <= faceModel(1) && H(m1,n1)>= faceModel(2))&&(S(m1,n1) <= faceModel(3) && S(m1,n1)>= faceModel(4)))
                binarizedFaceImage(m1,n1) = 255; % set to 255 if the value within the threshold
            else
                binarizedFaceImage(m1,n1) = 0;
            end
        end
    end
end


function FinalFaceImage=getMorphFace(binarizedFaceImage)
     %FinalFaceImage = bwmorph(binarizedFaceImage,'open');
     %FinalFaceImage = bwmorph(FinalFaceImage,'close');
     se=strel('square',4);  %Morphological structural element
     FinalFaceImage=imopen(binarizedFaceImage,se); % Morphological opening. 
     FinalFaceImage=imclose(FinalFaceImage,se);  %Morphological closing
end

function task1main()
    faceModel1 = getFaceModel1(); % get threshold
    faceModel2 = getFaceModel2();
    binarizedFaceImage1 = getBinarizedFaceImage1('face1.jpg',faceModel1); % test face1
    binarizedFaceImage2 = getBinarizedFaceImage2('face2.jpg',faceModel2); % tese face2
    FinalFaceImage1=getMorphFace(binarizedFaceImage1); % morphological operation 
    FinalFaceImage2=getMorphFace(binarizedFaceImage2); % morphological operation 
    figure('name','The first iamage');
    subplot(1,2,1);imshow(binarizedFaceImage1);
    title('Binarized Face Image');
    subplot(1,2,2);imshow(FinalFaceImage1);
    title('Final Face Image');
    figure('name','The second iamage');
    subplot(1,2,1);imshow(binarizedFaceImage2);
    title('Binarized Face Image');
    subplot(1,2,2);imshow(FinalFaceImage2);
    title('Final Face Image');
end

task1main();

end
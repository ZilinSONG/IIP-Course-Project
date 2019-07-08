function[]= Task2()

function [mean_img] = Mean_filter(img)   % Mean Filter
    m = ones(3);   % filter template with size of 3 * 3
    mf = m * 1/9;  % modify template
    img_r = img(:,:,1); % R channels
    img_g = img(:,:,2); % G channels
    img_b = img(:,:,3); % B channels
    mr = conv2(img_r, mf); % filter R channel
    mg = conv2(img_g, mf); % filter G channel
    mb = conv2(img_b, mf); % filter B channel
    mean_img = uint8(cat(3,mr,mg,mb)); %Regroup RGB image
    
end

function [Gas_img] = Gas_filter1(img)   % Gaussian Filter with delta = 1 and size 2 X 2
    gf = zeros(5); % filter template
    r = 2; %filter randius
    sum = 0;  % SUM of P(X,Y)^2
    delta = 1; % delta in guassian filter
    for row = -r:r
        for col = -r:r
            gf(row+3,col+3) = (exp(1)^(-(row^2 +col^2)/(2*delta^2)))/(2*pi*delta^2); %Calculate P(X,Y)^2
            sum = sum + gf(row+3,col+3);
        end
    end
    for row = -r:r
        for col = -r:r
            gf(row+3,col+3) = gf(row+3,col+3)/sum; %calculate pixel value
        end
    end
    img_r = img(:,:,1); % R channel
    img_g = img(:,:,2); % G channel
    img_b = img(:,:,3); % B channel
    gr = conv2(img_r, gf);
    gg = conv2(img_g, gf);
    gb = conv2(img_b, gf);
    Gas_img = uint8(cat(3,gr,gg,gb));
       
end
function [Gas_img] = Gas_filter2(img)  % Gaussian Filter with delta = 2 and size 7 X 7
    gf = zeros(7); % filter template
    r = 3; %filter randius
    sum = 0;  % SUM of P(X,Y)^2
    delta = 2; % delta in guassian filter
    for row = -r:r
        for col = -r:r
            gf(row+4,col+4) = (exp((-(row^2 +col^2)/(2*delta^2))))/(2*pi*delta^2);
            sum = sum + gf(row+4,col+4);
        end
    end
    for row = -r:r
        for col = -3:3
            gf(row+4,col+4) = gf(row+4,col+4)/sum;
        end
    end
    img_r = img(:,:,1); % R channel
    img_g = img(:,:,2);
    img_b = img(:,:,3);
    gr = conv2(img_r, gf);
    gg = conv2(img_g, gf);
    gb = conv2(img_b, gf);
    Gas_img = uint8(cat(3,gr,gg,gb)); % group RBG channels
       
end

function [image2] = med_filter(image1)  % Median Filter
    [row,col] = size(image1);
    image2 = zeros(row, col,'uint8'); % to obtain the processed image
    image_tem = zeros(row+2, col+2,'double'); % used to process the image
    for r = 1:row
        for c = 1:col
            image_tem(r+1,c+1) = image1(r,c); % assignment value
        end
    end
    % assign values to the for corners
    image_tem(1,1) = image_tem(2,2);
    image_tem(1,2+col) = image_tem(2,col);
    image_tem(2+row, 1) = image_tem(1+row,2);
    image_tem(2+row,2+col) = image_tem(row+1, col+1);
    % assign values to four sides
    for r = 1: row  % right and left side
        image_tem(r+1,1) = image_tem(r+1,2);
        image_tem(r+1,col+2) = image_tem(r+1,col+1);
    end
    for c = 1:col   % top and bottom side
        image_tem(1, c+1) = image_tem(2,c+1);
        image_tem(row+2,c+1) = image_tem(row+1, c+1);
    end
    for i = 2:row+1
        for j = 2:col+1
            t_image =image_tem(i-1:i+1,j-1:j+1); % the matrix to be processed this loop
            med1 = sort(t_image(:)); % sort median in this matrix
            med = med1(5);   % get median
            image2(i-1,j-1) = uint8(med);
        end
    end
end

function [Med_img] = Med_filter(img)   % Median Filter
    img_r = img(:,:,1); % R channel
    img_g = img(:,:,2); % G
    img_b = img(:,:,3); % B
    medr = med_filter(img_r); % fliter R channel
    medg = med_filter(img_g);
    medb = med_filter(img_b);
    Med_img = uint8(cat(3,medr,medg,medb));

end

function [image2] = ani_filter(image1)  % Anisotropic filter
        [row,col] = size(image1);
        image2 = zeros(row, col,'uint8'); % to obtain the processed image
        image_tem = zeros(row+2, col+2,'double'); % used to process the image
        for r = 1:row
            for c = 1:col
                image_tem(r+1,c+1) = image1(r,c); % assignment value
            end
        end
        % assign values to the for corners
        image_tem(1,1) = image_tem(2,2);
        image_tem(1,2+col) = image_tem(2,col);
        image_tem(2+row, 1) = image_tem(1+row,2);
        image_tem(2+row,2+col) = image_tem(row+1, col+1);
        % assign values to four sides
        for r = 1: row % right and left side
            image_tem(r+1,1) = image_tem(r+1,2);
            image_tem(r+1,col+2) = image_tem(r+1,col+1);
        end
        for c = 1:col  % top and bottom side
            image_tem(1, c+1) = image_tem(2,c+1);
            image_tem(row+2,c+1) = image_tem(row+1, c+1);
        end
        for i = 2:row+1
            for j = 2:col+1
                t_image =image_tem(i-1:i+1,j-1:j+1); % the matrix to be processed this loop
                D = max(max(t_image)) - min(min(t_image)); % the maximum possible difference
                D = double(D);
                delta_all = 1-abs((t_image-t_image(2,2)))/(D+10^(-10)); % delta materix(1-d/D)
                delta_all(2,2) = 0;  % set centre to 0
                delta = sum(delta_all);    % sum of s(p,q)
                delta_all = delta_all.*t_image; % get q*s(p,q)
                ani = sum(delta_all)/delta; % get result
                image2(i-1,j-1) = uint8(ani);
            end
        end
    end

function [Ani_img] = Ani_filter(img)  % Anisotropic filter
    img_r = img(:,:,1); % R channel
    img_g = img(:,:,2);
    img_b = img(:,:,3);
    anir = ani_filter(img_r); % filtered R channel
    anig = ani_filter(img_g);
    anib = ani_filter(img_b);
    Ani_img = uint8(cat(3,anir,anig,anib));
end

function[image2] = bil_filter(image1)  % Bilateral filter
    rad = 2;  % radiua
    sigmas = 1; % sigmas for space
    sigmar = 10; % rage sigmas
   [row,col] = size(image1);
   image2 = zeros(row, col,'uint8');
   image_tem = zeros(row+4, col+4,'double'); % used to process
   for r = 1:row
       for c = 1:col
           image_tem(r+2,c+2) = image1(r,c); % assignment value
       end
   end

   [x,y] = meshgrid(-2:2,-2:2);
   filter1 = exp(-(x.^2 + y.^2)/(2*sigmas^2))/(2*pi*sigmas^2); % range filter 
   for i = rad+1:row+rad
       for j = rad+1:col+rad
           t_image =image_tem(i-rad:i+rad,j-rad:j+rad); % the matrix processed this loop
           filter2 = exp(-(t_image-image_tem(i,j)).^2/(2*sigmar^2))/(2*pi*sigmar^2); % space filter
           fin_filter = filter1*filter2; % bilateral filter 
           res = image_tem(i-rad:i+2,j-2:j+2).*fin_filter;
           fin = sum(sum(res))/sum(sum(fin_filter)); % calculate and get value
           image2(i-1,j-1) = uint8(fin); % assign
        end
    end
   
    
end

function [bil_img] = Bil_filter(img)    % Bilateral filter
    img_r = img(:,:,1);
    img_g = img(:,:,2);
    img_b = img(:,:,3);
    bilr = bil_filter(img_r);
    bilg = bil_filter(img_g);
    bilb = bil_filter(img_b);
    bil_img = uint8(cat(3,bilr,bilg,bilb));
    imshow(bil_img)        
end

function [spn_image] = SP_noise(img)   % add salt & pepper filter
    img_r = double(img(:,:,1));
    img_g = double(img(:,:,2));
    img_b = double(img(:,:,3));
    [row, col] = size(img_r);
    for r = 1:row
        for c = 1:col
            if rand()<0.1 % if nosie
                if rand()<0.5 % white of black
                    img_r(r,c) = 0;
                    img_g(r,c) = 0;
                    img_b(r,c) = 0;
                else
                    img_r(r,c) = 255;
                    img_g(r,c) = 255;
                    img_b(r,c) = 255;
                end
            end
        end
    end
    spn_image = uint8(cat(3,img_r,img_g,img_b));
end

function [guan_image] = Gua_nosie(img)  % add gaussian noise
    img_r = img(:,:,1);
    img_g = img(:,:,2);
    img_b = img(:,:,3);
    s = 20; % sigma
    sizergb = size(img_r);
    gas_n = randn(sizergb)*s;
    gas_r = double(img_r)+gas_n;
    gas_g = double(img_g)+gas_n;
    gas_b = double(img_b)+gas_n;
    guan_image = uint8(cat(3,gas_r,gas_g,gas_b));
end

function task2main()
    img = imread('lena.jpg');
    dis_image1 = Gua_nosie(img);   % distorted image with Gaussian noise 考=20 
    dis_image2 = SP_noise(img);    % distorted image with 10% of salt & pepper noise
    dis_image3 = Gas_filter2(img); % distorted image with a 7℅7 Gaussian filter with 考 = 2
    figure('name','Distorted Image');
    subplot(2,2,1);imshow(img);
    title('Original Image');
    subplot(2,2,2);imshow(dis_image1);
    title('Gaussian noise with 考 =20');
    subplot(2,2,3);imshow(dis_image2);
    title('10% of salt & pepper noise');
    subplot(2,2,4);imshow(dis_image3);
    title('7℅7 Gaussian filter with 考 = 2');
    figure('name','Image with Gaussian noise');
    subplot(2,3,1);imshow(dis_image1);
    title('Image with Gaussian noise');
    subplot(2,3,2);imshow(Mean_filter(dis_image1));
    title('Mean filter');
    subplot(2,3,3);imshow(Gas_filter1(dis_image1));
    title('Gaussian filter');
    subplot(2,3,4);imshow(Med_filter(dis_image1));
    title('Median filter');
    subplot(2,3,5);imshow(Ani_filter(dis_image1)); 
    title('Anisotropic filter');
    subplot(2,3,6);imshow(Bil_filter(dis_image1));
    title('Bilateral filter ');
    figure('name','Image with S&P noise');
    subplot(2,3,1);imshow(dis_image2);
    title('Image with S&P noise');
    subplot(2,3,2);imshow(Mean_filter(dis_image2));
    title('Mean filter');
    subplot(2,3,3);imshow(Gas_filter1(dis_image2));
    title('Gaussian filter');
    subplot(2,3,4);imshow(Med_filter(dis_image2));
    title('Median filter');
    subplot(2,3,5);imshow(Ani_filter(dis_image2));
    title('Anisotropic filter');
    subplot(2,3,6);imshow(Bil_filter(dis_image2));
    title('Bilateral filter ');
    figure('name','Image with Gaussian filter');
    subplot(2,3,1);imshow(dis_image3);
    title('Image with Gaussian filter');
    subplot(2,3,2);imshow(Mean_filter(dis_image3));
    title('Mean filter');
    subplot(2,3,3);imshow(Gas_filter1(dis_image3));
    title('Gaussian filter');
    subplot(2,3,4);imshow(Med_filter(dis_image3));
    title('Median filter');
    subplot(2,3,5);imshow(Ani_filter(dis_image3));
    title('Anisotropic filter');
    subplot(2,3,6);imshow(Bil_filter(dis_image3));
    title('Bilateral filter ');
end

task2main();
end
clear all, clc, close all

%% 1.1

f = 5.405e9;
c = 3e8;
subs = 3;
infosFileMaster = h5info("data_and_functions/20190704.h5");
realPartImage = h5read(infosFileMaster.Filename, "/i_VV");
realPartImage = realPartImage(1:subs:end, 1:subs:end);
imagPartImage = h5read(infosFileMaster.Filename, "/q_VV");
imagPartImage = imagPartImage(1:subs:end, 1:subs:end);
distancesMaster = h5read(infosFileMaster.Filename, "/topoDistance");
distancesMaster = distancesMaster(1:subs:end, 1:subs:end);
master = single(realPartImage)+1j.*single(imagPartImage);
[rows, columns] = size(master);

absmaster = abs(master);
phasemaster = angle(master);
auxvar = sort(reshape(abs(master), [rows*columns, 1]));
treshold = auxvar(round(0.98*length(auxvar)));
absmaster(absmaster > treshold) = treshold;
clear auxvar; % now we have clipped the absolute values

load('data_and_functions/geocoding_infos.mat'); %here we load the geocoding infos
geocodingInfos = GecSubs(geocodingInfos, [subs subs]); %since they never change, we load them at the beginning and keep use them read only

filtersize = 9;
filter = ones(filtersize);
filter = filter/sum(sum(filter)); % we normalize the gaussian filter


master = absmaster.*exp(1j*phasemaster); %here we get back the complex matrix after clipping the peaks of the abs value

figure(1)
subplot(1, 2, 1)
imagesc(abs(master))
title('Abs. val. of master image')
subplot(1, 2, 2)
imagesc(angle(master))
title('Phase of master image')

figure(2)
imagesc(conv2(abs(master), filter))
title('Multilooked abs. val of master image')

masterold = master;

%here we geocode the multilooked image
fl = floor(length(filter)/2);
geocodedImage = geocodeRadarImage(conv2(abs(master(1+fl:end-fl, 1+fl:end-fl)), filter), geocodingInfos);

figure(3)
grid2image(geocodedImage, geocodingInfos.xref);
title('Geocoded master image')
%here we write the geotiff for the point 1.1
geotiffwrite('1-1geotiff.tif', geocodedImage, geocodingInfos.xref) 
clear geocodedImage
%% 1.2

% subs = 3;
infosFileMaster = h5info("data_and_functions/20190716.h5"); %%load the second image
realPartImage = h5read(infosFileMaster.Filename, "/i_VV");
realPartImage = realPartImage(1:subs:end, 1:subs:end);
imagPartImage = h5read(infosFileMaster.Filename, "/q_VV");
imagPartImage = imagPartImage(1:subs:end, 1:subs:end);
slave = single(realPartImage)+1j.*single(imagPartImage);
[rows, columns] = size(slave);

absmaster = abs(slave);
phasemaster = angle(slave);
auxvar = sort(reshape(abs(slave), [rows*columns, 1]));
treshold = auxvar(round(0.98*length(auxvar)));
absmaster(absmaster > treshold) = treshold;
clear auxvar; % now we have clipped the absolute values

slave = absmaster.*(exp(1j*phasemaster)); % here we get back the complex matrix after clipping the peaks of the abs value
clear absmaster phasemaster

slavenew = slave;

figure(4)
subplot(1, 2, 1)
imagesc(angle(masterold.*conj(slavenew))); % interferometric phase as it is, with fast fringes
title('Interferometric phase')
subplot(1, 2, 2)
imagesc(angle(conv2(masterold.*conj(slavenew), filter)));
title('Interferometric phase (multilooked)')

figure(5)
subplot(1, 2, 1)
imagesc(conv2(abs(master), filter))
title('Multilooked abs. val of master image')
subplot(1, 2, 2)
imagesc(angle(conv2(masterold.*conj(slavenew), filter)));
title('Interferometric phase (multilooked)')

infosFileMaster = h5info("data_and_functions/20190704.h5");
distancesMaster1 = h5read(infosFileMaster.Filename, "/topoDistance"); %read topoDistance of master image
distancesMaster1 = distancesMaster1(1:subs:end, 1:subs:end); %%subsample according to subs variable
masterold = masterold.*exp(1j.*(4*pi*f/c)*distancesMaster1); %%compensate for topography phase changes

infosFileMaster = h5info("data_and_functions/20190716.h5");
distancesMaster2 = h5read(infosFileMaster.Filename, "/topoDistance"); %read topoDistance of slave image
distancesMaster2 = distancesMaster2(1:subs:end, 1:subs:end);  %%subsample according to subs variable
slavenew = slavenew.*exp(1j.*(4*pi*f/c)*distancesMaster2); %%compensate for topography phase changes
% master old and slave new corrected with the distances

figure(6)
subplot(1, 2, 1)
imagesc(angle(masterold.*conj(slavenew)));
title('Int. phase (after compensation) - 1.2')
subplot(1, 2, 2)
imagesc(angle(conv2(masterold.*conj(slavenew), filter)));
title('Int. phase (after comp, multilooked) - 1.3')

% 1.2 QGIS section
interf = angle(masterold.*conj(slavenew));
geocodedImage = geocodeRadarImage(interf, geocodingInfos);
geotiffwrite("1-2interferometricPhaseDirty.tif", geocodedImage, geocodingInfos.xref)
clear interf
clear geocodedImage
%% 1.3
% 1.3 QGIS section
interf = angle(conv2(masterold.*conj(slavenew), filter)); 
interf = interf(1+fl:end-fl, 1+fl:end-fl);
geocodedImage = geocodeRadarImage(interf, geocodingInfos);
geotiffwrite("1-3interferometricPhaseClean.tif", geocodedImage, geocodingInfos.xref)
clear interf
clear geocodedImage
%% 2
interferogram = masterold.*conj(slavenew);

slavenewfilt = conv2(abs(slavenew(1+floor(length(filter)/2):end-floor(length(filter)/2), 1+floor(length(filter)/2):end-floor(length(filter)/2))), filter).*exp(i.*(angle(slavenew)));
masteroldfilt = conv2(abs(masterold(1+floor(length(filter)/2):end-floor(length(filter)/2), 1+floor(length(filter)/2):end-floor(length(filter)/2))), filter).*exp(i.*(angle(masterold)));

coherenceMap = masteroldfilt.*conj(slavenewfilt)./((masteroldfilt.*conj(masteroldfilt).*slavenewfilt.*conj(slavenewfilt)).^0.5);
coherenceMap = conv2(coherenceMap(1+fl:end-fl, 1+fl:end-fl), filter);
coherenceSubsampled = coherenceMap(1:20:end, 1:5:end);
unwrappedInterferogram = unwrap_IRLS(double(angle(coherenceSubsampled)),coherenceSubsampled, 100, [], 1);

%2.1
figure(7)
subplot(1, 2, 1)
imagesc(abs(coherenceMap));
title('Abs. val. of coh. map (multilooked)')
subplot(1, 2, 2)
imagesc(angle(coherenceMap));
title('Phase of coh. map (multilooked)')

% 2.1 QGIS section
geocodedImage = geocodeRadarImage(abs(coherenceMap), geocodingInfos);
geotiffwrite("2-1AbsCohMap.tif", geocodedImage, geocodingInfos.xref)
clear geocodedImage

%2.2
figure(8)
imagesc(unwrappedInterferogram);
title('Unwrapped Interferogram')

% 2.2 QGIS section
[Xq, Yq] = meshgrid(0:(size(unwrappedInterferogram, 2)/size(slavenewfilt,2)):size(unwrappedInterferogram, 2), 0:(size(unwrappedInterferogram, 1)/size(slavenewfilt,1)):size(unwrappedInterferogram, 1));
unint = interp2(unwrappedInterferogram, Xq, Yq);
unint = unint( 1:end-1, 1:end-1);
geocodedImage = geocodeRadarImage(unint, geocodingInfos);
geotiffwrite("2-2UnwrappedInterf.tif", geocodedImage, geocodingInfos.xref)
clear unint
clear geocodedImage

%% 3
gradientx = conv2(unwrappedInterferogram, [-1 0 1; -1 0 1; -1 0 1]);
gradienty = conv2(unwrappedInterferogram, [-1 -1 -1; 0 0 0 ; 1 1 1]);
gradient = sqrt(gradientx.^2 + gradienty.^2);

figure(9)
subplot(1, 2, 1)
imagesc(gradient);
title('Gradient of unwr. interf.')
subplot(1, 2, 2)
imagesc(conv2(movmax(gradient, [1 1]), filter));
title('Gradient of unwr. interf. (enhanced look)')

% 3.1 QGIS section
[Xq, Yq] = meshgrid(0:(size(gradient, 2)/size(slavenewfilt,2)):size(gradient, 2), 0:(size(gradient, 1)/size(slavenewfilt,1)):size(gradient, 1));
geoGrad = interp2(gradient, Xq, Yq);
geoGrad = geoGrad( 1:end-1, 1:end-1);
geocodedImage = geocodeRadarImage(geoGrad, geocodingInfos);
geotiffwrite("3-1Gradients.tif", geocodedImage, geocodingInfos.xref)
clear geoGrad
clear geocodedImage

figure(10)
improfile(unwrappedInterferogram*c/4/f/pi, [150 300], [100 350])
title('Profile of the ground movement')

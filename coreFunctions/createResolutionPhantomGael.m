function phantom = createResolutionPhantomGael(N, fov, diameters)

%N = [200 200 1];
%fov = [0.08 0.08 0.0001]; %fov.size;
%diameters = [0.012, 0.006, 0.012, 0.006]; % [m] diameter of the circles 

    fovCm = fov.* 100;
    diameters = 100.*diameters;

    kx = linspace(-fovCm(1)/2, fovCm(1)/2, N(1)); % cm -> m
    ky = linspace(-fovCm(2)/2, fovCm(2)/2, N(2));
    [x,y] = meshgrid(kx,ky);

    fovCm = fovCm * 0.9;
    phantom = zeros(size(x));

    signX = [-1, -1, 1, 1];
    signY = [-1, 1, -1, 1];
    circels = [0,0];

    for l=1:length(diameters)
        numberOfCircelsX = floor((fovCm(1)./2.0 + (1-0.6) .* diameters)./(2.*diameters(l)));
        numberOfCircelsY = floor((fovCm(2)./2.0 + (1-0.6) .* diameters)./(2.*diameters(l)));
        %numberOfCircelsY = floor(fovCm(2)./2.0./(2.*diameters(l)+0.6));

        for kx=0:numberOfCircelsX-1
            for ky=0:numberOfCircelsY-1

                circels(1) = signX(l)* 2 * (kx+0.6) * diameters(l);
                circels(2) = signY(l)* 2 * (ky+0.6) * diameters(l);
                
                dist = ((x+circels(1))./diameters(l).*2).^2 + ((y+circels(2))./diameters(l).*2).^2;
                
                phantom(dist <= 1) = 1;
            end
        end
    end
    phantom = phantom';    
end
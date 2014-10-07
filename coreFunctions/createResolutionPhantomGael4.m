function phantom = createResolutionPhantomGael4(N, thickness)
% Create a single square almost in the middle
% with a thickness of 'thickness' pixel

% N = [93 91];
% thickness = 2;

    phantom = zeros(N(1),N(2));
    i=1;
    for i=floor((N(1)-thickness)/2):floor((N(1)+thickness)/2)-1
        for j=floor((N(2)-thickness)/2):floor((N(2)+thickness)/2)-1
            phantom(i,j) = 1;
        end
    end
end
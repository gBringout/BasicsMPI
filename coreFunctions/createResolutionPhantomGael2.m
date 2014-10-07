function phantom = createResolutionPhantomGael2(N, thickness)

% N = [93 91];
% thickness = 3;

    phantom = zeros(N(1),N(2));
    i=1;
    while i<N(2)
        k=1;
        l=1;
        % Make the top of the cross
        while l<=thickness % on l line
            while k<N(1) %until the end
                for j=1:thickness % for j pixel
                    phantom(i,k) = 0;
                    k=k+1;
                end
                for j=1:thickness
                    phantom(i,k) = 1;
                    k=k+1;
                end
                for j=1:thickness
                    phantom(i,k) = 0;
                    k=k+1;
                end
                for j=1:thickness
                    phantom(i,k) = 0;
                    k=k+1;
                end
            end
            k=1;
            i=i+1;
            l=l+1;
        end
        % Middle of the coil
        k=1;
        l=1;
        % Make the top of the cross
        while l<=thickness % on l line
            while k<N(1) %until the end
                for j=1:thickness % for j pixel
                    phantom(i,k) = 1;
                    k=k+1;
                end
                for j=1:thickness
                    phantom(i,k) = 0;
                    k=k+1;
                end
                for j=1:thickness
                    phantom(i,k) = 1;
                    k=k+1;
                end
                for j=1:thickness
                    phantom(i,k) = 0;
                    k=k+1;
                end
            end
            k=1;
            i=i+1;
            l=l+1;
        end
        % Bottom of the cross
        k=1;
        l=1;
        % Make the top of the cross
        while l<=thickness % on l line
            while k<N(1) %until the end
                for j=1:thickness % for j pixel
                    phantom(i,k) = 0;
                    k=k+1;
                end
                for j=1:thickness
                    phantom(i,k) = 1;
                    k=k+1;
                end
                for j=1:thickness
                    phantom(i,k) = 0;
                    k=k+1;
                end
                for j=1:thickness
                    phantom(i,k) = 0;
                    k=k+1;
                end
            end
            k=1;
            i=i+1;
            l=l+1;
        end
        % That is an empty line
        k=1;
        l=1;
        % Make the top of the cross
        while l<=thickness % on l line
            while k<N(1) %until the end
                for j=1:thickness % for j pixel
                    phantom(i,k) = 0;
                    k=k+1;
                end
            end
            k=1;
            i=i+1;
            l=l+1;
        end
    end
    phantom = phantom(1:N(),1:N(2));
    %imagesc(phantom)  
end
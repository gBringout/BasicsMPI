function phantom = createResolutionPhantomGael(N)

% N = [93 91];
% thickness = 3;

    thickness = 5;
    phantom = zeros(N(1),N(2));
    i=1;
    while i<N(1)
        k=1;
        l=1;
        % Make the top of the cross
        while l<=thickness+1 % on l line
            while k<=N(2) % on k coloumn
                if mod(l,thickness+1)==0
                    phantom(i,k) = 0;
                    phantom(i,k+1) = 0;
                    phantom(i,k+2) = 0;
                    phantom(i,k+3) = 0;
                    phantom(i,k+4) = 0;
                    k=k+6;
                elseif mod(l,thickness+1)==1
                    phantom(i,k) = 0;
                    phantom(i,k+1) = 0;
                    phantom(i,k+2) = 1;
                    phantom(i,k+3) = 0;
                    phantom(i,k+4) = 0;
                    k=k+6;
                elseif mod(l,thickness+1)==2
                    phantom(i,k) = 0;
                    phantom(i,k+1) = 1;
                    phantom(i,k+2) = 1;
                    phantom(i,k+3) = 1;
                    phantom(i,k+4) = 0;
                    k=k+6;
                elseif mod(l,thickness+1)==3
                    phantom(i,k) = 1;
                    phantom(i,k+1) = 1;
                    phantom(i,k+2) = 0;
                    phantom(i,k+3) = 1;
                    phantom(i,k+4) = 1;
                    k=k+6;
                elseif mod(l,thickness+1)==4
                    phantom(i,k) = 0;
                    phantom(i,k+1) = 1;
                    phantom(i,k+2) = 1;
                    phantom(i,k+3) = 1;
                    phantom(i,k+4) = 0;
                    k=k+6;
                elseif mod(l,thickness+1)==5
                    phantom(i,k) = 0;
                    phantom(i,k+1) = 0;
                    phantom(i,k+2) = 1;
                    phantom(i,k+3) = 0;
                    phantom(i,k+4) = 0;
                    k=k+6;
                elseif mod(l,thickness+1)==6
                    phantom(i,k) = 0;
                    phantom(i,k+1) = 0;
                    phantom(i,k+2) = 0;
                    phantom(i,k+3) = 0;
                    phantom(i,k+4) = 0;
                    k=k+6;
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
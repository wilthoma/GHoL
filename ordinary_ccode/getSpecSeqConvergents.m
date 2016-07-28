function nn = getSpecSeqConvergents(D, DD, nev1, nev2, nev3, dims)
% outputs the dimensions of the convergents of the spectral sequence
% D, DD: differentials 1->2, 3->1
% nev1,2,3 : nr non-univalent vertices in spaces 1,2,3

% find masks
MAXFIL = max([nev1;nev2;nev3]); % max filtration level

masks = cell(MAXFIL+2,3);
for j=0:MAXFIL+1
   masks{j+1,1} = (nev1==j);
   masks{j+1,2} = (nev2==j);
   masks{j+1,3} = (nev3==j);
end
nn=[];

% find kernel and image dimensions at E2 level.
% (other levels ignored so far)
%masks
%D
%DD
for j=0:MAXFIL
    % find kernel dimension (within j-th graded component)
    if isempty(D)
        kerdim=sum(masks{j+1,1});
    else               
        if j>0
            theD1 = full(D(masks{j+1,2} | masks{j,2}, masks{j+1,1}));
            theD2 = full(D(masks{j+1,2} | masks{j,2}, masks{j,1}));
            n1=size(theD1,2);
            kerdim = n1-rank([theD1,theD2])+rank(theD2);
        else
            theD1 = full(D(masks{j+1,2} , masks{j+1,1}));
            n1=size(theD1,2);
            kerdim = n1-rank(theD1);
        end
    end
      
    % find the image dimension
    if isempty(DD)
        imdim=0;
    else        
        if j>=0
            theDD1 = full(DD(masks{j+2,1}, masks{j+2,3}));
            theDD2 = full(DD(masks{j+1,1}, masks{j+1,3}));
            theDD3 = full(DD(masks{j+1,1}, masks{j+2,3}));
            if isempty(theDD1)
                theker1 = eye(size(theDD3,2));
            else
                theker1 = null(theDD1);
            end

            imdim = rank([theDD3*theker1, theDD2]);
        %else
            %theDD2 = full(DD(masks{j+1,1}, masks{j+1,3}));
            %imdim = rank(theDD2);
        end
    end
    
    nn(j+1) = kerdim - imdim;
end
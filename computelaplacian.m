function L = computelaplacian(A, options)
    % compute laplacian matrix for directed graph
    if strcmpi(options, 'directed')
        d = size(A,1);
        dout = sum(A,2);
        P = A./repmat(dout,1,d);
        P(isnan(P))=1e-6;
        [~,pi] = perron(P, 'left');
        Pi = diag(pi);
        L = Pi - (Pi*P+P'*Pi)/2;
    % compute laplacian matrix for undirected graph
    else
        D = diag(sum(A,2));
        L = D - A;
    end
end